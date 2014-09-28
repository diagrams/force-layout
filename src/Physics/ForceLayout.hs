{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TemplateHaskell     #-}
{-# LANGUAGE TypeFamilies        #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  Physics.ForceLayout
-- Copyright   :  (c) 2011 Brent Yorgey
-- License     :  BSD-style (see LICENSE)
-- Maintainer  :  byorgey@cis.upenn.edu
--
-- A simple, Haskell-native simulator for doing force-directed layout,
-- /e.g./ of trees or graphs.
--
-- To use, just create an 'Ensemble' like so:
--
-- > import           Physics.ForceLayout
-- > import qualified Data.Map              as M
-- > import           Data.AffineSpace.Point
-- > import           Data.Default (def)
-- >
-- > e :: Ensemble (Double, Double)
-- > e = Ensemble [ (edges,    hookeForce 0.05 4)
-- >              , (allPairs, coulombForce 1)
-- >              ]
-- >              particleMap
-- >   where edges       = [(1,2), (2,3), (2,5), (3,5), (3,4), (4,5)]
-- >         allPairs    = [(x,y) | x <- [1..4], y <- [x+1..5]]
-- >         particleMap = M.fromList . zip [1..]
-- >                     . map (initParticle . P)
-- >                     $ [ (2.0, 3.1), (6.3, 7.2)
-- >                       , (0.3, 4.2), (1.6, -1.1)
-- >                       , (4.8, 2.9)
-- >                       ]
--
-- Then run a simulation using either 'simulate' (to get the list of
-- all intermediate states) or 'forceLayout' (to get only the ending
-- state):
--
-- > e' :: Ensemble (Double, Double)
-- > e' = forceLayout (def & damping     .~ 0.8
-- >                       & energyLimit .~ Just 0.001
-- >                       & stepLimit   .~ Nothing
-- >                  )
-- >                  e
--
-- See the diagrams-contrib package
-- (<http://github.com/diagrams/diagrams-contrib/>) for more
-- examples.
-----------------------------------------------------------------------------

module Physics.ForceLayout
       ( -- * Data structures

         Particle(..), pos, vel, force, fixed
       , initParticle

       , PID
       , Edge
       , Ensemble(..), interForces, particleForces, particles

         -- * Pre-defined interForces

       , hookeForce
       , coulombForce
       , distForce
       , gravitate

         -- * Running simulations

       , ForceLayoutOpts(..)
       , damping, energyLimit, stepLimit
       , simulate
       , forceLayout

         -- * Internals

       , ensembleStep
       , particleStep
       , recalcForces
       , kineticEnergy

       ) where

import           Data.Default.Class
-- import           Data.Foldable          (foldMap)
import qualified Data.Map               as M
-- import           Data.Monoid

import           Control.Lens
import           Control.Monad

import           Linear.Affine
import           Linear.Vector
import           Linear.Metric

------------------------------------------------------------
--  Particles
------------------------------------------------------------

-- | A particle has a current position, current velocity, and current
--   force acting on it.
data Particle v n = Particle { _pos   :: Point v n
                             , _vel   :: v n
                             , _force :: v n
                             , _fixed :: Bool
                             }
  deriving Eq

makeLenses ''Particle

-- | Create an initial particle at rest at a particular location.
initParticle :: (Additive v, Num n) => Point v n -> Particle v n
initParticle p = Particle p zero zero False

-- resetParticle :: (Additive v, Num n) => Particle v n -> Particle v n
-- resetParticle = initParticle . view pos

------------------------------------------------------------
--  Ensembles
------------------------------------------------------------

-- | Used to uniquely identify particles.
type PID = Int

-- | An edge is a pair of particle identifiers.
type Edge = (PID, PID)

-- | An @Ensemble@ is a physical configuration of particles.  It
--   consists of a mapping from particle IDs (unique integers) to
--   particles, and a list of interForces that are operative.  Each force
--   has a list of edges to which it applies, and is represented by a
--   function giving the force between any two points.
data Ensemble v n = Ensemble
  { _interForces :: [([Edge], Point v n -> Point v n -> v n)]
       -- ^ Forces acting between particles, towards the first. The second 
       --   particle receives and equal and opposiite force. This can be used 
       --   to simulate springs and repultion between particles. See 
       --   'hookeForce' and 'coulombForce'.

  , _particleForces :: [([PID], Point v n -> v n)]
       -- ^ Force acting on a particle based on its position. This can used to 
       --   get a particle to gravitate towards or away from something constant.

  , _particles :: M.Map PID (Particle v n)
       -- ^ Map from @PID@ to correnponding particles.
  }

makeLenses ''Ensemble

------------------------------------------------------------
--  Simulation internals
------------------------------------------------------------

-- | Simulate one time step for an entire ensemble, with the given
--   damping factor.
ensembleStep :: (Additive v, Num n) => n -> Ensemble v n -> Ensemble v n
ensembleStep d = over (particles . mapped) (particleStep d) . recalcForces

-- | Simulate one time step for a particle (assuming the force acting
--   on it has already been computed), with the given damping factor.
particleStep :: (Additive v, Num n) => n -> Particle v n -> Particle v n
particleStep d p = p &~ do
  v' <- vel <.= d *^ (p^.vel ^+^ p^.force)
  unless (p ^. fixed) $ pos %= (.+^ v')

-- particleStep d = stepPos . stepVel
--   where stepVel p = vel .~ (d *^ (p^.vel ^+^ p^.force)) $ p
--         stepPos p 
--          | p ^. fixed = p
--          | otherwise  = pos %~ (.+^ p^.vel) $ p

-- | Recalculate all the interForces acting in the next time step of an
--   ensemble.
recalcForces :: (Additive v, Num n) => Ensemble v n -> Ensemble v n
recalcForces = addParticleForces . addInterForces . zeroForces
  where
    zeroForces = particles . mapped . force .~ zero

    addInterForces (Ensemble ifs pfs ps) = Ensemble ifs pfs ps'
      where
        -- apply all the particle adjustments by folding over the list of 
        -- adjustments
        ps' = foldl (.) id (concatMap mkAdjusters ifs) ps

        -- makes a list of functions that adjust the particle map for the given 
        -- set of inter-particle forces
        mkAdjusters (edges,f) = map (applyForce f) edges

        -- apply a force to two particles that make an edge
        applyForce f (i1, i2) m
          = case (M.lookup i1 m, M.lookup i2 m) of
              (Just p1, Just p2) ->
                let f' = f (p1^.pos) (p2^.pos)
                in  ( M.adjust (force ^+^~ f') i1
                    -- . M.adjust (force %~ (^-^ f')) i2
                    . M.adjust (force ^-^~ f') i2
                    ) m
              _    -> m


    -- similar to addInterForces but each function only needs a single particle
    addParticleForces (Ensemble ifs pfs ps) = Ensemble ifs pfs ps'
      where
        ps' = foldl (.) id (concatMap mkAdjuster pfs) ps

        mkAdjuster (pIds, f) = map (applyForce f) pIds

        applyForce f = M.adjust (\p -> p & force ^+^~ f (p^.pos))

(^+^~) :: (Additive v, Num n) => ASetter s t (v n) (v n) -> v n -> s -> t
l ^+^~ x = over l (^+^ x)

(^-^~) :: (Additive v, Num n) => ASetter s t (v n) (v n) -> v n -> s -> t
l ^-^~ x = over l (^-^ x)


-- | Compute the total kinetic energy of an ensemble.
kineticEnergy :: (Metric v, Num n) => Ensemble v n -> n
kineticEnergy = sumOf (particles . folded . vel . to quadrance)

------------------------------------------------------------
--  Simulation
------------------------------------------------------------

-- | Options for customizing a simulation.
data ForceLayoutOpts n =
  FLOpts
  { _damping     :: n         -- ^ Damping factor to be
                              --   applied at each step.
                              --   Should be between 0 and 1.
                              --   The default is 0.8.
  , _energyLimit :: Maybe n   -- ^ Kinetic energy below which
                              --   simulation should stop.
                              --   If @Nothing@, pay no
                              --   attention to kinetic
                              --   energy.  The default is
                              --   @Just 0.001@.
  , _stepLimit   :: Maybe Int -- ^ Maximum number of
                              --   simulation steps.  If
                              --   @Nothing@, pay no
                              --   attention to the number of
                              --   steps.  The default is
                              --   @Just 1000@.
  }
  deriving Show

makeLenses ''ForceLayoutOpts

instance Fractional n => Default (ForceLayoutOpts n) where
  def = FLOpts
        { _damping     = 0.8
        , _energyLimit = Just 0.001
        , _stepLimit   = Just 1000
        }

-- | Simulate a starting ensemble according to the given options,
--   producing a list of all the intermediate ensembles.  Useful for,
--   /e.g./, making an animation.  Note that the resulting list could
--   be infinite, if a 'stepLimit' is not specified and either the
--   kinetic energy never falls below the specified threshold, or no
--   energy threshold is specified.
simulate :: (Metric v, Num n, Ord n)
         => ForceLayoutOpts n -> Ensemble v n -> [Ensemble v n]
simulate opts e
  = (e:)
  . takeWhile (maybe (const True) (<) (opts ^. energyLimit) . kineticEnergy)
  . maybe id take (opts ^. stepLimit)
  . drop 1
  . iterate (ensembleStep (opts ^. damping))
  $ e

-- | Run a simluation from a starting ensemble, yielding either the
--   first ensemble to have kinetic energy below the 'energyLimit' (if
--   given), or the ensemble that results after a number of steps
--   equal to the 'stepLimit' (if given), whichever comes first.
--   Otherwise @forceLayout@ will not terminate.
forceLayout :: (Metric v, Num n, Ord n)
            => ForceLayoutOpts n -> Ensemble v n -> Ensemble v n
forceLayout opts = last . simulate opts

------------------------------------------------------------
--  Standard interForces
------------------------------------------------------------

-- | @distForce f p1 p2@ computes the force between two points as a
--   multiple of the unit vector in the direction from @p1@ to @p2@,
--   given a function @f@ which computes the force's magnitude as a
--   function of the distance between the points.
distForce :: (Metric v, Floating n) => (n -> n) -> Point v n -> Point v n -> v n
distForce f p1 p2 = withLength (f (distance p1 p2)) (p2 .-. p1)
  where withLength s v = s *^ signorm v

-- | @hookeForce k l p1 p2@ computes the force on @p1@, assuming that
--   @p1@ and @p2@ are connected by a spring with equilibrium length @l@
--   and spring constant @k@.
hookeForce :: (Metric v, Floating n) => n -> n -> Point v n -> Point v n -> v n
hookeForce k l = distForce (\d -> k * (d - l))

-- | @coulombForce k@ computes the electrostatic repulsive force
--   between two charged particles, with constant of proportionality
--   @k@.
coulombForce :: (Metric v, Floating n) => n -> Point v n -> Point v n -> v n
coulombForce k = distForce (\d -> -k/(d*d))

-- | @gravitate k p1 p2@ computes the force of pulling @p2@ towards @p1@ with
--   spring constant. This can be usefull for hinting a particle to a
--   particular location without freezing it, or to pull all particles towards
--   the \"center\" so that particles don't drift off.
gravitate :: (Metric v, Floating n) => n -> Point v n -> Point v n -> v n
gravitate k = distForce (* negate k)

