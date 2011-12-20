{-# LANGUAGE DeriveDataTypeable
           , NoMonomorphismRestriction 
           , TemplateHaskell
           , TypeOperators
           , FlexibleContexts
           , ScopedTypeVariables
  #-}
module Physics.ForceLayout
       ( -- * Data structures
         
         Particle, pos, vel, force
       , initParticle
         
       , Edge
       , Ensemble, forces, particles
                           
         -- * Pre-defined forces
                           
       , hookeForce
       , coulombForce
       , distForce
                           
         -- * Running simulations
                           
       , ForceLayoutOpts(..)
       , simulate
       , forceLayout
         
         -- * Internals
         
       , ensembleStep
       , particleStep
       , recalcForces
       , kineticEnergy
         
       ) where

import           Control.Monad
import           Control.Newtype
import           Data.AffineSpace
import           Data.AffineSpace.Point
import           Data.Foldable                   (foldMap)
import qualified Data.Foldable     as F
import           Data.Label                      (mkLabels)
import qualified Data.Label as L
import qualified Data.Map          as M
import           Data.Maybe
import           Data.Monoid
import           Data.VectorSpace         hiding (Sum)
import           Prelude

------------------------------------------------------------
--  Particles
------------------------------------------------------------

-- | A particle has a current position, current velocity, and current
--   force acting on it.
data Particle v = Particle { _pos   :: Point v
                           , _vel   :: v
                           , _force :: v
                           }
  deriving (Eq, Show)
           
mkLabels [''Particle]

-- | Create an initial particle at rest at a particular location.
initParticle :: AdditiveGroup v => Point v -> Particle v
initParticle p = Particle p zeroV zeroV

------------------------------------------------------------
--  Ensembles
------------------------------------------------------------

-- | An edge is a pair of particle identifiers.
type Edge = (Int, Int)

-- | An @Ensemble@ is a physical configuration of particles.  It
--   consists of a mapping from particle IDs (unique integers) to
--   particles, and a list of forces that are operative.  Each force
--   has a list of edges to which it applies, and is represented by a
--   function giving the force between any two points.
data Ensemble v = Ensemble { _forces    :: [([Edge], Point v -> Point v -> v)]
                           , _particles :: M.Map Int (Particle v)
                           }

mkLabels [''Ensemble]

------------------------------------------------------------
--  Simulation internals
------------------------------------------------------------

-- | Simulate one time step for an entire ensemble, with the given
--   damping factor.
ensembleStep :: VectorSpace v => ForceLayoutOpts v -> Ensemble v -> Ensemble v
ensembleStep opts = L.modify particles (M.map (particleStep (damping opts))) . recalcForces

-- | Simulate one time step for a particle (assuming the force acting
--   on it has already been computed), with the given damping factor.
particleStep :: VectorSpace v => Scalar v -> Particle v -> Particle v
particleStep d = stepPos . stepVel
  where stepVel p = L.set vel (d *^ (L.get vel p ^+^ L.get force p)) p
        stepPos p = L.modify pos (.+^ L.get vel p) p

-- | Recalculate all the forces acting in the next time step of an
--   ensemble.
recalcForces :: forall v. AdditiveGroup v => Ensemble v -> Ensemble v
recalcForces = calcForces . zeroForces
  where zeroForces = L.modify particles . M.map $ L.set force zeroV
        calcForces (Ensemble fs ps) 
          = Ensemble fs
            (ala Endo foldMap (concatMap (\(es, f) -> (map (mkForce f) es)) fs) ps)
        mkForce :: (Point v -> Point v -> v) -> Edge -> M.Map Int (Particle v) -> M.Map Int (Particle v)
        mkForce f (i1, i2) m
          = case (M.lookup i1 m, M.lookup i2 m) of
              (Just p1, Just p2) -> 
                ( M.adjust (L.modify force (^+^ f (L.get pos p1) (L.get pos p2))) i1
                . M.adjust (L.modify force (^-^ f (L.get pos p1) (L.get pos p2))) i2)
                m
              _                  -> m

-- | Compute the total kinetic energy of an ensemble.
kineticEnergy :: (InnerSpace v, Num (Scalar v)) => Ensemble v -> Scalar v
kineticEnergy = ala Sum F.foldMap . fmap (magnitudeSq . L.get vel) . L.get particles

------------------------------------------------------------
--  Simulation
------------------------------------------------------------

-- | Options for customizing simulation.
data ForceLayoutOpts v = 
  FLOpts 
  { damping     :: Scalar v           -- ^ Damping factor to be
                                      --   applied at each step.
                                      --   Should be between 0 and 1.
  , energyLimit :: Maybe (Scalar v)   -- ^ Kinetic energy below which
                                      --   simulation should stop.
  , stepLimit   :: Maybe Int          -- ^ Maximum number of
                                      --   simulation steps.
  }

-- | Simulate a starting ensemble according to the given options,
--   producing a list of all the intermediate ensembles.  Useful for,
--   /e.g./, making an animation.  Note that if a 'stepLimit' is not
--   specified, the resulting list could be infinite, if the kinetic
--   energy never falls below the specified threshold (or if no energy
--   threshold is given).
simulate :: (InnerSpace v, Ord (Scalar v), Num (Scalar v)) 
         => ForceLayoutOpts v -> Ensemble v -> [Ensemble v]
simulate opts e
  = (e:)
  . takeWhile (maybe (const True) (>) (energyLimit opts) . kineticEnergy)
  . maybe id take (stepLimit opts)
  . drop 1
  . iterate (ensembleStep opts)
  $ e   

-- | Run a simluation from a starting ensemble, yielding either the
--   first ensemble to have kinetic energy below the 'energyLimit' (if
--   given), or the ensemble that results after a number of steps
--   equal to the 'stepLimit' (if given).  Otherwise @forceLayout@
--   will not terminate.
forceLayout :: (InnerSpace v, Ord (Scalar v), Num (Scalar v)) 
            => ForceLayoutOpts v -> Ensemble v -> Ensemble v
forceLayout opts = last . simulate opts

------------------------------------------------------------
--  Standard forces
------------------------------------------------------------

-- | @distForce f p1 p2@ computes the force between two points as a
--   multiple of the unit vector in the direction from @p1@ to @p2@,
--   given a function @f@ which computes the force's magnitude as a
--   function of the distance between the points.
distForce :: (InnerSpace v, Floating (Scalar v)) 
          => (Scalar v -> Scalar v) -> Point v -> Point v -> v
distForce f p1 p2 = withLength (f (distance p1 p2)) (p2 .-. p1)
  where withLength s v = s *^ normalized v

-- | @hookeForce k l p1 p2@ computes the force on @p1@, assuming that
--   @p1@ and @p2@ are connected by a spring with equilibrium length @l@
--   and spring constant @k@.
hookeForce :: (InnerSpace v, Floating (Scalar v)) 
           => Scalar v -> Scalar v -> Point v -> Point v -> v
hookeForce k l = distForce (\d -> k * (d - l))

-- | @coulombForce k@ computes the electrostatic repulsive force
--   between two charged particles, with constant of proportionality
--   @k@.
coulombForce :: (InnerSpace v, Floating (Scalar v)) 
             => Scalar v -> Point v -> Point v -> v
coulombForce k = distForce (\d -> -k * 1/(d*d))