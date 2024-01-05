import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import *

# Start with merged particles (pruned + packed)
genForTau3Mu = finalGenParticles.clone(
  src = cms.InputTag("mergedGenParticles"),
  select = cms.vstring(
	"drop *",
        "++keep+ ( abs(pdgId)==13  || abs(pdgId)==15  || abs(pdgId)==11 || abs(pdgId)==211 || abs(pdgId)==321 || abs(pdgId)==12 || abs(pdgId)==14 || abs(pdgId)==16 || abs(pdgId)==431 || abs(pdgId)==511 || abs(pdgId)==521)",  #keep particles, + their daughters, ++ their mothers and grandmothers
   )
)

genParticleForTau3MuTable = genParticleTable.clone(
  src = cms.InputTag("genForTau3Mu"),
  variables = cms.PSet(
      genParticleTable.variables,
      vx = Var("vx()", float, doc="x coordinate of the production vertex position, in cm", precision=10),
      vy = Var("vy()", float, doc="y coordinate of the production vertex position, in cm", precision=10),
      vz = Var("vz()", float, doc="z coordinate of the production vertex position, in cm", precision=10),
  )
)

genParticleT3MSequence = cms.Sequence(genForTau3Mu)
genParticleT3MTables = cms.Sequence(genParticleForTau3MuTable)

