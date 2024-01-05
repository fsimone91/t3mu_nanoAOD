from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.met_cff import *

## from PhysicsTools.Tau3muNANO.trgbits_cff import * # modified
## 
## ## for gen and trigger muon
## from PhysicsTools.Tau3muNANO.genparticlesT3m_cff import * # define new
## from PhysicsTools.Tau3muNANO.particlelevelT3m_cff import * # define new
## from PhysicsTools.Tau3muNANO.triggerObjectsTau3Mu_cff import * #define new
## from PhysicsTools.Tau3muNANO.muonsTau3mu_cff import * # define new 

from PhysicsTools.HFt3muNANO.trimu_cff import * #define new
from PhysicsTools.HFt3muNANO.dimu_cff import * #define new
from PhysicsTools.HFt3muNANO.muon_cff import * #define new
from PhysicsTools.HFt3muNANO.gen_cff import * #define new
from PhysicsTools.HFt3muNANO.particlelevelT3m_cff import * # define new
from PhysicsTools.HFt3muNANO.l1TrigObj_cff import * # define new

vertexTable.svSrc = cms.InputTag("slimmedSecondaryVertices")

nanoSequence = cms.Sequence(nanoMetadata + 
                            cms.Sequence(vertexTask) + cms.Sequence(metTablesTask) +           
                            cms.Sequence(globalTablesTask) + cms.Sequence(vertexTablesTask)
                           )
nanoSequenceMC = cms.Sequence(particleLevelT3mSequence + genParticleT3MSequence + cms.Sequence(metMCTask)
                             + cms.Sequence(globalTablesMCTask) + cms.Sequence(genWeightsTableTask) + genParticleT3MTables + lheInfoTable
                           )

def nanoAOD_customizeTables(process):
  #  process.nanoSequence = cms.Sequence( process.nanoSequence + muonT3mSequence + muonT3mTables)
    process.nanoSequence = cms.Sequence( process.nanoSequence + cms.Sequence(l1TablesTask) )
    return process

def nanoAOD_customizeHFt3mu(process):
    process.nanoHFt3muSequence = cms.Sequence( Tau3MuSequence + Tau3MuTableSequence + DimuonSequence + DimuonTableSequence + muonSequence + muonTableSequence )
    return process

def nanoAOD_customizeMC(process):
  for name, path in process.paths.iteritems():
    path.insert(0, nanoSequenceMC)

## from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
## def nanoAOD_customizeMC(process):
##     for name, path in process.paths.iteritems():
##         # replace all the non-match embedded inputs with the matched ones
##         massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
##         #massSearchReplaceAnyInputTag(path, 'tracksX:SelectedTracks', 'tracksXMCMatchEmbedded')
## 
##         # modify the path to include mc-specific info
##         path.insert(0, nanoSequenceMC)
##         path.replace(process.muonT3mSequence, process.muonT3mMC)
##         #path.replace(process.tracksXSequence, process.tracksXMC)
