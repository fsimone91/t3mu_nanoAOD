from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('isPreECALleakage',True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Data taken during ECAL leakage"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

#options.setDefault('maxEvents', -1)
options.setDefault('maxEvents', 10000)
options.setDefault('tag', 'test')
options.parseArguments()

# global tags:
# MC pre ECAL leakage:  124X_mcRun3_2022_realistic_v12
# MC post ECAL leakage: 124X_mcRun3_2022_realistic_postEE_v1
# 2022 ABCDE ReReco : 124X_dataRun3_v15
# 2022 FG Prompt : 124X_dataRun3_PromptAnalysis_v2
if not options.isMC :
    globaltag = '130X_dataRun3_PromptAnalysis_v1'
else :
    globaltag = '130X_mcRun3_2022_realistic_v5' if options.isPreECALleakage else '130X_mcRun3_2022_realistic_postEE_v6'


if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['tau3muNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['xFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles :
    if options.isMC :
        # signal channel
        options.inputFiles = ['/store/mc/Run3Summer22EEMiniAODv3/DstoTau_Tauto3Mu_3MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2550000/0168c233-ed5c-4ed7-91ea-531362f5eef1.root'] if options.isPreECALleakage else \
                            ['/store/mc/Run3Summer22EEMiniAODv3/WtoTauNu_Tauto3Mu_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/975d40c3-629d-41e5-8887-cb34ca21e308.root']
        # control channel
        #options.inputFiles = ['/store/mc/Run3Summer22MiniAODv3/DstoPhiPi_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v2/2810000/0da9edba-f8b9-4e0c-8be1-282cdd2b5685.root'] if options.isPreECALleakage else \
        #                     ['/store/mc/Run3Summer22EEMiniAODv3/DstoPhiPi_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/00589525-be33-4abd-af78-428bb9ace158.root']


    else :
        # data
        options.inputFiles = [
          'root://xrootd-cms.infn.it///store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/355/863/00000/389f9ca1-f590-4691-b7f2-41e0146a8a79.root',
          'root://xrootd-cms.infn.it///store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/355/862/00000/fc972444-ec73-42ea-897c-f2b918fbee7a.root',
          'root://xrootd-cms.infn.it///store/data/Run2022C/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/355/872/00000/2d46ff3d-123e-439e-a2ab-a71627d06e46.root'] if options.isPreECALleakage else \
                             ['root://xrootd-cms.infn.it///store/data/Run2022G/ParkingDoubleMuonLowMass5/MINIAOD/PromptReco-v1/000/362/437/00000/645c4ac6-81f3-4c9d-9e29-21934c2ff83c.root']

annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

#from Configuration.StandardSequences.Eras import eras
process = cms.Process('Tau3muNANO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.HFt3muNANO.nanoTau3Mu_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)
process.options = cms.untracked.PSet(
    #Rethrow 
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(options.wantSummary),
)

# process.nanoMetadata.strings.tag = annotation

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.HFt3muNANO.nanoTau3Mu_cff import *
process = nanoAOD_customizeHFt3mu(process)
process = nanoAOD_customizeTables(process)
#process = PrepMuonCustomNanoAOD(process)

# Path and EndPath definitions
process.nanoAOD_TauTo3mu_step = cms.Path(process.nanoSequence + process.nanoHFt3muSequence)

## MC specific
if options.isMC:
  nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_TauTo3mu_step,
                                process.endjob_step,
                                process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_TauTo3mu_step'
                                   )
)


process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
