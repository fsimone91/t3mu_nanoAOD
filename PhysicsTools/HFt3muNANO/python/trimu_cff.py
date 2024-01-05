import FWCore.ParameterSet.Config as cms

from PhysicsTools.HFt3muNANO.common_cff import *

########## inputs preparation ################

Paths=["HLT_DoubleMu4_3_LowMass", "HLT_DoubleMu3_TkMu_DsTau3Mu_v*", "HLT_DoubleMu3_Trk_Tau3mu_v*", "HLT_DoubleMu4_LowMass_Displaced", "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v*"]

#Trigger bit requirement
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
HLTskim = hlt.hltHighLevel.clone()
HLTskim.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
HLTskim.HLTPaths = Paths
HLTskim.andOr = cms.bool( True )
HLTskim.throw = cms.bool( False )

# Tau -> 3mu candidate
looseMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string('pt > 2 &&  abs(eta)<2.4 && (innerTrack().isNonnull) && (charge!=0) && (innerTrack().hitPattern().numberOfValidPixelHits()>0)'), 
    filter = cms.bool(True)                                
)

threeMuonsFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("looseMuons"),
    minNumber = cms.uint32(3),
    #filter = cms.bool(True)
)

threeMuonsCand = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    cut = cms.string('(abs(charge)=1)'),       
    decay = cms.string("looseMuons looseMuons looseMuons")
) 

threeMuonsCandFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("threeMuonsCand"),
    minNumber = cms.uint32(1),
)

muonTripletForTau3Mu = cms.EDProducer('TripletProducer',
    src = cms.InputTag('threeMuonsCand'),
    #packedCandidatesSrc = cms.InputTag('packedPFCandidates'),
    #transientTracksSrc = cms.InputTag(''),
    # selection definition
    muonSelection = cms.string('((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
)

# PV refitter
PVrefitterForTau3Mu = cms.EDProducer('PVertexProducer',
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriplets'),
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    packedCandidatesSrc = cms.InputTag('packedPFCandidates'),
    offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
    #transientTracksSrc = cms.InputTag(''),
)

################################### Tables #####################################

Tau3MuTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('PVrefitterForTau3Mu', 'RefittedPV'),
    cut = cms.string("mass>0 && mass<5"), #can be tightened
    name = cms.string("triplet"),
    doc = cms.string("Triplet variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        mu1_pt  = ufloat("mu1_pt"),
        mu1_eta = ufloat("mu1_eta"),
        mu1_phi = ufloat("mu1_phi"),
        mu2_pt  = ufloat("mu2_pt"),
        mu2_eta = ufloat("mu2_eta"),
        mu2_phi = ufloat("mu2_phi"),
        mu3_pt  = ufloat("mu3_pt"),
        mu3_eta = ufloat("mu3_eta"),
        mu3_phi = ufloat("mu3_phi"),
        #mu1_idx = uint("mu1_idx"),
        #mu2_idx = uint("mu2_idx"),
        #mu3_idx = uint("mu3_idx"),
        mu1_charge = uint("mu1_charge"),
        mu2_charge = uint("mu2_charge"),
        mu3_charge = uint("mu3_charge"),
        mass = ufloat("mass"),
        charge = uint("charge"),
        kalfit_valid      = ubool( "kalfit_valid",       doc = "Kalman SV fit, isValid boolean"),  
        kalfit_vtx_x   = ufloat("kalfit_vtx_x",    doc = "Kalman SV fit, vtx coordinate x"), 
        kalfit_vtx_y   = ufloat("kalfit_vtx_y",    doc = "Kalman SV fit, vtx coordinate y"), 
        kalfit_vtx_z   = ufloat("kalfit_vtx_z",    doc = "Kalman SV fit, vtx coordinate z"), 
        kalfit_vtx_prob   = ufloat("kalfit_vtx_prob",    doc = "Kalman SV fit, vtx probability"), 
        kalfit_vtx_chi2   = ufloat("kalfit_vtx_chi2",    doc = "Kalman SV fit, chi2"),
        kalfit_vtx_dof    = uint("kalfit_vtx_dof",       doc = "Kalman SV fit, degrees of freedom"),
        kalfit_mass       = ufloat("kalfit_mass",        doc = "Kalman SV fit, refitted 3mu mass"), 
        kalfit_massErr    = ufloat("kalfit_massErr",     doc = "Kalman SV fit, refitted 3mu mass error"), 
        kalfit_mu1_pt  = ufloat("kalfit_mu1_pt",   doc = "refitted track kinematics"),
        kalfit_mu1_eta = ufloat("kalfit_mu1_eta",  doc = "refitted track kinematics"),
        kalfit_mu1_phi = ufloat("kalfit_mu1_phi",  doc = "refitted track kinematics"),
        kalfit_mu2_pt  = ufloat("kalfit_mu2_pt",   doc = "refitted track kinematics"),
        kalfit_mu2_eta = ufloat("kalfit_mu2_eta",  doc = "refitted track kinematics"),
        kalfit_mu2_phi = ufloat("kalfit_mu2_phi",  doc = "refitted track kinematics"),
        kalfit_mu3_pt  = ufloat("kalfit_mu3_pt",   doc = "refitted track kinematics"),
        kalfit_mu3_eta = ufloat("kalfit_mu3_eta",  doc = "refitted track kinematics"),
        kalfit_mu3_phi = ufloat("kalfit_mu3_phi",  doc = "refitted track kinematics"),
        PVprelim_valid      = ubool( "PVprelim_valid",       doc = "PV before refit, isValid boolean"),  
        PVprelim_ntracks    = uint( "PVprelim_ntracks",      doc = "PV before refit, number of associated tracks"),  
        PVprelim_x   = ufloat("PVprelim_x",    doc = "PV before refit, vtx coordinate x"), 
        PVprelim_y   = ufloat("PVprelim_y",    doc = "PV before refit, vtx coordinate y"), 
        PVprelim_z   = ufloat("PVprelim_z",    doc = "PV before refit, vtx coordinate z"), 
        PVprelim_prob   = ufloat("PVprelim_prob",    doc = "PV before refit, vtx probability"), 
        PVprelim_chi2   = ufloat("PVprelim_chi2",    doc = "PV before refit, chi2"),
        PVprelim_dof    = uint("PVprelim_dof",       doc = "PV before refit, degrees of freedom"),
        PVrefitted_valid      = ubool( "PVrefitted_valid",       doc = "PV after refit, isValid boolean"),  
        PVrefitted_x   = ufloat("PVrefitted_x",    doc = "PV after refit, vtx coordinate x"), 
        PVrefitted_y   = ufloat("PVrefitted_y",    doc = "PV after refit, vtx coordinate y"), 
        PVrefitted_z   = ufloat("PVrefitted_z",    doc = "PV after refit, vtx coordinate z"), 
        PVrefitted_prob   = ufloat("PVrefitted_prob",    doc = "PV after refit, vtx probability"), 
        PVrefitted_chi2   = ufloat("PVrefitted_chi2",    doc = "PV after refit, chi2"),
        PVrefitted_dof    = uint("PVrefitted_dof",       doc = "PV after refit, degrees of freedom"),
        PVrefitted_ntracks= uint("PVrefitted_ntracks",   doc = "PV after refit, number of associated tracks"),  
    )
)

########################### Sequencies  ############################

Tau3MuSequence = cms.Sequence(
    (HLTskim * looseMuons * threeMuonsFilter * threeMuonsCand * threeMuonsCandFilter * muonTripletForTau3Mu + PVrefitterForTau3Mu)
)

Tau3MuTableSequence = cms.Sequence( Tau3MuTable )
