import FWCore.ParameterSet.Config as cms
from PhysicsTools.HFt3muNANO.common_cff import *

########## inputs preparation ################

Path2022=["HLT_DoubleMu4_3_LowMass", "HLT_DoubleMu3_TkMu_DsTau3Mu_v", "HLT_DoubleMu4_LowMass_Displaced"]
#W channel:
#Path2022=["HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1","HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15","HLT_DoubleMu4_3_LowMass"]

# Tau -> 3mu
diMuonForTau3Mu = cms.EDProducer('DimuonProducer',
    src = cms.InputTag('slimmedMuons'),
    packedCandidatesSrc = cms.InputTag('packedPFCandidates'),
    #transientTracksSrc = cms.InputTag(''),
    # selection definition
    muonSelection = cms.string('((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
)


################################### Tables #####################################

DimuonTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('diMuonForTau3Mu', 'SelectedDimuons'),
    cut = cms.string("mass<10"),
    name = cms.string("dimu"),
    doc = cms.string("Dimuon OS variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        mu1_pt = ufloat("mu1_pt"),
        mu1_eta = ufloat("mu1_eta"),
        mu1_phi = ufloat("mu1_phi"),
        mu2_pt = ufloat("mu2_pt"),
        mu2_eta = ufloat("mu2_eta"),
        mu2_phi = ufloat("mu2_phi"),
        mu1_idx = uint("mu1_idx"),
        mu2_idx = uint("mu2_idx"),
        mu1_charge = uint("mu1_charge"),
        mu2_charge = uint("mu2_charge"),
        mass = ufloat("mass"),
        charge = uint("charge"),
        kalfit_valid      = ubool( "kalfit_valid",       doc = "Kalman SV fit, isValid boolean"),  
        kalfit_vtx_prob   = ufloat("kalfit_vtx_prob",    doc = "Kalman SV fit, vtx probability"), 
        kalfit_vtx_chi2dof= ufloat("kalfit_vtx_chi2dof", doc = "Kalman SV fit, normalised chi2"),
        kalfit_mass       = ufloat("kalfit_mass",        doc = "Kalman SV fit, refitted 2mu mass"), 
        kalfit_massErr    = ufloat("kalfit_massErr",     doc = "Kalman SV fit, refitted 2mu mass error"), 
        kalfit_mu1_pt  = ufloat("kalfit_mu1_pt",   doc = "refitted track kinematics"),
        kalfit_mu1_eta = ufloat("kalfit_mu1_eta",  doc = "refitted track kinematics"),
        kalfit_mu1_phi = ufloat("kalfit_mu1_phi",  doc = "refitted track kinematics"),
        kalfit_mu2_pt  = ufloat("kalfit_mu2_pt",   doc = "refitted track kinematics"),
        kalfit_mu2_eta = ufloat("kalfit_mu2_eta",  doc = "refitted track kinematics"),
        kalfit_mu2_phi = ufloat("kalfit_mu2_phi",  doc = "refitted track kinematics"),
    )
)

########################### Sequencies  ############################

DimuonSequence = cms.Sequence(
    (diMuonForTau3Mu)
)

DimuonTableSequence = cms.Sequence( DimuonTable )
