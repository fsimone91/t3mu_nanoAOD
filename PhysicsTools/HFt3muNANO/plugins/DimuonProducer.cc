#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"

//#include "PhysicsTools/HFt3muNANO/interface/KinematicFitResult.h"
#include "PhysicsTools/HFt3muNANO/interface/Utils.h"

#include <TLorentzVector.h>
#include <TVector.h>
#include <TVector3.h>
#include <TMatrix.h>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>


// 
// DimuonProducer is designed for HF tau->3mu analysis
// 
typedef std::vector<reco::TransientTrack> TransientTrackCollection;

namespace{
  const float muon_mass_     = 0.10565837;
  const float kaon_mass_     = 0.493677;
  const float mass_err_      = 1.6e-5;
  const float pion_mass_     = 0.139570;
  const float jpsi_mass_     = 3.0969;
  const float proton_mass_   = 0.9382;

  struct KalmanFitResult{
    TransientVertex refitVertex;
    std::vector<reco::TransientTrack> refitTtracks;
    
    bool valid() const {
      return refitVertex.isValid() and refitVertex.hasRefittedTracks() and refitTtracks.size()>1;
    }
  
    float refit_mass() const
    {
      if ( not valid() ) return -1.0;
      reco::Track track1 = refitTtracks.at(0).track();
      reco::Track track2 = refitTtracks.at(1).track();

      TLorentzVector d1, d2, mother;
      d1.SetPxPyPzE(track1.px(), track1.py(), track1.pz(), sqrt(pow(track1.p(), 2.0) + pow(muon_mass_, 2.0)));
      d2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), sqrt(pow(track2.p(), 2.0) + pow(muon_mass_, 2.0)));
      mother = d1 + d2;
      return mother.M();
    }

    float refit_massErr() const
    {
      if ( not valid() ) return -1.0;
      //WIP
      return -1.0;
    }
    
    float chi2() const
    {
      if ( not valid() ) return -1.0;
      return refitVertex.totalChiSquared();
    }
    
    float ndof() const
    {
      return refitVertex.degreesOfFreedom();
    }
    
    float vtxProb() const
    {
      if ( not valid() ) return -1.0;
      return TMath::Prob((double)refitVertex.totalChiSquared(), int(rint(refitVertex.degreesOfFreedom())));
    }
    
  };
}

using namespace std;

///////////////////////////////////////////////////////////////////////////
/////                             P L U G I N
/////////////////////////////////////////////////////////////////////////////

class DimuonProducer : public edm::stream::EDProducer<> {
  
  public:
  
    typedef std::vector<pat::Muon> MuonCollection;
    //typedef std::vector<reco::TransientTrack> TransientTrackCollection;
    typedef std::vector<pat::PackedCandidate> PackedCandidatesCollection;

    // constructor
    explicit DimuonProducer(const edm::ParameterSet &cfg);
    ~DimuonProducer() override {};

  
  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);
    const StringCutObjectSelector<pat::Muon> muon_selection_;    

    const edm::EDGetTokenT<MuonCollection> src_;
    //const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
    const edm::EDGetTokenT<PackedCandidatesCollection> pkdCand_src_;

    // Transient Track Builder
    utils::ESTokenHandle<TransientTrackBuilder, TransientTrackRecord> m_transientTrackBuilder;

    const bool debug = false;

    pat::CompositeCandidate 
    getDimuon(const edm::Event& iEvent,
	       const edm::Ptr<pat::Muon> mu1,
               const edm::Ptr<pat::Muon> mu2) const;

    KalmanFitResult
    refitter(const edm::Event& iEvent,
             const reco::TransientTrack& cand1,
             const reco::TransientTrack& cand2);

  //  KinematicFitResult 
  //  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
  //                            std::vector<float> masses) const;

    void
    addFitInfo(pat::CompositeCandidate& cand,
               const KalmanFitResult& fit,
               std::string name );
};

DimuonProducer::DimuonProducer(const edm::ParameterSet& iConfig):
   muon_selection_{iConfig.getParameter<std::string>("muonSelection")},
   src_{consumes<MuonCollection>( iConfig.getParameter<edm::InputTag>("src") )},
   //ttracks_src_{consumes<TransientTrackCollection>( iConfig.getParameter<edm::InputTag>("transientTracksSrc") )},
   pkdCand_src_{consumes<PackedCandidatesCollection>( iConfig.getParameter<edm::InputTag>("packedCandidatesSrc") )},
   m_transientTrackBuilder{consumesCollector(), "TransientTrackBuilder"}
   {
     produces<pat::CompositeCandidateCollection>("SelectedDimuons");
   }


void DimuonProducer::produce(edm::Event &evt, const edm::EventSetup& iSetup) {

  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  //edm::Handle<TransientTrackCollection> ttracks;
  //evt.getByToken(ttracks_src_, ttracks);

  edm::Handle<PackedCandidatesCollection> PFcands;
  evt.getByToken(pkdCand_src_, PFcands);

  m_transientTrackBuilder.getFromES(iSetup);

  // output collection
  auto dimuons  = std::make_unique<pat::CompositeCandidateCollection>();

  // build 2mu candidates
  auto nPFcands = PFcands->size();
  auto nMuons = muons->size();

  cout << "n muons in event" << nMuons << endl;
  // loop on muons to build candidates
  // muon_selection_ is a configurable pre-selection
  // like ('pt > 0.5 &&  abs(eta)<2.4 && (innerTrack().isNonnull) && (charge!=0) && (innerTrack().hitPattern().numberOfValidPixelHits()>0)'), 
  if ( nMuons > 1 ){
    for ( unsigned int i = 0; i < nMuons; ++i ) {
      const edm::Ptr<pat::Muon> l1_ptr(muons, i);
      if(!muon_selection_(*l1_ptr)) continue;
      for ( unsigned int j = i+1; j < nMuons; ++j ) {
        const edm::Ptr<pat::Muon> l2_ptr(muons, j);
        if(!muon_selection_(*l2_ptr)) continue;
        //cut on total charge, looking for OS pairs
        if( std::abs(l1_ptr->charge() + l2_ptr->charge() )!= 0 ) continue; 
        pat::CompositeCandidate dimuCand;
        dimuCand = getDimuon(evt, l1_ptr, l2_ptr); 
        dimuCand.addUserInt("mu1_idx", i );
        dimuCand.addUserInt("mu2_idx", j );
        cout << "dimuon in event" << endl;

        // refit SV
        KalmanFitResult sv;
        reco::TransientTrack ttrk1=m_transientTrackBuilder->build( *l1_ptr->bestTrack() );
        reco::TransientTrack ttrk2=m_transientTrackBuilder->build( *l2_ptr->bestTrack() );
        sv = refitter(evt, ttrk1, ttrk2);
        if(!sv.valid()) continue;

        addFitInfo(dimuCand, sv, "kalfit");
        // push in the event
        dimuons->push_back(dimuCand);

      }
    }
  }

  evt.put(std::move(dimuons),  "SelectedDimuons");
}

pat::CompositeCandidate
DimuonProducer::getDimuon(const edm::Event& iEvent,
	                    const edm::Ptr<pat::Muon> mu1,
                            const edm::Ptr<pat::Muon> mu2)
const {
  pat::CompositeCandidate tauCand;

  tauCand.setP4(mu1->p4() + mu2->p4());
  tauCand.setCharge(mu1->charge() + mu2->charge());
  tauCand.addUserInt("charge", tauCand.charge());

  // Muon charge
  tauCand.addUserInt("mu1_charge", mu1->charge() );
  tauCand.addUserInt("mu2_charge", mu2->charge() );

  // Muon kinematics
  tauCand.addUserFloat( "mu1_pt",  mu1->pt() );
  tauCand.addUserFloat( "mu1_eta", mu1->eta() );
  tauCand.addUserFloat( "mu1_phi", mu1->phi() );
  tauCand.addUserFloat( "mu2_pt",  mu2->pt() );
  tauCand.addUserFloat( "mu2_eta", mu2->eta() );
  tauCand.addUserFloat( "mu2_phi", mu2->phi() );

  // Invariant mass before refit
  TLorentzVector m1, m2, tau;
  m1.SetPxPyPzE(mu1->px(), mu1->py(), mu1->pz(), mu1->energy() );
  m2.SetPxPyPzE(mu2->px(), mu2->py(), mu2->pz(), mu2->energy() );
  tau = m1 + m2 ;
  tauCand.addUserFloat( "mass", tau.M() );

  // Use UserCands as they should not use memory but keep the Ptr itself
  tauCand.addUserCand("mu1", mu1 );
  tauCand.addUserCand("mu2", mu2 );    

  return tauCand;
}

KalmanFitResult
DimuonProducer::refitter(const edm::Event& iEvent,
	                  const reco::TransientTrack& trk1,
                          const reco::TransientTrack& trk2)
{
  std::vector<float> masses;
  std::vector<reco::TransientTrack> ttrks_prefit;
  ttrks_prefit.push_back(trk1);
  ttrks_prefit.push_back(trk2);

  //refit 2mu vertex
  KalmanFitResult result;
  KalmanVertexFitter SVfitter (true);
  TransientVertex SV = SVfitter.vertex(ttrks_prefit);
  vector <reco::TransientTrack> ttrks_postfit = SV.refittedTracks(); 
  result.refitVertex = SV;
  result.refitTtracks = ttrks_postfit;
  return result;
}


void
DimuonProducer::addFitInfo( pat::CompositeCandidate& cand, const KalmanFitResult& fit, std::string name )
{
  cand.addUserInt(   name+"_valid",       fit.valid() );
  cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
  cand.addUserFloat( name+"_vtx_chi2dof", fit.chi2()>0?fit.chi2()/fit.ndof():-1);
  cand.addUserFloat( name+"_mass",        fit.refit_mass() );
  cand.addUserFloat( name+"_massErr",     fit.refit_massErr() );

  reco::Track trk1 = fit.refitTtracks.at(0).track();
  reco::Track trk2 = fit.refitTtracks.at(1).track();

  cand.addUserFloat( name+"_mu1_pt",      trk1.pt() );
  cand.addUserFloat( name+"_mu1_eta",     trk1.eta() );
  cand.addUserFloat( name+"_mu1_phi",     trk1.phi() );
  cand.addUserFloat( name+"_mu2_pt",      trk2.pt() );
  cand.addUserFloat( name+"_mu2_eta",     trk2.eta() );
  cand.addUserFloat( name+"_mu2_phi",     trk2.phi() );
  // cand.addUserFloat( name+"_lxy",         fit.lxy );
  // cand.addUserFloat( name+"_sigLxy",      fit.sigLxy );
  // cand.addUserFloat( name+"_cosAlphaXY",  fit.cosAlpha );
  // cand.addUserFloat( name+"_sipBS",       fit.ipSigBS );
  // cand.addUserFloat( name+"_sipPV",       fit.ipSigPV );
  // cand.addUserFloat( name+"_pt",          fit.p3().perp() );
  // cand.addUserFloat( name+"_eta",         fit.p3().eta() );
  // cand.addUserFloat( name+"_phi",         fit.p3().phi() );
}    

DEFINE_FWK_MODULE(DimuonProducer);
