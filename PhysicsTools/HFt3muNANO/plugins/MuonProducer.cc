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
// MuonProducer is designed for HF tau->3mu analysis
// 
typedef std::vector<reco::TransientTrack> TransientTrackCollection;

namespace{
  const float muon_mass_     = 0.10565837;
  const float kaon_mass_     = 0.493677;
  const float mass_err_      = 1.6e-5;
  const float pion_mass_     = 0.139570;
  const float jpsi_mass_     = 3.0969;
  const float proton_mass_   = 0.9382;

}

using namespace std;

///////////////////////////////////////////////////////////////////////////
/////                             P L U G I N
/////////////////////////////////////////////////////////////////////////////

class MuonProducer : public edm::stream::EDProducer<> {
  
  public:
  
    typedef std::vector<pat::Muon> MuonCollection;
    //typedef std::vector<reco::TransientTrack> TransientTrackCollection;
    typedef std::vector<pat::PackedCandidate> PackedCandidatesCollection;

    // constructor
    explicit MuonProducer(const edm::ParameterSet &cfg);
    ~MuonProducer() override {};

  
  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);
    const StringCutObjectSelector<pat::Muon> muon_selection_;    

    const edm::EDGetTokenT<MuonCollection> src_;
    //const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;

    const bool debug = false;

};

MuonProducer::MuonProducer(const edm::ParameterSet& iConfig):
   muon_selection_{iConfig.getParameter<std::string>("muonSelection")},
   src_{consumes<MuonCollection>( iConfig.getParameter<edm::InputTag>("src") )}
   {
     produces<std::vector<pat::Muon>>("SelectedMuons");
   }


void MuonProducer::produce(edm::Event &evt, const edm::EventSetup& iSetup) {

  // input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  // output collection
  auto selmuons  = std::make_unique<std::vector<pat::Muon>>();

  auto nMuons = muons->size();

  // loop on muons to build candidates
  // muon_selection_ is a configurable pre-selection
  // like ('pt > 0.5 &&  abs(eta)<2.4 && (innerTrack().isNonnull) && (charge!=0) && (innerTrack().hitPattern().numberOfValidPixelHits()>0)'), 
  if ( nMuons > 2 ){
    for ( unsigned int i = 0; i < nMuons; ++i ) {
      const edm::Ptr<pat::Muon> l1_ptr(muons, i);
      if(!muon_selection_(*l1_ptr)) continue;

      //space here for customisation (i.e. computing MVA ID from Run2)

      // push in the event
      selmuons->push_back(*l1_ptr);

    }
  }

  evt.put(std::move(selmuons),  "SelectedMuons");
}

DEFINE_FWK_MODULE(MuonProducer);
