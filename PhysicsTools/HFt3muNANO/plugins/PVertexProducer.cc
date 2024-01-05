
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

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

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

#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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
// PVertexProducer is designed for HF tau->3mu analysis
// 
typedef std::vector<reco::TransientTrack> TransientTrackCollection;

using namespace std;

///////////////////////////////////////////////////////////////////////////
/////                             P L U G I N
/////////////////////////////////////////////////////////////////////////////

class PVertexProducer : public edm::stream::EDProducer<> {
  
  public:
  
    typedef std::vector<pat::Muon> MuonCollection;
    //typedef std::vector<reco::TransientTrack> TransientTrackCollection;
    typedef std::vector<pat::PackedCandidate> PackedCandidatesCollection;

    // constructor
    explicit PVertexProducer(const edm::ParameterSet &cfg);
    ~PVertexProducer() override {};

  
  private:

    virtual void produce(edm::Event&, const edm::EventSetup&);

    const edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > src_;
    const edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
    //const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
    const edm::EDGetTokenT<PackedCandidatesCollection> pkdCand_src_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpot_src_;

    // Transient Track Builder
    utils::ESTokenHandle<TransientTrackBuilder, TransientTrackRecord> m_transientTrackBuilder;

    const bool debug = false;


    bool tracksMatchByDeltaR(const reco::Track* trk1, const reco::Track* trk2);
    bool tracksMatchByDeltaR2(const reco::TransientTrack trk1, const reco::Track* trk2);
    void removeTracks3(vector<reco::TransientTrack> &pvTracks, const std::vector<reco::Track*> svTracks);

    void
    addVertexInfo(TransientVertex vtx,
                  pat::CompositeCandidate& cand,
                  std::string name );

    void
    addDisplacementInfo(reco::Vertex v1, 
                        VertexState v2,
                        pat::CompositeCandidate& cand,
                        std::string name );
};

PVertexProducer::PVertexProducer(const edm::ParameterSet& iConfig):
   src_{consumes<edm::View<pat::CompositeCandidate> >( iConfig.getParameter<edm::InputTag>("src") )}, //takes SelectedTriplets as input
   vertexToken_{consumes<edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) },
   pkdCand_src_{consumes<PackedCandidatesCollection>( iConfig.getParameter<edm::InputTag>("packedCandidatesSrc") )},
   beamSpot_src_{consumes<reco::BeamSpot>( iConfig.getParameter<edm::InputTag>("offlineBeamSpot") )},
   m_transientTrackBuilder{consumesCollector(), "TransientTrackBuilder"}
   {
     produces<pat::CompositeCandidateCollection>("RefittedPV");
   }


void PVertexProducer::produce(edm::Event &evt, const edm::EventSetup& iSetup) {

  // input
  edm::Handle<edm::View<pat::CompositeCandidate> > triplets;
  evt.getByToken(src_, triplets);
  
  //edm::Handle<TransientTrackCollection> ttracks;
  //evt.getByToken(ttracks_src_, ttracks);

  edm::Handle< edm::View<reco::Vertex> >vertices;
  evt.getByToken(vertexToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found

  edm::Handle<PackedCandidatesCollection> PFcands;
  evt.getByToken(pkdCand_src_, PFcands);

  edm::Handle<reco::BeamSpot> beamSpots;
  evt.getByToken(beamSpot_src_, beamSpots);

  m_transientTrackBuilder.getFromES(iSetup);

  // output collection
  auto modified_triplets  = std::make_unique<pat::CompositeCandidateCollection>();

  auto nTriplets = triplets->size();
  auto nPV = vertices->size();

  cout << "n triplet candidates in event" << nTriplets << endl;
  cout << "n PV =" << nPV << endl;

  if ( ! beamSpots.isValid() ) return; // skip event if no BS
  const reco::BeamSpot& beamSpot = *beamSpots.product();

  std::vector<uint> vtxKeys;
  vector<pat::PackedCandidate> MyPFcands;
  uint kk = 0;

  // first select only PV associated to some PF candidates
  for (std::vector<pat::PackedCandidate>::const_iterator cand = PFcands->begin(); cand != PFcands->end(), kk!= PFcands->size(); ++cand, ++kk) {
    if (cand->charge()==0 || cand->vertexRef().isNull() ) continue;
    int key = cand->vertexRef().key();
    int quality = cand->pvAssociationQuality();
    if(cand->fromPV(cand->vertexRef().key())<2) continue;
    if(cand->fromPV(cand->vertexRef().key())==2 && quality!=pat::PackedCandidate::UsedInFitLoose) continue;
    
    if ( !(cand->bestTrack()) ) continue;
    //if (quality != pat::PackedCandidate::UsedInFitTight)  continue;
    //if (quality != pat::PackedCandidate::UsedInFitLoose)  continue;
    vtxKeys.push_back(key);
    MyPFcands.push_back(*cand);
  }

  if( !(vtxKeys.size()>0 && vertices->size()>0) ) return;

  // sort and clean list
  sort( vtxKeys.begin(), vtxKeys.end() );
  vtxKeys.erase( unique(vtxKeys.begin(), vtxKeys.end() ), vtxKeys.end() );

  vector<vector<pat::PackedCandidate>> AssoCandToVtx;
  //loop on list of selected primavy vertices to store associated PF candidates
  for(uint i=0; i<vtxKeys.size(); i++){
    vector<pat::PackedCandidate> tmp_cand;

    for (vector<pat::PackedCandidate>::const_iterator myc = MyPFcands.begin(); myc != MyPFcands.end(); ++myc) {
        //cout<<"cand vtx="<<myc->vertexRef().key()<<" cand pt="<<myc->pt()<<endl;
        if(myc->vertexRef().key()==vtxKeys.at(i)) {
            tmp_cand.push_back(*myc);
        }
    }
    AssoCandToVtx.push_back(tmp_cand);
  }//loop vtxKeys
 
  vector<vector<reco::TransientTrack>> transTracksAssoToVtx;
  // extract transient track from associated PF cand
  for(uint i=0; i<AssoCandToVtx.size(); i++){
    std::auto_ptr<reco::TrackCollection> newTrackCollection = std::auto_ptr<reco::TrackCollection>(new reco::TrackCollection);
    std::vector<reco::TransientTrack> transTracks;
  
    for(vector<pat::PackedCandidate>::const_iterator c = AssoCandToVtx.at(i).begin(); c != AssoCandToVtx.at(i).end(); ++c)
        newTrackCollection->push_back(*(c->bestTrack()));

    for (std::vector<reco::Track>::const_iterator iter = newTrackCollection->begin(); iter != newTrackCollection->end(); ++iter){
        reco::TransientTrack tt = m_transientTrackBuilder->build(*iter);
        transTracks.push_back(tt);
    }
    transTracksAssoToVtx.push_back(transTracks);
  }//loop AssoCandToVtx

  // loop on triplet candidates
  for(edm::View<pat::CompositeCandidate>::const_iterator candIt=triplets->begin(); candIt!=triplets->end(); ++candIt){

    const auto * c1 = candIt->daughter(0);
    const pat::Muon *mu1 = dynamic_cast<const pat::Muon *>(c1);
    
    const auto * c2 = candIt->daughter(1);
    const pat::Muon *mu2 = dynamic_cast<const pat::Muon *>(c2);
    
    const auto * c3 = candIt->daughter(2);
    const pat::Muon *mu3 = dynamic_cast<const pat::Muon *>(c3);
    cout<<"PV producer mu1 pt="<<mu1->pt()<<" m2="<<mu2->pt()<<" m3="<<mu3->pt()<<endl;

    const reco::TransientTrack transientTrack1=m_transientTrackBuilder->build( mu1->innerTrack() );
    const reco::TransientTrack transientTrack2=m_transientTrackBuilder->build( mu2->innerTrack() );
    const reco::TransientTrack transientTrack3=m_transientTrackBuilder->build( mu3->innerTrack() );
    reco::Track Track1 =transientTrack1.track();
    reco::Track Track2 =transientTrack2.track();
    reco::Track Track3 =transientTrack3.track();
    vector<reco::Track*> SVTrackRef;

    // SVTrackRef contains transient tracks to be subtracted from PV
    SVTrackRef.push_back(&Track1);
    SVTrackRef.push_back(&Track2);
    SVTrackRef.push_back(&Track3);

    // now select PV based on pointing angle
    double dphi_pv = -1.0;
    uint primaryvertex_index=0;
    uint selVtxId = 0;

    math::Error<3>::type covMatrix;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        covMatrix[i][j] = candIt->userFloat( "kalfit_vtx_err_"+std::to_string(i)+std::to_string(j) );
      }
    }

    reco::Vertex SV = reco::Vertex(candIt->vertex(), covMatrix, candIt->userFloat("kalfit_vtx_chi2"), candIt->userInt("kalfit_vtx_dof"), candIt->numberOfDaughters() );
    TLorentzVector tauVector;
    tauVector.SetPtEtaPhiM(candIt->pt(), candIt->eta(), candIt->phi(), candIt->mass());

    // loop on vertices
    for(uint vtxIndex=0; vtxIndex<vertices->size(); vtxIndex++ ){
      // loop on vtxkeys
      for(uint k=0; k<vtxKeys.size(); k++){
        if(vtxKeys[k]==vtxIndex){
          // build 3D vector distance between tau and PV
          TVector3 Dv3D_reco(SV.x() - (*vertices)[vtxIndex].x(), SV.y() - (*vertices)[vtxIndex].y(), SV.z() - (*vertices)[vtxIndex].z());
          // compute pointing angle
          double Cosdphi_3D = Dv3D_reco.Dot(tauVector.Vect())/(Dv3D_reco.Mag()*tauVector.Vect().Mag());
          if(Cosdphi_3D>dphi_pv){
             dphi_pv = Cosdphi_3D;
             primaryvertex_index=vtxIndex;
             selVtxId=k;
          }
        }
      }
    }// loop on vertices

    cout<<"Max Cosdphi_3D="<<dphi_pv<<" selVtxId="<<selVtxId<<" primaryvertex_index="<<primaryvertex_index<<endl;
    cout<<"Closest PV index before refit: "<<primaryvertex_index<<" x="<<(*vertices)[primaryvertex_index].x()<<" y="<<(*vertices)[primaryvertex_index].y()<<" z="<<(*vertices)[primaryvertex_index].z()<<endl;

    std::vector<reco::TransientTrack> pvTracks_original;
    std::map<const reco::Track*, reco::TransientTrack> pvTrackMap_refit;

    TransientVertex PVertex_prelim;
    bool valid_PV;

    // preliminary fit of PV if it has at least 2 associated tracks
    if(transTracksAssoToVtx.at(selVtxId).size()>1){
      valid_PV = true;
      KalmanVertexFitter PVertex_fitter_prelim(true);
      PVertex_prelim = PVertex_fitter_prelim.vertex(transTracksAssoToVtx.at(selVtxId));
    }
    else valid_PV = false;

    pat::CompositeCandidate candIt_mod = *candIt;
    TransientVertex PVertex;

    if(valid_PV && PVertex_prelim.isValid()){

      // fill pre-refit PV information!
      addVertexInfo(PVertex_prelim, candIt_mod, "PVprelim");
      candIt_mod.addUserInt("PVprelim_ntracks", transTracksAssoToVtx.at(selVtxId).size());

      // remove SV tracks from PV
      removeTracks3(transTracksAssoToVtx.at(selVtxId), SVTrackRef);

      // refit PV
      KalmanVertexFitter PVertex_fitter(true);
      PVertex = PVertex_fitter.vertex(transTracksAssoToVtx.at(selVtxId));

      if(PVertex.isValid()){
        // fill refit PV information!
        addVertexInfo(PVertex, candIt_mod, "PVrefitted");
        candIt_mod.addUserInt("PVrefitted_ntracks", transTracksAssoToVtx.at(selVtxId).size());
      }
    }else{ continue; }

    cout<<" SV coordinates "<<SV.position()<<endl;
    cout<<" RefittedPV done"<<endl;

    // computing displacements and errors
    VertexState PVstate(PVertex.position(), PVertex.positionError());
    VertexState BSstate(beamSpot);

    addDisplacementInfo(SV, PVstate, candIt_mod, "PV");
    addDisplacementInfo(SV, BSstate, candIt_mod, "BS");

    // push in the event
    modified_triplets->push_back(candIt_mod);
  }//  loop on triplet candidates
  evt.put(std::move(modified_triplets),  "RefittedPV");
}

bool PVertexProducer::tracksMatchByDeltaR(const reco::Track* trk1, const reco::Track* trk2)
{
  if ( reco::deltaR(*trk1, *trk2) < 1.e-2 && trk1->charge() == trk2->charge() ) return true;
  else return false;
}

bool PVertexProducer::tracksMatchByDeltaR2(const reco::TransientTrack trk1, const reco::Track* trk2)
{
  if ( reco::deltaR(trk1.track(), *trk2) < 1.e-2 && trk1.track().charge() == trk2->charge() ) return true;
  else return false;
}

void PVertexProducer::removeTracks3(vector<reco::TransientTrack> &pvTracks, const std::vector<reco::Track*> svTracks)
{
  for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin(); svTrack != svTracks.end(); ++svTrack ){
    for(uint f=0;f<pvTracks.size(); f++){
      if ( tracksMatchByDeltaR2(pvTracks.at(f), *svTrack) ) {
        pvTracks.erase(pvTracks.begin()+f);
        break;
      }
    }
  }
}
    


void
PVertexProducer::addVertexInfo(TransientVertex vtx, pat::CompositeCandidate& cand, std::string name )
{
  cand.addUserInt(   name+"_valid",  vtx.isValid() );
  cand.addUserFloat( name+"_chi2",   vtx.totalChiSquared() );
  cand.addUserInt( name+"_dof",      vtx.degreesOfFreedom() );

  GlobalPoint vtxPos(vtx.position());
  cand.addUserFloat( name+"_x",      vtxPos.x() );
  cand.addUserFloat( name+"_y",      vtxPos.y() );
  cand.addUserFloat( name+"_z",      vtxPos.z() );

  float vtxProb = TMath::Prob((double)vtx.totalChiSquared(), int(rint(vtx.degreesOfFreedom())));
  cand.addUserFloat( name+"_prob",   vtxProb );

  // std::vector<double> PV_cov;
  // for (int i = 0; i < 3; i++) {
  //   for (int j = i; j < 3; j++) {
  //     PV_cov.push_back(PVertex.positionError().matrix()[i][j]);
  //   }
  // }
}    

void
PVertexProducer::addDisplacementInfo(reco::Vertex v1, VertexState v2, pat::CompositeCandidate& cand, std::string name){

  VertexDistanceXY vert2D;
  VertexDistance3D vert3D;

  cout<<"distance "<<vert2D.distance(v1, v2).value()<<endl;
  cout<<"error "<<vert2D.distance(v1, v2).error()<<endl;
  cand.addUserFloat( name+"_XY_dist",  vert2D.distance(v1, v2).value() );
  cand.addUserFloat( name+"_XY_err",   vert2D.distance(v1, v2).error() );
  cand.addUserFloat( name+"_XY_sig",   vert2D.distance(v1, v2).significance() );
  cand.addUserFloat( name+"_3D_dist",  vert3D.distance(v1, v2).value() );
  cand.addUserFloat( name+"_3D_err",   vert3D.distance(v1, v2).error() );
  cand.addUserFloat( name+"_3D_sig",   vert3D.distance(v1, v2).significance() );
}

DEFINE_FWK_MODULE(PVertexProducer);
