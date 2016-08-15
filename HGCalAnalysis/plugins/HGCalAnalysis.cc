// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"



#include "TTree.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalMultiCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "RecoNtuples/HGCalAnalysis/interface/AEvent.h"
#include "RecoNtuples/HGCalAnalysis/interface/AObData.h"

#include <string>
#include <map>

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalAnalysis(const edm::ParameterSet&);
  ~HGCalAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsHE;
  edm::EDGetTokenT<reco::CaloClusterCollection> _clusters;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;

  TTree                     *tree;
  AEvent                    *event;
  ARecHitCollection         *arhc;
  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  hgcal::RecHitTools         recHitTools;
};



HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  pre(iConfig.getParameter<double>("depthClusteringCone"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");


  if(detector=="both"){
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    _recHitsHE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 1;
  }else if(detector=="EE"){
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 2;
  }else{
    _recHitsHE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _clusters = consumes<reco::CaloClusterCollection>(edm::InputTag("imagingClusterHGCal"));
    algo = 3;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));

  edm::Service<TFileService> fs;
  fs->make<TH1F>("totale", "totale", 100, 0, 5.);
  tree = new TTree("hgc","Analysis");
  arhc = new ARecHitCollection();
  event = new AEvent();
  tree->Branch("event","AEvent",&event,16000,99);
  tree->Branch("rechits","ARecHitCollection",&arhc,16000,0);
}
HGCalAnalysis::~HGCalAnalysis()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  arhc->clear();

  edm::ESHandle<HGCalGeometry> geoHandleEE;
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandleEE);
  const HGCalGeometry& hgcGeoEE = *geoHandleEE;

  edm::ESHandle<HGCalGeometry> geoHandleHE;
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",geoHandleHE);
  const HGCalGeometry& hgcGeoHE = *geoHandleHE;

  recHitTools.getEventSetup(iSetup);

  int npart = 0;
  int nhit  = 0;

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleHE;
  Handle<reco::CaloClusterCollection> clusterHandle;

  std::vector<HGCalMultiCluster> multiClusters;
  iEvent.getByToken(_clusters,clusterHandle);
  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  float vx = 0.;
  float vy = 0.;
  float vz = 0.;

  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }


  //make a map detid-rechit
  std::map<HGCalDetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsHE,recHitHandleHE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      const HGCRecHitCollection& rechitsHE = *recHitHandleHE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[HGCalDetId(rechitsEE[i].detid())] = &rechitsEE[i];
      }
      for(unsigned int i = 0; i < rechitsHE.size(); i++){
	hitmap[HGCalDetId(rechitsHE[i].detid())] = &rechitsHE[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[HGCalDetId(rechitsEE[i].detid())] = &rechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsHE,recHitHandleHE);
      const HGCRecHitCollection& rechitsHE = *recHitHandleHE;
      for(unsigned int i = 0; i < rechitsHE.size(); i++){
	hitmap[HGCalDetId(rechitsHE[i].detid())] = &rechitsHE[i];
      }
      break;
    }
  default:
    break;
  }




  const reco::CaloClusterCollection &clusters = *clusterHandle;
  multiClusters = pre.makePreClusters(clusters);
  unsigned int cluster_index = 0;
  for(unsigned int i = 0; i < multiClusters.size(); i++){
//      int cl2dSeed = 0;
    for(HGCalMultiCluster::component_iterator it = multiClusters[i].begin();
	it!=multiClusters[i].end(); it++){
//      if(it->energy() > (it+cl2dSeed)->energy()) cl2dSeed = it - multiClusters[i].begin();

      const std::vector< std::pair<DetId, float> > &hf = it->hitsAndFractions();
      int ncoreHit = 0;
      int layer = 0;
//      int rhSeed = 0;

      for(unsigned int j = 0; j < hf.size(); j++){
	//here we loop over detid/fraction pairs
	ncoreHit += int(hf[j].second);
	int flags = 0x0;
	if (hf[j].second>0. && hf[j].second<1.)
	  flags = 0x1;
	else if(hf[j].second<0.)
	  flags = 0x3;
	else if(hf[j].second==0.)
	  flags = 0x2;
	const HGCRecHit *hit = hitmap[hf[j].first];
	layer = HGCalDetId(hf[j].first).layer();

    const GlobalPoint position = recHitTools.getPosition(hf[j].first);
    const unsigned int wafer = recHitTools.getWafer(hf[j].first);
    const unsigned int cell  = recHitTools.getCell(hf[j].first);
    const double cellThickness = recHitTools.getSiThickness(hf[j].first);
    const bool isHalfCell = recHitTools.isHalfCell(hf[j].first);
    const double eta = recHitTools.getEta(position, vz);
    const double phi = recHitTools.getPhi(position);
    const double pt = recHitTools.getPt(position, hit->energy(), vz);

//	if(hit->energy() > hitmap[hf[rhSeed].first]->energy()) rhSeed = j;

	if(layer < 29){
	  nhit++;
	}
    arhc->push_back(ARecHit(layer, wafer, cell,
                position.x(), position.y(), position.z(),
                eta,phi,pt,
                hit->energy(), hit->time(), cellThickness,
                isHalfCell, flags, cluster_index));

      }
      cluster_index++;
    }
  }
  event->set(iEvent.run(),iEvent.id().event(),nhit,
	     vx,vy,vz);
  tree->Fill();
}

void
HGCalAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalAnalysis::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
