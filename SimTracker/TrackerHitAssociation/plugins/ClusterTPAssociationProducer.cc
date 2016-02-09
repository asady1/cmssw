#include <memory>
#include <vector>
#include <utility>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociationList.h"

class ClusterTPAssociationProducer : public edm::global::EDProducer<>
{
public:
  typedef std::vector<OmniClusterRef> OmniClusterCollection;

  explicit ClusterTPAssociationProducer(const edm::ParameterSet&);
  ~ClusterTPAssociationProducer();

private:
  virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  template <typename T>
  std::vector<std::pair<uint32_t, EncodedEventId> >
  getSimTrackId(const edm::Handle<edm::DetSetVector<T> >& simLinks, const DetId& detId, uint32_t channel) const;

  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink> > sipixelSimLinksToken_;
  edm::EDGetTokenT<edm::DetSetVector<StripDigiSimLink> > sistripSimLinksToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClustersToken_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > stripClustersToken_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
};

ClusterTPAssociationProducer::ClusterTPAssociationProducer(const edm::ParameterSet & cfg)
  : sipixelSimLinksToken_(consumes<edm::DetSetVector<PixelDigiSimLink> >(cfg.getParameter<edm::InputTag>("pixelSimLinkSrc"))),
    sistripSimLinksToken_(consumes<edm::DetSetVector<StripDigiSimLink> >(cfg.getParameter<edm::InputTag>("stripSimLinkSrc"))),
    pixelClustersToken_(consumes<edmNew::DetSetVector<SiPixelCluster> >(cfg.getParameter<edm::InputTag>("pixelClusterSrc"))),
    stripClustersToken_(consumes<edmNew::DetSetVector<SiStripCluster> >(cfg.getParameter<edm::InputTag>("stripClusterSrc"))),
    trackingParticleToken_(consumes<TrackingParticleCollection>(cfg.getParameter<edm::InputTag>("trackingParticleSrc")))
{
  produces<ClusterTPAssociationList>();
}

ClusterTPAssociationProducer::~ClusterTPAssociationProducer() {
}
		
void ClusterTPAssociationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& es) const {
  auto clusterTPList = std::make_unique<ClusterTPAssociationList>();
 
  // Pixel DigiSimLink
  edm::Handle<edm::DetSetVector<PixelDigiSimLink> > sipixelSimLinks;
  //  iEvent.getByLabel(_pixelSimLinkSrc, sipixelSimLinks);
  iEvent.getByToken(sipixelSimLinksToken_,sipixelSimLinks);

  // SiStrip DigiSimLink
  edm::Handle<edm::DetSetVector<StripDigiSimLink> > sistripSimLinks;
  iEvent.getByToken(sistripSimLinksToken_,sistripSimLinks);

  // Pixel Cluster
  edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelClusters;
  bool foundPixelClusters = iEvent.getByToken(pixelClustersToken_,pixelClusters);

  // Strip Cluster
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > stripClusters;
  bool foundStripClusters = iEvent.getByToken(stripClustersToken_,stripClusters);


  // TrackingParticle
  edm::Handle<TrackingParticleCollection>  TPCollectionH;
  iEvent.getByToken(trackingParticleToken_,TPCollectionH);

  // prepare temporary map between SimTrackId and TrackingParticle index
  std::map<std::pair<size_t, EncodedEventId>, TrackingParticleRef> mapping;
  for (TrackingParticleCollection::size_type itp = 0;
                                             itp < TPCollectionH.product()->size(); ++itp) {
    TrackingParticleRef trackingParticle(TPCollectionH, itp);

    // SimTracks inside TrackingParticle
    EncodedEventId eid(trackingParticle->eventId());
    //size_t index = 0;
    for (std::vector<SimTrack>::const_iterator itrk  = trackingParticle->g4Track_begin(); 
                                               itrk != trackingParticle->g4Track_end(); ++itrk) {
      std::pair<uint32_t, EncodedEventId> trkid(itrk->trackId(), eid);
      //std::cout << "creating map for id: " << trkid.first << " with tp: " << trackingParticle.key() << std::endl;
      mapping.insert(std::make_pair(trkid, trackingParticle));
    }
  }

  if ( foundPixelClusters ) {
    // Pixel Clusters 
    for (edmNew::DetSetVector<SiPixelCluster>::const_iterator iter  = pixelClusters->begin(); 
                                                            iter != pixelClusters->end(); ++iter) {
      uint32_t detid = iter->id(); 
      DetId detId(detid);
      edmNew::DetSet<SiPixelCluster> link_pixel = (*iter);
      for (edmNew::DetSet<SiPixelCluster>::const_iterator di  = link_pixel.begin(); 
	   di != link_pixel.end(); ++di) {
	const SiPixelCluster& cluster = (*di);
	edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> c_ref = 
	  edmNew::makeRefTo(pixelClusters, di);
	
	std::set<std::pair<uint32_t, EncodedEventId> > simTkIds; 
	for (int irow = cluster.minPixelRow(); irow <= cluster.maxPixelRow(); ++irow) {
	  for (int icol = cluster.minPixelCol(); icol <= cluster.maxPixelCol(); ++icol) {
	    uint32_t channel = PixelChannelIdentifier::pixelToChannel(irow, icol);
	    std::vector<std::pair<uint32_t, EncodedEventId> > trkid(getSimTrackId<PixelDigiSimLink>(sipixelSimLinks, detId, channel));
	    if (trkid.size()==0) continue; 
	    simTkIds.insert(trkid.begin(),trkid.end());
	  }
	}
	for (std::set<std::pair<uint32_t, EncodedEventId> >::const_iterator iset  = simTkIds.begin(); 
	     iset != simTkIds.end(); iset++) {
	  auto ipos = mapping.find(*iset);
	  if (ipos != mapping.end()) {
	    //std::cout << "cluster in detid: " << detid << " from tp: " << ipos->second.key() << " " << iset->first << std::endl;
	    clusterTPList->push_back(std::make_pair(OmniClusterRef(c_ref), ipos->second));
	  }
	}
      }
    }
  }

  if ( foundStripClusters ) {
    // Strip Clusters
    for (edmNew::DetSetVector<SiStripCluster>::const_iterator iter  = stripClusters->begin(false), eter = stripClusters->end(false); 
	 iter != eter; ++iter) {
      if (!(*iter).isValid()) continue;
      uint32_t detid = iter->id();  
      DetId detId(detid);
      edmNew::DetSet<SiStripCluster> link_strip = (*iter);
      for (edmNew::DetSet<SiStripCluster>::const_iterator di  = link_strip.begin(); 
	   di != link_strip.end(); di++) {
	const SiStripCluster& cluster = (*di);
	edm::Ref<edmNew::DetSetVector<SiStripCluster>, SiStripCluster> c_ref = 
	  edmNew::makeRefTo(stripClusters, di);
	
	std::set<std::pair<uint32_t, EncodedEventId> > simTkIds; 
	int first  = cluster.firstStrip();     
	int last   = first + cluster.amplitudes().size();
	
	for (int istr = first; istr < last; ++istr) {
	  std::vector<std::pair<uint32_t, EncodedEventId> > trkid(getSimTrackId<StripDigiSimLink>(sistripSimLinks, detId, istr));
	  if (trkid.size()==0) continue; 
	  simTkIds.insert(trkid.begin(),trkid.end());
	}
	for (std::set<std::pair<uint32_t, EncodedEventId> >::const_iterator iset  = simTkIds.begin(); 
	     iset != simTkIds.end(); iset++) {
	  auto ipos = mapping.find(*iset);
	  if (ipos != mapping.end()) {
	    //std::cout << "cluster in detid: " << detid << " from tp: " << ipos->second.key() << " " << iset->first << std::endl;
	    clusterTPList->push_back(std::make_pair(OmniClusterRef(c_ref), ipos->second));
	  } 
	}
      } 
    }
  }

  std::sort(clusterTPList->begin(), clusterTPList->end(), clusterTPAssociationListGreater);
  iEvent.put(std::move(clusterTPList));
}

template <typename T>
std::vector<std::pair<uint32_t, EncodedEventId> >
//std::pair<uint32_t, EncodedEventId>
ClusterTPAssociationProducer::getSimTrackId(const edm::Handle<edm::DetSetVector<T> >& simLinks,
                                            const DetId& detId, uint32_t channel) const 
{
  //std::pair<uint32_t, EncodedEventId> simTrkId;
  std::vector<std::pair<uint32_t, EncodedEventId> > simTrkId;
  auto isearch = simLinks->find(detId);
  if (isearch != simLinks->end()) {
    // Loop over DigiSimLink in this det unit
    edm::DetSet<T> link_detset = (*isearch);
    for (typename edm::DetSet<T>::const_iterator it  = link_detset.data.begin(); 
                                                 it != link_detset.data.end(); ++it) {
      if (channel == it->channel()) {
        simTrkId.push_back(std::make_pair(it->SimTrackId(), it->eventId()));
      } 
    }
  }
  return simTrkId;
}
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ClusterTPAssociationProducer);
