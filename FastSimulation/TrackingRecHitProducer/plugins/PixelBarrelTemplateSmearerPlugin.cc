
// SiPixel Gaussian Smearing
#include "FastSimulation/TrackingRecHitProducer/interface/PixelTemplateSmearerBase.h"
#include "FastSimulation/TrackingRecHitProducer/interface/TrackingRecHitAlgorithmFactory.h"
#include "FastSimulation/TrackingRecHitProducer/interface/TrackingRecHitProduct.h"
#include "FastSimulation/TrackingRecHitProducer/interface/PixelResolutionHistograms.h"

// Geometry
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

// Famos
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"
#include "FastSimulation/Utilities/interface/SimpleHistogramGenerator.h"

// STL

// ROOT
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>


using namespace std;

class PixelBarrelTemplateSmearerPlugin:
    public PixelTemplateSmearerBase
{
    public:
        explicit PixelBarrelTemplateSmearerPlugin(
            const std::string& name,
            const edm::ParameterSet& config,
            edm::ConsumesCollector& consumesCollector
        );
        virtual ~PixelBarrelTemplateSmearerPlugin();

    private:
        void initialize();
};


PixelBarrelTemplateSmearerPlugin::PixelBarrelTemplateSmearerPlugin(
    const std::string& name,
    const edm::ParameterSet& config,
    edm::ConsumesCollector& consumesCollector
):
    PixelTemplateSmearerBase(name,config,consumesCollector)
{
    isForward = false;
    
    //--- Resolution file names.
    theBigPixelResolutionFileName = config.getParameter<string>( "BigPixelBarrelResolutionFile" );
    theEdgePixelResolutionFileName =  config.getParameter<string>( "EdgePixelBarrelResolutionFile" );
    theRegularPixelResolutionFileName = config.getParameter<string>( "RegularPixelBarrelResolutionFile" );

    //--- Merging info.
    theMergingProbabilityFileName = config.getParameter<string>( "MergingProbabilityBarrelFile" );
    theMergingProbabilityFile = new TFile( edm::FileInPath( theMergingProbabilityFileName ).fullPath().c_str()  ,"READ");
    theMergedPixelResolutionXFileName = config.getParameter<string>( "MergedPixelBarrelResolutionXFile" );
    theMergedPixelResolutionXFile = new TFile( edm::FileInPath( theMergedPixelResolutionXFileName ).fullPath().c_str()  ,"READ");
    theMergedPixelResolutionYFileName = config.getParameter<string>( "MergedPixelBarrelResolutionYFile" );
    theMergedPixelResolutionYFile = new TFile( edm::FileInPath( theMergedPixelResolutionYFileName ).fullPath().c_str()  ,"READ");

    initialize();
    
    if (!SiPixelTemplate::pushfile(templateId, thePixelTemp_))
    {
        throw cms::Exception("PixelBarrelTemplateSmearerPlugin:")<<"SiPixel Barrel Template Not Loaded Correctly!"<<endl;
    }

}


PixelBarrelTemplateSmearerPlugin::~PixelBarrelTemplateSmearerPlugin()
{
}


void PixelBarrelTemplateSmearerPlugin::initialize()
{
   rescotAlpha_binMin = -0.2;
   rescotAlpha_binWidth = 0.08 ;
   rescotAlpha_binN = 5;
   rescotBeta_binMin = -5.5;
   rescotBeta_binWidth = 1.0;
   rescotBeta_binN = 11;
   resqbin_binMin = 0;
   resqbin_binWidth = 1;
   resqbin_binN = 4;

  theRegularPixelResolutions =
    new PixelResolutionHistograms( theRegularPixelResolutionFileName.c_str(), "", 1, 0,
                                      1.0,     //  const double cotbetaBinWidth = 0.15;
                                     -5.5,     //  const double cotbetaLowEdge  = 0.0 ;
                                      11,      //  const int    cotbetaBins     = 4;
                                      0.08,    //  const double cotalphaBinWidth        = 0.02 ;
                                      -0.2,    //  const double cotalphaLowEdge = -0.04 ;
                                      5,       //  const int    cotalphaBins    = 4;
                                      1,       //  const int    qbinWidth       = 1;
                                      4);      //  const int    qbins           = 4;

  theBigPixelResolutions =
    new PixelResolutionHistograms( theBigPixelResolutionFileName.c_str(), "", 1, 2,
                                      1.0,     //  const double cotbetaBinWidth = 0.15;
                                     -5.5,     //  const double cotbetaLowEdge  = 0.0 ;
                                      11,      //  const int    cotbetaBins     = 4;
                                      0.08,    //  const double cotalphaBinWidth        = 0.02 ;
                                      -0.2,    //  const double cotalphaLowEdge = -0.04 ;
                                      5,       //  const int    cotalphaBins    = 4;
                                      1,       //  const int    qbinWidth       = 1;
                                      4);      //  const int    qbins           = 4;

  theEdgePixelResolutions =
    new PixelResolutionHistograms( theEdgePixelResolutionFileName.c_str(), "", 1, 1,
                                      1.0,     //  const double cotbetaBinWidth = 0.15;
                                     -5.5,     //  const double cotbetaLowEdge  = 0.0 ;
                                      11,      //  const int    cotbetaBins     = 4;
                                      0.08,    //  const double cotalphaBinWidth        = 0.02 ;
                                      -0.2,    //  const double cotalphaLowEdge = -0.04 ;
                                      5,       //  const int    cotalphaBins    = 4;
                                      1,       //  const int    qbinWidth       = 1;
                                      4);      //  const int    qbins           = 4;
  
}


DEFINE_EDM_PLUGIN(
    TrackingRecHitAlgorithmFactory,
    PixelBarrelTemplateSmearerPlugin,
    "PixelBarrelTemplateSmearerPlugin"
);
