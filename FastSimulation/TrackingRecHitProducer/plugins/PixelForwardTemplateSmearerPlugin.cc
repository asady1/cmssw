/** PixelForwardTemplateSmearerPlugin.cc
 * ---------------------------------------------------------------------
 * Description:  see PixelForwardTemplateSmearerPlugin.h
 * Authors:  R. Ranieri (CERN), M. Galanti
 * History: Oct 11, 2006 -  initial version
 * 
 * New Pixel Resolution Parameterization 
 * Introduce SiPixelTemplate Object to Assign Pixel Errors
 * by G. Hu
 * ---------------------------------------------------------------------
 */

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
//#include <TAxis.h>

//#define FAMOS_DEBUG

const double microntocm = 0.0001;
using namespace std;


class PixelForwardTemplateSmearerPlugin:
    public PixelTemplateSmearerBase
{
    public:
          explicit PixelForwardTemplateSmearerPlugin(
            const std::string& name,
            const edm::ParameterSet& config,
            edm::ConsumesCollector& consumesCollector
          );
          virtual ~PixelForwardTemplateSmearerPlugin();

    private:
        void initialize();
};


PixelForwardTemplateSmearerPlugin::PixelForwardTemplateSmearerPlugin( 
    const std::string& name,
    const edm::ParameterSet& config,
    edm::ConsumesCollector& consumesCollector 
):
    PixelTemplateSmearerBase(name,config,consumesCollector)
{

    isForward = true;

    //--- Resolution file names.
    theEdgePixelResolutionFileName = config.getParameter<string>( "EdgePixelForwardResolutionFile" );
    theBigPixelResolutionFileName = config.getParameter<string>( "BigPixelForwardResolutionFile" );
    theRegularPixelResolutionFileName = config.getParameter<string>( "RegularPixelForwardResolutionFile" );

    //--- Hit merging info.
    theMergingProbabilityFileName = config.getParameter<string>( "MergingProbabilityForwardFile" );
    theMergingProbabilityFile =new TFile( edm::FileInPath( theMergingProbabilityFileName ).fullPath().c_str()  ,"READ");
    theMergedPixelResolutionXFileName = config.getParameter<string>( "MergedPixelForwardResolutionXFile" );
    theMergedPixelResolutionXFile = new TFile( edm::FileInPath( theMergedPixelResolutionXFileName ).fullPath().c_str()  ,"READ");
    theMergedPixelResolutionYFileName = config.getParameter<string>( "MergedPixelForwardResolutionYFile" );
    theMergedPixelResolutionYFile = new TFile( edm::FileInPath( theMergedPixelResolutionYFileName ).fullPath().c_str()  ,"READ");

    initialize();
    
    
    
    if( ! SiPixelTemplate::pushfile(templateId, thePixelTemp_) )
    {
        throw cms::Exception("PixelForwardTemplateSmearerPlugin:") <<"SiPixel Forward Template Not Loaded Correctly!";
    }
}


PixelForwardTemplateSmearerPlugin::~PixelForwardTemplateSmearerPlugin()
{
}


void PixelForwardTemplateSmearerPlugin::initialize()
{
    rescotAlpha_binMin = 0.1;
    rescotAlpha_binWidth = 0.1;
    rescotAlpha_binN = 4;
    rescotBeta_binMin = 0.;
    rescotBeta_binWidth = 0.15;
    rescotBeta_binN = 4;
    resqbin_binMin = 0;
    resqbin_binWidth = 1;
    resqbin_binN = 4;

  theRegularPixelResolutions =
    new PixelResolutionHistograms( theRegularPixelResolutionFileName.c_str(), "", 0, 0,
                                      0.15,    //  const double cotbetaBinWidth = 0.15;
                                      0.0,     //  const double cotbetaLowEdge  = 0.0 ;
                                      4,       //  const int    cotbetaBins     = 4;
                                      0.1,     //  const double cotalphaBinWidth        = 0.02 ;
                                      0.1,     //  const double cotalphaLowEdge = -0.04 ;
                                      4,       //  const int    cotalphaBins    = 4;
                                      1,     //  const int    qbinWidth       = 1;
                                      4);      //  const int    qbins           = 4;

  theBigPixelResolutions =
    new PixelResolutionHistograms( theBigPixelResolutionFileName.c_str(), "", 0, 2,
                                      0.15,    //  const double cotbetaBinWidth = 0.15;
                                      0.0,     //  const double cotbetaLowEdge  = 0.0 ;
                                      4,       //  const int    cotbetaBins     = 4;
                                      0.1,     //  const double cotalphaBinWidth        = 0.02 ;
                                      0.1,     //  const double cotalphaLowEdge = -0.04 ;
                                      4,       //  const int    cotalphaBins    = 4;
                                      0,     //  const int    qbinWidth       = 1;
                                      0);      //  const int    qbins           = 4;

  theEdgePixelResolutions =
    new PixelResolutionHistograms( theEdgePixelResolutionFileName.c_str(), "", 0, 1,
                                      0.15,    //  const double cotbetaBinWidth = 0.15;
                                      0.0,     //  const double cotbetaLowEdge  = 0.0 ;
                                      4,       //  const int    cotbetaBins     = 4;
                                      0.1,     //  const double cotalphaBinWidth        = 0.02 ;
                                      0.1,     //  const double cotalphaLowEdge = -0.04 ;
                                      4,       //  const int    cotalphaBins    = 4;
                                      1,     //  const int    qbinWidth       = 1;
                                      4);      //  const int    qbins           = 4;

}

DEFINE_EDM_PLUGIN(
    TrackingRecHitAlgorithmFactory,
    PixelForwardTemplateSmearerPlugin,
    "PixelForwardTemplateSmearerPlugin"
);
