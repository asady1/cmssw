#ifndef FastSimulation_TrackingRecHitProducer_PixelResolutionHistograms_h
#define FastSimulation_TrackingRecHitProducer_PixelResolutionHistograms_h 1

class TFile;
class TH1F;
class RandomEngineAndDistribution;
class SimpleHistogramGenerator;

#define COTBETA_HIST_MAX  30
#define COTALPHA_HIST_MAX 20
#define QBIN_HIST_MAX      4

class PixelResolutionHistograms {
public:

  //--- Constructor to use when generating resolution histograms.  
  //    We make empty histograms (which we own), but generator pointers 
  //    remain null.
  //
  PixelResolutionHistograms( unsigned int detType,     // where am I?
			     unsigned int ifEdgeOrBig, // is this a regular, edge, or big histogram set
			     double cotbetaBinWidth,   // cot(beta) : bin width,
			     double cotbetaLowEdge,    //           : low endpoint,
			     int	   cotbetaBins,       //           : # of bins
			     double cotalphaBinWidth,  // cot(alpha): bin width, 
			     double cotalphaLowEdge,   //           : low endpoint,
			     int	cotalphaBins,         //           : # of bins
			     int	qbinWidth,            // qbin bin width
			     int	qbins 		  );  // # of qbin bins
  
  //--- Constructor to use when reading the histograms from a file (e.g. when 
  //    inside a running FastSim job).  We get the histograms from a
  //    ROOT file, and we do *not* own them.  But we do own the
  //    generators.
  //
  PixelResolutionHistograms( const char * filename,    // ROOT file to open, with histograms
			     const char * rootdir,     // ROOT dir, "" if none
			     unsigned int detType,     // where am I?
			     unsigned int ifEdgeOrBig, // is this a regular, edge, or big histogram set
			     double cotbetaBinWidth,   // cot(beta) : bin width,
			     double cotbetaLowEdge,    //           : low endpoint,
			     int	   cotbetaBins,       //           : # of bins
			     double cotalphaBinWidth,  // cot(alpha): bin width, 
			     double cotalphaLowEdge,   //           : low endpoint,
			     int	cotalphaBins,         //           : # of bins
			     int	qbinWidth,            // qbin bin width
			     int	qbins 		  );  // # of qbin bins
  

  //--- Destructor (virtual, just in case)
  virtual ~PixelResolutionHistograms();

  //--- Status after construction (esp.loading from file). Non-zero if there
  //    were problems.
  inline int status() { return status_ ; }

  //--- Fill one entry in one resolution histogram.  Use when making histograms.
  int Fill( double dx, double dy,     // the difference wrt true hit 
	    double cotalpha, double cotbeta,  // cotangent of local angles
	    int qbin,                 // Qbin = category for how much charge we have
	    int nxpix, int nypix );   // how long the clusters are in x and y


  //--- Get generators, for resolution in X and Y.  Use in FastSim.
  const SimpleHistogramGenerator * getGeneratorX( double cotalpha,
						  double cotbeta, 
						  int qbin,
						  bool singlex,
						  bool edge);

  const SimpleHistogramGenerator * getGeneratorY( double cotalpha,
						  double cotbeta, 
						  int qbin,
						  bool singley,
						  bool edge);


 private:
  //Where am I?
  const unsigned int detType_    ;  // 1 for barrel, 0 for forward
  //What type of histogram am I?
  const unsigned int ifEdgeOrBig_; // 0 for regular, 1 for edge, 2 for big
  //Resolution binning
  const double 	cotbetaBinWidth_  ;
  const double 	cotbetaLowEdge_	  ;
  const int	cotbetaBins_	  ;
  const double	cotalphaBinWidth_ ;
  const double	cotalphaLowEdge_  ;
  const int	cotalphaBins_	  ;
  const int	qbinWidth_	  ;
  const int	qbins_  	  ;

  // Resolution histograms.  I (Petar) tried to dynamically allocate
  // these histograms, but all possible implementations were somewhat
  // complicated, which would make the code harder to understand,
  // debug, and thus support in the long term.  Since we are here only
  // booking pointers of histograms, we will instead book larger
  // matrices, and leave them partially empty.  But who cares -- the
  // wasted memory of a few hundred of null pointers is negligible.
  //
  // The constructor will then fill only the first cotbetaBins_ x
  // cotalphaBins_ x qbinBins_ histograms in the matrix, and we'll
  // ignore the rest.
  //
  TH1F *  resMultiPixelXHist_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ][ QBIN_HIST_MAX ] ;
  TH1F * resSinglePixelXHist_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ];
  TH1F *  resMultiPixelYHist_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ][ QBIN_HIST_MAX ];
  TH1F * resSinglePixelYHist_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ];

  // File with histograms to load.
  TFile * file_ ;
  
  // Status of loading.  Check if there were errors.
  int     status_ ;       

  // Identical binning and parameterization for FastSim generators.
  SimpleHistogramGenerator *  resMultiPixelXGen_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ][ QBIN_HIST_MAX ] ;
  SimpleHistogramGenerator * resSinglePixelXGen_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ];
  SimpleHistogramGenerator *  resMultiPixelYGen_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ][ QBIN_HIST_MAX ];
  SimpleHistogramGenerator * resSinglePixelYGen_ [ COTBETA_HIST_MAX ][ COTALPHA_HIST_MAX ];
  
};
#endif
