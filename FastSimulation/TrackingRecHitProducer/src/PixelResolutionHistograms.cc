#define OUTSIDE_CMSSW_RELEASE 0

#if OUTSIDE_CMSSW_RELEASE
#include "PixelResolutionHistograms.h"

//--- Include a dummy implementation of SimpleHistogramGenerator.
class TH1F;
class SimpleHistogramGenerator {
public:
  SimpleHistogramGenerator(TH1F * hist) : hist_(hist) {};
private:
  TH1F * hist_;
};

#else
//--- The real thing.
#include "FastSimulation/TrackingRecHitProducer/interface/PixelResolutionHistograms.h"
#include "FastSimulation/Utilities/interface/RandomEngineAndDistribution.h"
#include "FastSimulation/Utilities/interface/SimpleHistogramGenerator.h"

#endif


// Generic C stuff
#include <math.h>       /* floor */
#include <iostream>
#include <string>

// ROOT
#include <TFile.h>
#include <TH1F.h>

// Global definitions 
const float cmtomicron = 10000.0;


//------------------------------------------------------------------------------
//  Constructor: Books the FastSim Histograms, given the input parameters
//  which are provided as arguments. These variables are then const inside
//  the class. (That is, once we make the histograms, we can't change the
//  definition of the binning.)
//------------------------------------------------------------------------------
PixelResolutionHistograms::
PixelResolutionHistograms( unsigned int detType, 
			   unsigned int ifEdgeOrBig,
			   double 	cotbetaBinWidth,
			   double 	cotbetaLowEdge,
			   int	cotbetaBins,
			   double	cotalphaBinWidth,
			   double	cotalphaLowEdge,
			   int	cotalphaBins,
			   int	qbinWidth,
			   int	qbins )
  : detType_          ( detType ),
    ifEdgeOrBig_      ( ifEdgeOrBig ),
    cotbetaBinWidth_  ( cotbetaBinWidth ),
    cotbetaLowEdge_   ( cotbetaLowEdge ), 
    cotbetaBins_      ( cotbetaBins ),
    cotalphaBinWidth_ ( cotalphaBinWidth  ),
    cotalphaLowEdge_  ( cotalphaLowEdge ),
    cotalphaBins_     ( cotalphaBins ),
    qbinWidth_        ( qbinWidth ),
    qbins_            ( qbins ),
    resMultiPixelXHist_(), resSinglePixelXHist_(),  // all to nullptr
    resMultiPixelYHist_(), resSinglePixelYHist_(),  // all to nullptr
    file_(), 
    status_(0)
{
  Char_t histo[200];
  Char_t title[200];
  
  TFile output( "PixelBarrelResolution.root", "recreate" );
  //Resolution binning
  // const double 	cotbetaBinWidth = 1.0;
  // const double 	cotbetaLowEdge	= -11.5 ;
  // const int	cotbetaBins	= 23;
  // const double	cotalphaBinWidth	= 0.08 ;
  // const double	cotalphaLowEdge	= -0.36 ;
  // const int	cotalphaBins	= 9;
  // const int	qbinWidth	= 1;
  // const int	qbins		= 4;

  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	//
	sprintf( histo, "hx%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );  //information of bits of histogram names
	//--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha, qbins
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//
	resMultiPixelXHist_ [ ii ][ jj ][ kk ] = new TH1F(histo, title, 1000, -0.05, 0.05);

	sprintf( histo, "hy%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//
	resMultiPixelYHist_ [ ii ][ jj ][ kk ] = new TH1F(histo, title, 1000, -0.05, 0.05);
      }
    }
  }


  for( int ii=0; ii<cotbetaBins_; ii++) {
    for( int jj=0; jj<cotalphaBins_; jj++) {

      sprintf( histo, "hx%d0%02d%d", detType_, ii+1, jj+1 );  //information of bits of histogram names
      //first bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //
      resSinglePixelXHist_ [ ii ][ jj ] = new TH1F(histo, title, 1000, -0.05, 0.05);

      sprintf( histo, "hy%d0%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //
      resSinglePixelYHist_ [ ii ][ jj ] = new TH1F(histo, title, 1000, -0.05, 0.05);

    }
  }

}



//------------------------------------------------------------------------------
//  Another constructor: load the histograms from one file.
//     filename = full path to filename
//     rootdir  = ROOT directory inside the file
//
//  The other parameters are the same (needed later) and must correspond
//  to the histograms we are loading from the file.
//------------------------------------------------------------------------------
PixelResolutionHistograms::
PixelResolutionHistograms( const char * filename, 
			   const char * rootdir,
			   unsigned int detType,
			   unsigned int ifEdgeOrBig,
			   double 	cotbetaBinWidth,
			   double 	cotbetaLowEdge,
			   int	cotbetaBins,
			   double	cotalphaBinWidth,
			   double	cotalphaLowEdge,
			   int	cotalphaBins,
			   int	qbinWidth,
			   int	qbins )
  : detType_          ( detType ),
    ifEdgeOrBig_      ( ifEdgeOrBig  ),
    cotbetaBinWidth_  ( cotbetaBinWidth ),
    cotbetaLowEdge_   ( cotbetaLowEdge ), 
    cotbetaBins_      ( cotbetaBins ),
    cotalphaBinWidth_ ( cotalphaBinWidth  ),
    cotalphaLowEdge_  ( cotalphaLowEdge ),
    cotalphaBins_     ( cotalphaBins ),
    qbinWidth_        ( qbinWidth ),
    qbins_            ( qbins ),
    resMultiPixelXHist_(), resSinglePixelXHist_(),  // all to nullptr
    resMultiPixelYHist_(), resSinglePixelYHist_(),  // all to nullptr
    file_(),
    status_(0)
{
  Char_t histo[200];       // the name of the histogram
  Char_t title[200];       // histo title, for debugging and sanity checking (compare inside file)
  TH1F * tmphist = nullptr;  // cache for histo pointer

  //--- Open the file for reading.
  file_ = new TFile( filename  ,"READ");
  if ( !file_ ) {
    status_ = 1;
    return;
  }
  if((ifEdgeOrBig_ == 0)||(ifEdgeOrBig_ == 2 and detType_ == 1)){ //regular histograms, or barrel big histograms (both have multi and single pixle histos
    //multi pixel histograms
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	//
	sprintf( histo, "hx%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );  //information of bits of histogram names
	//--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha, qbins
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
 	//std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title) 
	//		  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 2;
	  return;
	}
	//std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
	//	  << std::endl;
	resMultiPixelXHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelXGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );

	sprintf( histo, "hy%d1%02d%d%d", detType_, ii+1, jj+1, kk+1 );
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
	//	  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 3;
	  return;
	}
	//std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
	//	  << std::endl;
	resMultiPixelYHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelYGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );
      }
    }
  }
  
  for( int ii=0; ii<cotbetaBins_; ii++) {
    for( int jj=0; jj<cotalphaBins_; jj++) {

      //single pixel histogram
      sprintf( histo, "hx%d0%02d%d", detType_, ii+1, jj+1 );  //information of bits of histogram names
      //--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 X",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //      std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
      //	<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 4;
	  return;
      }
      //      std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
      //		<< std::endl;
      resSinglePixelXHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelXGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );


      sprintf( histo, "hy%d0%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 Y",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
      //		<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 5;
	  return;
      }
      //std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
      //		<< std::endl;
      resSinglePixelYHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelYGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );

    }
  }
  }
  if(ifEdgeOrBig_ == 1){ //edge histograms, they have only single histograms but do have qbin
      for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
      
	sprintf( histo, "hx%d0%02d%d%d", detType_, ii+1, jj+1, kk+1 );  //information of bits of histogram names
	//--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha, qbins
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 X",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title) 
	//	  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 2;
	  return;
	}
	//std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
	//	  << std::endl;
	resMultiPixelXHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelXGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );


	sprintf( histo, "hy%d0%02d%d%d", detType_, ii+1, jj+1, kk+1 );
	sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f qbin %d npixel>1 Y",
		 cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
		 cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_,
		 kk+1 );
	//std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
	//	  << std::endl;
	//
	tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
	if ( !tmphist ) {
	  status_ = 3;
	  return;
	}
	//std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
	//	  << std::endl;
	resMultiPixelYHist_ [ ii ][ jj ][ kk ] = tmphist;
	resMultiPixelYGen_  [ ii ][ jj ][ kk ] = new SimpleHistogramGenerator( tmphist );
      }
    }
      }
  }
  if(ifEdgeOrBig_ == 2){
    if (detType_ == 0){ //big forward histograms, which only have single histograms
  for( int ii=0; ii<cotbetaBins_; ii++) {
    for( int jj=0; jj<cotalphaBins_; jj++) {

      sprintf( histo, "hx%d0%02d%d", detType_, ii+1, jj+1 );  //information of bits of histogram names
      //--- First bit 1/0 barrel/forward, second 1/0 multi/single, cotbeta, cotalpha
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 X",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
      //	<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 4;
	  return;
      }
      //std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
      //		<< std::endl;
      resSinglePixelXHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelXGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );


      sprintf( histo, "hy%d0%02d%d", detType_, ii+1, jj+1 );
      sprintf( title, "cotbeta %.1f-%.1f cotalpha %.2f-%.2f npixel=1 Y",
	       cotbetaLowEdge_ + ii*cotbetaBinWidth_ , cotbetaLowEdge_ + (ii+1)*cotbetaBinWidth_,
	       cotalphaLowEdge_ +jj*cotalphaBinWidth_, cotalphaLowEdge_ +(jj+1)*cotalphaBinWidth_ );
      //std::cout << "Debug: histo=" << std::string(histo) << " title=" << std::string(title)  
      //	<< std::endl;
      //
      tmphist = (TH1F*) file_->Get( Form( "%s%s" , rootdir, histo ) );
      if ( !tmphist ) {
	  status_ = 5;
	  return;
      }
      //std::cout << "Debug: found histo with title = " << std::string( tmphist->GetTitle() ) 
      //	<< std::endl;
      resSinglePixelYHist_ [ ii ][ jj ] = tmphist;
      resSinglePixelYGen_  [ ii ][ jj ] = new SimpleHistogramGenerator( tmphist );

      
    }
  }
    }}
}




//------------------------------------------------------------------------------
//  Destructor.  Use file_ pointer to tell whether we loaded the histograms
//  from a file (and do not own them), or we built them ourselves and thus need
//  to delete them.
//------------------------------------------------------------------------------
PixelResolutionHistograms::~PixelResolutionHistograms()
{
  //--- Delete histograms, but only if we own them. If 
  //--- they came from a file, let them be.
  //
  if ( file_ ) {
    file_->Close();
    delete file_ ;
    file_ = 0;
  }
  else {
    //--- We made the histograms, so delete them all.
    for( int ii=0; ii<cotbetaBins_; ii++ ) {
      for( int jj=0; jj<cotalphaBins_; jj++ ) {
	for( int kk=0; kk<qbins_; kk++ ) {
	  delete resMultiPixelXHist_ [ ii ][ jj ][ kk ];
	  delete resMultiPixelYHist_ [ ii ][ jj ][ kk ];
	}
      }
    }
    for( int ii=0; ii<cotbetaBins_; ii++ ) {
      for( int jj=0; jj<cotalphaBins_; jj++ ) {
	delete resSinglePixelXHist_ [ ii ][ jj ];
	delete resSinglePixelYHist_ [ ii ][ jj ]; 
	//	delete qbinHist_ [ ii ][ jj ];
      }
    }
  } // else



  //--- Delete FastSim generators. (It's safe to delete a nullptr.)
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      for( int kk=0; kk<qbins_; kk++ ) {
	delete resMultiPixelXGen_ [ ii ][ jj ][ kk ];
	delete resMultiPixelYGen_ [ ii ][ jj ][ kk ];
      }
    }
  }
  for( int ii=0; ii<cotbetaBins_; ii++ ) {
    for( int jj=0; jj<cotalphaBins_; jj++ ) {
      delete resSinglePixelXGen_ [ ii ][ jj ];
      delete resSinglePixelYGen_ [ ii ][ jj ]; 
      //   delete qbinGen_ [ ii ][ jj ];
    }
  }

}


  



//------------------------------------------------------------------------------
//  Fills the appropriate FastSim histograms.
//  Returns 0 if the relevant histogram(s) were found and filled, 1 if not.
//------------------------------------------------------------------------------
int
PixelResolutionHistograms::
Fill( double dx, double dy, double cotalpha, double cotbeta, 
      int qbin, int nxpix, int nypix ) 
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  //  std::cout << "Debug (getGeneratorY): "<<icotalpha<<" "<<icotbeta<<" "<<iqbin<< std::endl;
  if( icotalpha >= 0 && icotalpha < cotalphaBins_ && icotbeta >= 0 && icotbeta < cotbetaBins_ ) {
    if( nxpix == 1 )
      resSinglePixelXHist_[icotbeta][icotalpha]->Fill(dx/cmtomicron);
    else
      resMultiPixelXHist_[icotbeta][icotalpha][iqbin]->Fill(dx/cmtomicron);
    if( nypix == 1 )
      resSinglePixelYHist_[icotbeta][icotalpha]->Fill(dy/cmtomicron);
    else
      resMultiPixelYHist_[icotbeta][icotalpha][iqbin]->Fill(dy/cmtomicron);
  }
  
  return 0;  
}



  //--- Get generators, for resolution in X and Y.  Use in FastSim.
const SimpleHistogramGenerator * 
PixelResolutionHistograms::
getGeneratorX( double cotalpha, double cotbeta, int qbin, bool single, bool edge )
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  //handling underflow and overflow
  if ( icotalpha < 0 )
    icotalpha = 0;
  
  if ( icotalpha >= cotalphaBins_ ) 
    icotalpha = cotalphaBins_ - 1;
  if ( icotbeta < 0 )
    icotbeta = 0;

  if ( icotbeta >= cotbetaBins_ ) 
    icotbeta = cotbetaBins_ - 1;
  // std::cout << "Debug (getGeneratorX): "<<icotalpha<<" "<<icotbeta<<" "<<iqbin<< std::endl;

  if( single ){
    if (edge) //edge histograms have cotbeta, cotalpha, and qbin but are labeled as single
      return resMultiPixelXGen_[icotbeta][icotalpha][iqbin];
    else
      return resSinglePixelXGen_[icotbeta][icotalpha];
  }
  else
    return resMultiPixelXGen_[icotbeta][icotalpha][iqbin]; 
}

const SimpleHistogramGenerator * 
PixelResolutionHistograms::
getGeneratorY( double cotalpha, double cotbeta, int qbin, bool single, bool edge )
{
  int icotalpha, icotbeta, iqbin ;
  icotalpha = (int)floor( (cotalpha - cotalphaLowEdge_)/cotalphaBinWidth_ ) ;
  icotbeta  = (int)floor( (cotbeta - cotbetaLowEdge_) /cotbetaBinWidth_ ) ;
  iqbin = qbin > 2 ? 3 : qbin;
  //
  //handling underflow and overflow
  if ( icotalpha < 0 )
    icotalpha = 0;

  if ( icotalpha >= cotalphaBins_ ) 
    icotalpha = cotalphaBins_ - 1;

  if ( icotbeta < 0 )
    icotbeta = 0;

  if ( icotbeta >= cotbetaBins_ ) 
    icotbeta = cotbetaBins_ - 1;
  //std::cout << "Debug (getGeneratorY): "<<icotalpha<<" "<<icotbeta<<" "<<iqbin<< std::endl;

    if( single ){
    if (edge) //edge histograms have cotbeta, cotalpha, and qbin but are labeled as single
	return resMultiPixelYGen_[icotbeta][icotalpha][iqbin];
      else
	return resSinglePixelYGen_[icotbeta][icotalpha];
    }
   else
      return resMultiPixelYGen_[icotbeta][icotalpha][iqbin];
}
