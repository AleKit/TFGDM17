/* \file PPP.C
   \brief root macro to make CTA public performance pages

   main input: WP Phys-style root file with IRFs

   \changelog
   2015-01-21  first version and submission to svn

   \author Gernot Maier

 */


#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraphSmooth.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMarker.h"
#include "TROOT.h"
#include "TText.h"

//#include "CTASensitivityRequirements.C"

using namespace std;

class PlotCTAPerformancePlots
{

  private:

    string             fVersionNumber;             // plot a version text (e.g. date)

    bool               bPlotRequirements;          // plot CTA requirements
    bool               bPlotOtherInstruments;      // plot non-CTA performance curves
    bool               bPlotVERITAS;               // plot un-published VERITAS performance curves

    string             fPlotPrintName;             // print canvases to pdf

    double             fEnergy_TeV_min;            // min energy (plotting range) [TeV]
    double             fEnergy_TeV_max;            // min energy (plotting range) [TeV]
    double             fEnergy_TeV_Eres_min;       // min energy for energy resolution [TeV]
    double             fEnergy_TeV_Erec_oneSided;  // set energy from which the one sided (upper) part
    // is used only for the calculation of the energy resolution

    // plotting style
    int                fCTAColor;
    double             fCTALineWidth;
    int                fOtherInstrumentColor;

    int                fHAWC1yColor;
    int                fHAWC5yColor;
    int                fMAGICColor;
    int                fVERITASColor;
    int                fHESSColor;

    // private functions
    TGraphErrors*      getHistogramFromFile( string iIRF_root_file, string iHistogramName, 
        string iTitle,
        bool iGetErrorBarX = false, bool iGetErrorBarY = false );

    TGraphAsymmErrors *readTextFile( string iTextFile );

    TCanvas*           plot_canvas( string iName, string iYaxis, bool iYLog, 
        double iYmin, double iYmax, 
        int iCanvasX = 900, int iCanvasY = 600,
        bool EnergyAxisIsTrueEnergy = false );
    void printCanvas( TCanvas *c, string iTitle = "" );

    void writeASCII_File( TGraphErrors *g, string iTitle, string iFile, double iRebin = -1 );

  public:

    PlotCTAPerformancePlots();
    ~PlotCTAPerformancePlots() {}

    void plot_and_print_AngularResolution( string iIRF_root_file, bool bNewCanvas = true,
        int iLineStyle = 1,
        double x_marker = 0., double y_marker = 0.,
        string iTextMarker = "", bool i_bPrint = false,
        string iASCIIFILE = "" );

    void plot_and_print_Backgroundrate( string iIRF_root_file, bool bNewCanvas = true, 
        int iLineStyle = 1, 
        double x_marker = 0., double y_marker = 0.,
        string iTextMarker = "", int iColor = 1,
        bool i_bPrint = false, bool i_bSquDegree = true,
        string iASCIIFILE = "" );

    void plot_and_print_DifferentialFluxSensitivity( string iIRF_root_file, bool bNewCanvas = true,
        int iLineStyle = 1,
        double x_marker = 0., double y_marker = 0.,
        double iTextAngle = 0., string iTextMarker = "",
        bool i_bPrint = false, string iASCIIFILE = "" );

    void plot_and_print_EffectiveArea( string iIRF_root_file, bool bEffAreaWithoutDirectionCut,
        bool bNewCanvas = true, int iLineStyle = 1,
        double x_marker = 0., double y_marker = 0.,
        string iTextMarker = "", bool i_bPrint = false,
        string iASCIIFILE = "" );

    void plot_and_print_EnergyResolution( string iIRF_root_file, bool bNewCanvas = true, 
        int iLineStyle = 1,
        bool bUseMedian = true,
        double x_marker = 0., double y_marker = 0.,
        string iTextMarker = "", bool i_bPrint = false,
        bool bXaxisIsEtrue = true, string iASCII_File = "" );

    void plot_offAxisSensitivity( bool iSouth = true, string iASCIIFILE = "" );

    void plotPerformanceComparison();

    void plotPerformanceNorth( int iFileType = 1, bool iPlotEnergyResolution = true, bool bXaxisIsEtrue = false );

    void plotPerformanceSouth( int iFileType = 1, bool iPlotEnergyResolution = true, bool bXaxisIsEtrue = false );

    void setCTAPlottingStyle( int iColor = 1, int iLineWidth = 3. )
    {
      fCTAColor = iColor;
      fCTALineWidth = iLineWidth;
    }


    void setEnergyRange( double iE_TeV_min = 0.01, double iE_TeV_max = 300. )
    {
      fEnergy_TeV_min = iE_TeV_min;
      fEnergy_TeV_max = iE_TeV_max;
    }

    void setEnergyForOneSidedIntervalInEnergyReconstruction( double i_TeV_Erec_oneSided = 0.05, double iEnergy_TeV_Eres_min = 0.02 )
    {
      fEnergy_TeV_Erec_oneSided = i_TeV_Erec_oneSided;
      fEnergy_TeV_Eres_min = iEnergy_TeV_Eres_min;
    }

    void setOtherInstrumentColor( int iColor = 14 )                     // color of other instruments
    { 
      fOtherInstrumentColor = iColor;
    }

    void setHAWCColor( int i1y = 14, int i5y = 14 )
    {
      fHAWC1yColor = i1y;
      fHAWC5yColor = i5y;
    }

    void setHESSColor( int iHESS = 14 )
    {
      fHESSColor = iHESS;
    }

    void setVERITASColor( int iVTS = 14 )
    {
      fVERITASColor = iVTS;
    }

    void setMAGICColor( int iMAGIC = 14 )
    {
      fMAGICColor = iMAGIC;
    }

    void setPrintName( string iName = "" )
    {
      fPlotPrintName = iName;
    }

    void setPlotOtherInstruments( bool iOtherInstruments = false )
    {
      bPlotOtherInstruments = iOtherInstruments;
    }

    void setPlotCTARequirements( bool iPlotRequirements = false )
    {
      bPlotRequirements = iPlotRequirements;
    }

    void setPlotVERITAS( bool iVERITAS = false )
    {
      bPlotVERITAS = iVERITAS;
    }

    void setPlotVersionNumber( string iVersionNumber = "www.cta-observatory.org (2015-05-11)" )
    {
      fVersionNumber = iVersionNumber;
    }
};


//////////////////////////////////////////////////////////////////


PlotCTAPerformancePlots::PlotCTAPerformancePlots()
{

  // default setters
  setCTAPlottingStyle();
  setPlotCTARequirements();
  setOtherInstrumentColor();
  setPlotOtherInstruments();
  setPlotVERITAS();
  setEnergyRange();
  setPrintName();
  setPlotVersionNumber();
  setEnergyForOneSidedIntervalInEnergyReconstruction();
}

void PlotCTAPerformancePlots::writeASCII_File( TGraphErrors *i_g, string iTitle, string iFile, double iRebin )
{
  if( !i_g ) return;

  TGraphErrors *g = i_g;
  // rebin to 5 bins per decade
  if( iRebin > 0. )
  {
    g = new TGraphErrors( 1 );
    double x_min = 0.;
    double x_max = 0.;
    double y = 0.;
    i_g->GetPoint( 0, x_min, y );
    i_g->GetPoint( i_g->GetN()-1, x_max, y );
    int z = 0;
    int n_max = (int)(4. / iRebin);
    for( int i = 0; i < n_max; i++ )
    {
      double x = -1.8 + i * iRebin;
      if( x >= x_min && x <= x_max )
      {
        g->SetPoint( z, x, i_g->Eval( x ) );
        z++;
      }
    }
  }

  ofstream ofile;
  ofile.open( iFile.c_str() );

  double x = 0.;
  double xe = 0.;
  double y = 0.;
  double ye = 0.;

  ofile << endl;
  ofile << "CTA Performance - Northern site" << endl;
  ofile << "for further details, see " << fVersionNumber << endl;
  ofile << "=======================================================" << endl;
  if( iTitle.size() > 0 )
  {
    ofile << iTitle << endl;
    ofile << "=======================================================" << endl;
  }
  ofile << endl;
  if( iTitle.find( "Sensitivity" ) != string::npos )
  {
    ofile << "E_min (TeV) \t E_max (TeV) \t S (erg / cm^2 / s)" << endl;
    ofile << "----------- \t------------ \t-----------------" << endl;
    for( int i = 0; i < g->GetN(); i++ )
    {
      g->GetPoint( i, x, y );
      xe = g->GetErrorX( i );
      ye = g->GetErrorY( i );
      ofile << setprecision( 3 );
      ofile << TMath::Power( 10., x-xe ) << "\t\t";
      ofile << TMath::Power( 10., x+xe ) << "\t\t";
      ofile << y;
      ofile << endl;
    }
  }
  else if( iTitle.find( "Effec" ) != string::npos )
  {
    ofile << "E (TeV) \t Effective" << endl;
    ofile << "        \t Area (m^2)" << endl;
    ofile << "------- \t------------------" << endl;
    for( int i = 0; i < g->GetN(); i++ )
    {
      g->GetPoint( i, x, y );
      ofile << setprecision( 3 ); 
      ofile << TMath::Power( 10., x ) << "\t\t ";
      ofile << setprecision( 3 );
      ofile << y;
      ofile << endl;
    }
  }
  else if( iTitle.find( "background" ) != string::npos )
  {
    if( iTitle.find( "deg" ) != string::npos )
    {
      ofile << "E_min (TeV) \t E_max (TeV) \t Bck Rate (Hz/deg^2)" << endl;
    }
    else
    {
      ofile << "E_min (TeV) \t E_max (TeV) \t Bck Rate (Hz)" << endl;
    }
    ofile << "----------- \t------------ \t-----------------" << endl;
    for( int i = 0; i < g->GetN(); i++ )
    {
      g->GetPoint( i, x, y );
      xe = g->GetErrorX( i );
      ye = g->GetErrorY( i );
      ofile << setprecision( 3 );
      ofile << TMath::Power( 10., x-xe ) << "\t\t";
      ofile << TMath::Power( 10., x+xe ) << "\t\t";
      ofile << setprecision( 2 );
      ofile << y << " +- " << ye;
      ofile << endl;
    }
  }
  else
  {
    if( iTitle.find( "Angular" ) != string::npos )
    {
      ofile << "E (TeV) \t Angular" << endl;
      ofile << "        \t resolution (deg)" << endl;
    }
    else
    {
      ofile << "E (TeV) \t Energy" << endl;
      ofile << "        \t resolution" << endl;
    }
    ofile << "------- \t------------------" << endl;
    for( int i = 0; i < g->GetN(); i++ )
    {
      g->GetPoint( i, x, y );
      ofile << setprecision( 3 ); 
      ofile << TMath::Power( 10., x ) << "\t\t ";
      ofile << setprecision( 2 );
      ofile << y;
      ofile << endl;
    }
  }

  ofile.close();
}

void PlotCTAPerformancePlots::printCanvas( TCanvas *c, string iTitle )
{
  if( c && fPlotPrintName.size() > 0 )
  {
    c->Update();

    char htitle[500];
    sprintf( htitle, "%s-%s", fPlotPrintName.c_str(), iTitle.c_str() );
    if( bPlotOtherInstruments )
    {
      sprintf( htitle, "%s-OtherInstruments", htitle );
    }
    sprintf( htitle, "%s.pdf", htitle );
    c->Print( htitle );
  }
}

/*

   read values from a text file, expect the following format:

   EMin_GeV  EMin_GeV  Flux_CU
 */	
TGraphAsymmErrors *PlotCTAPerformancePlots::readTextFile( string iTextFile )
{
  ifstream is; 
  is.open( iTextFile.c_str(), ifstream::in );
  if( !is )
  {   
    return 0;
  }   
  TGraphAsymmErrors *iG = new TGraphAsymmErrors( 1 );
  string iLine = "";
  double x_1 = 0.;
  double x_2 = 0.;
  double y = 0.;
  int z = 0;
  while( is >> x_1 >> x_2 >> y )
  {   
    //         istringstream is_stream( iLine );
    //	 is_stream >> x_1;
    //	 is_stream >> x_2;
    //	 is_stream >> y;
    // fill graph in units of TeV
    iG->SetPoint( z, 0.5*(x_1+x_2)/1.e3, y );
    iG->SetPointEXhigh( z, (x_2 - 0.5*(x_1+x_2) )/1.e3 );
    iG->SetPointEXlow( z, (0.5*(x_1+x_2) - x_1 ) /1.e3 );
    z++;
  }   
  is.close();

  return iG;
}

/*

   read a IRF histogram from the CTA WP Phys root file

 */
TGraphErrors* PlotCTAPerformancePlots::getHistogramFromFile( string iIRF_root_file, string iHistogramName, string iTitle, 
    bool iGetErrorBarX, bool iGetErrorBarY )
{
  TFile *f = new TFile( iIRF_root_file.c_str() );
  if( f->IsZombie() )
  {
    cout << "error opening file " << iIRF_root_file << endl;
    return 0;
  }
  TH1F *h = (TH1F*)f->Get( iHistogramName.c_str() );
  if( h )
  {
    TGraphErrors *g = new TGraphErrors( 1 );
    int z = 0;
    cout << endl;
    cout << "reading from file " << iIRF_root_file << endl;
    cout << "======================================================================" << endl;
    cout << "Energy [TeV] \t" << iTitle << endl;
    cout << "----------------------------------------------------------------------" << endl;
    for( int i = 1; i <= h->GetNbinsX(); i++ )
    {
      if( h->GetBinContent( i ) > 0 )
      {
        g->SetPoint( z, h->GetXaxis()->GetBinCenter( i ), h->GetBinContent( i ) );
        //	       cout << setprecision( 3 ) << TMath::Power( 10., h->GetXaxis()->GetBinCenter( i ) ) << "\t\t";
        //	       cout << h->GetBinContent( i ) << endl;
        double i_x_error = 0.;
        double i_y_error = 0.;
        if( iGetErrorBarX )
        {
          i_x_error = 0.5*h->GetBinWidth( i );
        }
        if( iGetErrorBarY )
        {
          i_y_error = h->GetBinError( i );
        }
        g->SetPointError( z, i_x_error, i_y_error );
        z++;
      }
    }
    g->SetLineWidth( fCTALineWidth );
    return g;
  }
  return 0;
}

/*

   plot a canvas with axis, etc.

   (handles log Energy axis on a log scale

 */
TCanvas* PlotCTAPerformancePlots::plot_canvas( string iName, string iYaxis, bool iYLog,
    double iYmin, double iYmax,
    int iCanvasX, int iCanvasY,
    bool EnergyAxisIsTrueEnergy )
{
  // canvas
  TCanvas *c = new TCanvas( iName.c_str(), iYaxis.c_str(), 10, 10, iCanvasX, iCanvasY );
  c->SetGridx( 0 );
  c->SetGridy( 0 );
  c->SetLeftMargin( 0.15 );
  c->SetTopMargin( 0.08 );
  c->SetBottomMargin( 0.12 );

  // null histogram
  TH1D* hnull = new TH1D( "hnullSens", "", 10, log10(fEnergy_TeV_min), log10(fEnergy_TeV_max) );
  hnull->SetStats( 0 );
  if( EnergyAxisIsTrueEnergy )
  {
    hnull->SetXTitle( "Energy E_{T} (TeV)" );
  }
  else
  {
    hnull->SetXTitle( "Energy E_{R} (TeV)" );
  }
  hnull->SetYTitle( iYaxis.c_str() );
  hnull->SetMinimum( iYmin );
  hnull->SetMaximum( iYmax );
  hnull->Draw( "AH" );
  c->Update();

  // X-(Energy)-axis
  TGaxis* x1 = new TGaxis( c->GetUxmin(), c->GetUymin(), c->GetUxmax(), c->GetUymin(), fEnergy_TeV_min, fEnergy_TeV_max, 510, "G" );
  x1->SetTextFont( 42 );
  x1->SetLabelFont( 42 );
  x1->SetTitle( hnull->GetXaxis()->GetTitle() );
  x1->SetLabelOffset( 0.003 );
  x1->SetTitleOffset( 1.25 );
  x1->Draw();
  // upper ticks for energy axis
  TGaxis* x2 = new TGaxis( c->GetUxmin(), c->GetUymax(), c->GetUxmax(), c->GetUymax(), fEnergy_TeV_min, fEnergy_TeV_max, 510, "-UG" );
  x2->Draw();

  double iYTitleOffset = 1.2;
  if( (double)iCanvasX / (double)iCanvasY < 1.3 )
  {
    iYTitleOffset = 1.3;
  }

  if( iYLog )
  {
    TGaxis* y1 = new TGaxis( c->GetUxmin(), c->GetUymin(), c->GetUxmin(), c->GetUymax(), hnull->GetMinimum(), hnull->GetMaximum(), 510, "G" );
    y1->SetTitle( hnull->GetYaxis()->GetTitle() );
    y1->SetTextFont( 42 );
    y1->SetLabelFont( 42 );
    y1->SetTitleOffset( iYTitleOffset );
    y1->Draw();
    TGaxis* y2 = new TGaxis( c->GetUxmax(), c->GetUymin(), c->GetUxmax(), c->GetUymax(), hnull->GetMinimum(), hnull->GetMaximum(), 510, "+GU" );
    y2->SetTitleOffset( iYTitleOffset );
    y2->Draw();
  }
  else
  {
    TGaxis* y1 = new TGaxis( c->GetUxmin(), c->GetUymin(), c->GetUxmin(), c->GetUymax(), hnull->GetMinimum(), hnull->GetMaximum(), 510, "" );
    y1->SetTextFont( 42 );
    y1->SetLabelFont( 42 );
    y1->SetTitle( hnull->GetYaxis()->GetTitle() );
    y1->SetTitleOffset( iYTitleOffset );
    y1->Draw();

    TGaxis* y2 = new TGaxis( c->GetUxmax(), c->GetUymin(), c->GetUxmax(), c->GetUymax(), hnull->GetMinimum(), hnull->GetMaximum(), 510, "+U" );
    y2->SetTitleOffset( iYTitleOffset );
    y2->Draw();
  }
  if( fVersionNumber.size() > 0 )
  {
    TText *iV = new TText();
    iV->SetNDC( true );
    iV->SetTextFont( 42 );
    iV->SetTextAngle( 90. );
    iV->SetTextSize( iV->GetTextSize() * 0.5 );
    iV->DrawTextNDC( 0.92, 0.2, fVersionNumber.c_str() );
  }

  return c;
}

/*

   plot differential sensitivity


 */
void PlotCTAPerformancePlots::plot_and_print_DifferentialFluxSensitivity( string iIRF_root_file, bool bNewCanvas, 
    int iLineStyle, double x_marker, double y_marker,
    double iTextAngle, string iTextMarker,
    bool i_bPrint, string iASCIIFILE )
{
  TCanvas *c = 0;
  if( bNewCanvas )
  {
    // good for 50/5/0.5 hour plots
    //        c = plot_canvas( "diffFluxSensitivity", "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})", true, 2.5e-14, 2.e-10, 900, 600 );
    // good for 50h comparision (default)
    c = plot_canvas( "diffFluxSensitivity", "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})", true, 2.5e-14, 4.e-11, 900, 600 );
    // good for 5h comparision
    //        c = plot_canvas( "diffFluxSensitivity", "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})", true, 2.5e-13, 4.e-11, 900, 600 );
    c->SetLogy( true );
  }
  else
  {
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "diffFluxSensitivity" );
  }

  // set this to false if we don't want to plot CTA
  bool iPlotCTA = true;
  TGraphErrors *DiffSens = 0;
  if( iPlotCTA )
  {
    DiffSens = getHistogramFromFile( iIRF_root_file, "DiffSens", "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})" );
    if( !DiffSens )
    {
      cout << "error retrieving histogram DiffSens from " << iIRF_root_file << endl;
      return;
    }

    TGraphErrors *DiffSensC = getHistogramFromFile( iIRF_root_file, "DiffSens", "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})", true );
    DiffSensC->SetLineWidth( 1 );
    DiffSensC->SetLineColor( fCTAColor );
    DiffSensC->Draw( "p" );
    DiffSens->SetLineStyle( iLineStyle );
    DiffSens->SetLineColor( fCTAColor );
    DiffSens->Draw( "l" );

    if( iASCIIFILE.size() > 0 )
    {
      writeASCII_File( DiffSensC, "E^{2} x Flux Sensitivity (erg cm^{-2} s^{-1})", iASCIIFILE );
    }


    if( iTextMarker.size() > 0 )
    {
      TText *iL = new TText( x_marker, y_marker, iTextMarker.c_str() );
      iL->SetTextColor( fCTAColor );
      iL->SetTextFont( 42 );
      iL->SetTextSize( iL->GetTextSize()*0.9 );
      iL->SetTextAngle( iTextAngle );
      iL->Draw();
    }

    // differential sensitivity text
    TText *iDF = new TText( -1.7, 3.5e-14, "Differential sensitivity (5 bins per decade in energy)" );
    iDF->SetTextFont( 42 );
    iDF->SetTextSize( iDF->GetTextSize()*0.6);
    iDF->Draw();
  }

  if( bPlotOtherInstruments && bNewCanvas )
  {
    double x = 0.;
    double y = 0.;
    double y_m = 0.;
    double x_bin = 0.;
    int z = 0;
    // MAGIC
    TGraphAsymmErrors *magic_hegraCrab = readTextFile( "Instrument/MAGIC_DiffSensCU.dat" );
    // TGraph *magic_magicCrab = new TGraph( "Instrument/MAGIC_DiffSensCU.dat" );
    TGraphAsymmErrors *magic_magicCrab = readTextFile( "Instrument/MAGIC_DiffSensCU.dat" );
    if( !magic_magicCrab ) return;
    magic_magicCrab->Print();
    for( int i = 0; i < magic_hegraCrab->GetN(); i++ )
    {
      magic_magicCrab->GetPoint( i, x, y );
      double x_low  = x - magic_magicCrab->GetErrorXlow( i );
      double x_high = x + magic_magicCrab->GetErrorXhigh( i ); 
      x_bin = x_high - x_low;
      if( x > 0. )
      {
        // use HEGRA Crab spectrum (Aharonian et al 2004.)
        y_m = y * 2.83e-11 * TMath::Power( x, -2.62 );
        magic_hegraCrab->SetPoint( z, log10( x ), y_m * x * x * 1.e12 * TMath::Qe() / 1.e-7 );
        magic_hegraCrab->SetPointEXhigh( z, 0. );
        magic_hegraCrab->SetPointEXlow( z, 0. );
        cout << "HEGRA Crab at " << x << " TeV: " << y_m/1.e-12 << endl;
        // use MAGIC curve power law
        y_m = y * 3.39e-11 * TMath::Power( x, -2.51-0.21*log10( x ) );
        magic_magicCrab->SetPoint( z, log10( x ), y_m * x * x * 1.e12 * TMath::Qe() / 1.e-7 );
        magic_magicCrab->SetPointEXhigh( z, 0. );
        magic_magicCrab->SetPointEXlow( z, 0. );
        cout << "MAGIC Crab at " << x << " TeV: " << y_m/1.e-12 << endl;
        cout << "Ratio MAGIC/HEGRA at " << x << " TeV: " << (3.39e-11 * TMath::Power( x, -2.51-0.21*log10( x ) )) /  (2.83e-11 * TMath::Power( x, -2.62 )) << endl;
        z++;
      }
    }
    magic_hegraCrab->SetLineColor( fMAGICColor );
    magic_hegraCrab->SetLineStyle( 4 );
    //	magic_hegraCrab->Draw( "c" );

    magic_magicCrab->SetLineColor( fMAGICColor );
    magic_magicCrab->SetLineStyle( 9 );
    magic_magicCrab->Draw( "c" );

    TText *iLMagic = new TText( -0.27, 1.20e-12, "MAGIC 50 h" );
    iLMagic->SetTextColor( 1 );
    iLMagic->SetTextSize( iLMagic->GetTextSize()*0.5 );
    iLMagic->SetTextAngle( -15. );
    iLMagic->SetTextFont( 42 );
    iLMagic->SetTextColor( fMAGICColor );
    iLMagic->Draw();


    ///////////////////////////////////////
    // VERITAS sensitivities
    // (best of curve from 
    // http://veritas.sao.arizona.edu/about-veritas-mainmenu-81/veritas-specifications-mainmenu-111)

    // (in Crab Units)
    TGraph *VTS = new TGraph( "Instrument/VERITAS_V6_std_50hr_5sigma_VERITAS2014_DiffSens.dat" );
    VTS->SetLineColor( fVERITASColor );
    VTS->SetLineStyle( 9 );
    VTS->SetLineWidth( 3 );

    z = 0;
    for( int i = 0; i < VTS->GetN(); i++ )
    {
      VTS->GetPoint( i, x, y );
      // VTS Crab spectrum (unpublished; DO NOT USE FOR ANYTHING)
      y_m = y * 3.269e-11 * TMath::Power( x, -2.474 - 0.191*log10( x ) ); 
      VTS->SetPoint( z, log10( x ), y_m * x * x * 1.e12 * TMath::Qe() / 1.e-7 );
      cout << "VTS Crab at " << x << " TeV: " << y_m/1.e-12 << endl;
      z++;
    }
    VTS->Draw( "c" );
    TText *iLVTS = new TText( 0.9, 9.e-13, "VERITAS 50 h" );
    iLVTS->SetTextColor( 1 );
    iLVTS->SetTextSize( iLVTS->GetTextSize()*0.55 );
    iLVTS->SetTextAngle( 35. );
    iLVTS->SetTextFont( 42 );
    iLVTS->SetTextColor( fVERITASColor );
    iLVTS->Draw();

    ///////////////////////////////////////
    // HESS sensitivities (50h)
    // (ICRC 2015)
    vector< TGraph * > HESS;
    HESS.push_back( new TGraph( "Instrument/HESS_August2015_CT15_Combined_Std.dat" ) );
    HESS.push_back( new TGraph( "Instrument/HESS_August2015_CT14_Std.dat" ) );
    HESS.push_back( new TGraph( "Instrument/HESS_August2015_CT15_Stereo_Std.dat" ) );
    HESS.push_back( new TGraph( "Instrument/HESS_August2015_CT5_Mono_Std.dat" ) );
    for( unsigned int t = 0; t < HESS.size(); t++ )
    {
      if( !HESS[t] )
      {
        continue;
      }
      HESS[t]->SetLineColor( fHESSColor );
      HESS[t]->SetLineStyle( 4 );
      HESS[t]->SetLineWidth( 3 );

      z = 0;
      for( int i = 0; i < HESS[t]->GetN(); i++ )
      {
        HESS[t]->GetPoint( i, x, y );
        if( y < 0 ) 
        {
          continue;
        }
        x = TMath::Power( 10., x );
        HESS[t]->SetPoint( z, log10( x ), y * x * x * 1.e12 * TMath::Qe() / 1.e-7 );
        z++;
      }

      HESS[t]->Draw( "c" );
    }
    TText *iLHESS = new TText( -1.10, 2.1e-11, "H.E.S.S. 50 h" );
    iLHESS->SetTextColor( 1 );
    iLHESS->SetTextSize( iLHESS->GetTextSize()*0.5 );
    iLHESS->SetTextAngle( 285. );
    iLHESS->SetTextFont( 42 );
    iLHESS->SetTextColor( fHESSColor );
    iLHESS->Draw();

    ///////////////////////////////////////
    // HAWC sensitivities (1y and 5 y)
    // Note: calculated per quarter decade
    // from Fig 4 in Abeysekara et al 2013 (astro-ph/1306.5800)
    TGraph *HAWC5y = new TGraph( "Instrument/HAWC300_5y_QuarterDecade_DiffSens.dat" );
    HAWC5y->SetLineColor( fHAWC5yColor );
    HAWC5y->SetLineStyle( 5 );
    HAWC5y->SetLineWidth( 3 );
    for( int i = 0; i < HAWC5y->GetN(); i++ )
    {
      HAWC5y->GetPoint( i, x, y );
      HAWC5y->SetPoint( i, log10( x ) - 3., y );
    }
    TGraph *HAWC1y = new TGraph( "Instrument/HAWC300_1y_QuarterDecade_DiffSens.dat" );
    HAWC1y->SetLineColor( fHAWC1yColor );
    HAWC1y->SetLineStyle( 6 );
    HAWC1y->SetLineWidth( 3 );
    for( int i = 0; i < HAWC1y->GetN(); i++ )
    {
      HAWC1y->GetPoint( i, x, y );
      HAWC1y->SetPoint( i, log10( x ) - 3., y );
    }
    HAWC5y->Draw( "c" );
    HAWC1y->Draw( "c" );
    TText *iLHAWC1y = new TText( 1.3, 2.98e-12, "HAWC 1 year" );
    iLHAWC1y->SetTextSize( iLHAWC1y->GetTextSize()*0.5 );
    iLHAWC1y->SetTextAngle( 40. );
    iLHAWC1y->SetTextFont( 42 );
    iLHAWC1y->SetTextColor( fHAWC1yColor );
    iLHAWC1y->Draw();
    TText *iLHAWC5y = new TText( 1.3, 1.10e-12, "HAWC 5 year" );
    iLHAWC5y->SetTextSize( iLHAWC5y->GetTextSize()*0.5 );
    iLHAWC5y->SetTextAngle( 40. );
    iLHAWC5y->SetTextFont( 42 );
    iLHAWC5y->SetTextColor( fHAWC5yColor );
    iLHAWC5y->Draw();
  }
  ///////////////////////////////////////
  // plot sensitivitiy requirements

  if( bPlotRequirements && DiffSens )
  {
    CTASensitivityRequirements diffsens_req;
    TGraph *diffsens_req_graph = new TGraph( 1 );
    TGraph *diffsens_goal_graph = new TGraph( 1 );
    double x = 0.;
    double y = 0.;
    for( int i = 0; i < DiffSens->GetN(); i++ )
    {
      DiffSens->GetPoint( i, x, y );
      if( iIRF_root_file.find( "_50h" ) != string::npos )
      {
        if( iIRF_root_file.find( "Aar" ) != string::npos )
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req50_E2erg_south( TMath::Power( 10., x ) ) );
          diffsens_goal_graph->SetPoint( i, x, diffsens_req.Flux_goal50_E2erg_south( TMath::Power( 10., x ) ) );
        }
        else
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req50_E2erg_north( TMath::Power( 10., x ) ) );
          diffsens_goal_graph->SetPoint( i, x, diffsens_req.Flux_goal50_E2erg_north( TMath::Power( 10., x ) ) );
        }
      }
      else if( iIRF_root_file.find( "_5h" ) != string::npos )
      {
        if( iIRF_root_file.find( "Aar" ) != string::npos )
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req5_E2erg_south( TMath::Power( 10., x ) ) );
        }
        else
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req5_E2erg_north( TMath::Power( 10., x ) ) );
        }
      }
      else if( iIRF_root_file.find( "_0.5h" ) != string::npos )
      {
        if( iIRF_root_file.find( "Aar" ) != string::npos )
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req05_E2erg_south( TMath::Power( 10., x ) ) );
        }
        else
        {
          diffsens_req_graph->SetPoint( i, x, diffsens_req.Flux_req05_E2erg_north( TMath::Power( 10., x ) ) );
        }
      }
    }
    diffsens_req_graph->SetLineColor( 2 );
    diffsens_req_graph->SetLineStyle( iLineStyle );
    diffsens_req_graph->Draw( "l" );
    if( diffsens_goal_graph->GetN() > 1 )
    {
      diffsens_goal_graph->SetLineColor( 3 );
      diffsens_goal_graph->SetLineStyle( iLineStyle );
      diffsens_goal_graph->Draw( "l" );
    }	 
  }


  if( i_bPrint ) printCanvas( c, "DifferentialSensitivity" );
}

/*

   68\% angular resolution

 */
void PlotCTAPerformancePlots::plot_and_print_AngularResolution( string iIRF_root_file, bool bNewCanvas,
    int iLineStyle,
    double x_marker, double y_marker,
    string iTextMarker, bool i_bPrint, string iASCIIFile )
{
  TCanvas *c = 0;
  if( bNewCanvas )
  {
    c = plot_canvas( "AngResolution", "Angular Resolution (#circ)", false, 0.00, 0.25, 700, 600, false );
    // TMP        c = plot_canvas( "AngResolution", "Angular Resolution (#circ)", false, 0.00, 0.18, 700, 600, false );
  }
  else
  {
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "AngResolution" );
  }


  TGraphErrors *AngRes = getHistogramFromFile( iIRF_root_file, "AngRes", "Angular Resolution (deg)" );
  if( !AngRes )
  {
    cout << "error retrieving histogram AngRes from " << iIRF_root_file << endl;
    return;
  }

  AngRes->SetLineStyle( iLineStyle );
  AngRes->SetLineColor( fCTAColor );
  AngRes->Draw( "c" );

  if( iASCIIFile.size() > 0 )
  {
    writeASCII_File( AngRes, "Angular Resolution (68\% containment radius)", iASCIIFile );
  }

  TText *iL = new TText( x_marker, y_marker, iTextMarker.c_str() );
  iL->SetTextColor( fCTAColor );
  iL->SetTextFont( 42 );
  iL->SetTextSize( iL->GetTextSize()*0.9 );
  iL->SetTextAngle( 0. );
  iL->Draw();

  if( bPlotOtherInstruments )
  {
    double x = 0.;
    double y = 0.;
    int z = 0;
    //////////////////////////////////////////////////////////
    // HAWC

    // based on optimal bin
    TGraph *hawc = new TGraph( "Instrument/HAWC_AngRes68.dat" );
    for( int i = 0; i < hawc->GetN(); i++ )
    {
      hawc->GetPoint( i, x, y );
      if( x > 0. )
      {
        hawc->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    hawc->SetLineColor( fHAWC5yColor );
    hawc->SetLineStyle( 4 );
    hawc->Draw( "c" );
    // based on 2D Gaussian
    z = 0;
    TGraph *hawc2D = new TGraph( "Instrument/HAWC_AngRes2DGaussian.dat" );
    for( int i = 0; i < hawc->GetN(); i++ )
    {
      hawc2D->GetPoint( i, x, y );
      if( x > 0. )
      {
        hawc2D->SetPoint( z, log10( x ), y*1.52 );  // conversion factor from 2D sigma to 68% value
        z++;
      }
    }
    hawc2D->SetLineColor( 4 );
    hawc2D->SetLineStyle( 4 );
    hawc2D->Draw( "c" );


    TText *iLhawc = new TText( 1.6, 0.125, "HAWC" );
    iLhawc->SetTextSize( iLhawc->GetTextSize()*0.5 );
    iLhawc->SetTextColor( fHAWC5yColor );
    iLhawc->SetTextFont( 42 );
    iLhawc->Draw();


    //////////////////////////////////////////////////////////
    // LAT Pass 8
    TGraph *latp8 = new TGraph( "Instrument/LAT_pass8_AngRes68.dat" );
    z = 0;
    for( int i = 0; i < latp8->GetN(); i++ )
    {
      latp8->GetPoint( i, x, y );
      if( x > 0. )
      {
        latp8->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    latp8->SetLineStyle( 5 );
    latp8->SetLineColor( fOtherInstrumentColor );
    latp8->Draw( "c" );
    TText *iLlatp8 = new TText( -1.9, 0.08, "Fermi LAT Pass 8" );
    iLlatp8->SetTextSize( iLlatp8->GetTextSize()*0.5 );
    iLlatp8->SetTextColor( fOtherInstrumentColor );
    iLlatp8->SetTextFont( 42 );
    iLlatp8->Draw();
    //////////////////////////////////////////////////////////
    // MAGIC
    TGraph *magic = new TGraph( "Instrument/MAGIC_AngRes68.dat" );
    z = 0;
    for( int i = 0; i < magic->GetN(); i++ )
    {
      magic->GetPoint( i, x, y );
      if( x > 0. )
      {
        magic->SetPoint( z, log10( x/1.e3 ), y );
        z++;
      }
    }
    magic->SetLineStyle( 5 );
    magic->SetLineColor( fMAGICColor );
    magic->Draw( "c" );
    TText *iLMagic = new TText( -0.93, 0.154, "MAGIC" );
    iLMagic->SetTextSize( iLMagic->GetTextSize()*0.5 );
    iLMagic->SetTextColor( fMAGICColor );
    iLMagic->SetTextFont( 42 );
    iLMagic->Draw();
    //////////////////////////////////////////////////////////
    // VERITAS
    TGraph *VERITAS = new TGraph( "Instrument/VERITAS_V6_AngRes68.dat" );
    z = 0;
    for( int i = 0; i < VERITAS->GetN(); i++ )
    {
      VERITAS->GetPoint( i, x, y );
      if( x > 0. )
      {
        VERITAS->SetPoint( z, log10( x/1.e3 ), y );
        z++;
      }
    }
    VERITAS->SetLineStyle( 9 );
    VERITAS->SetLineColor( fVERITASColor );
    VERITAS->Draw( "c" );
    TText *iLVERITAS = new TText( -0.56, 0.11, "VERITAS" );
    iLVERITAS->SetTextSize( iLVERITAS->GetTextSize()*0.5 );
    iLVERITAS->SetTextColor( fVERITASColor );
    iLVERITAS->SetTextFont( 42 );
    iLVERITAS->Draw();
  }
  if( bPlotRequirements )
  {
    ////////////////////////////////////////////////////////////
    // requirements
    // from MAN-PO/121004, version 2.5, June 14, 2013

    // from Figure 4
    TGraph *cta_required = new TGraph( "Instrument/CTA_RequiredAngRes.dat" );
    double x = 0.;
    double y = 0.;
    int z = 0; 
    for( int i = 0; i < cta_required->GetN(); i++ )
    {
      cta_required->GetPoint( i, x, y );
      if( x > 0. )
      {
        cta_required->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    cta_required->SetLineColor( 2 );
    cta_required->SetLineStyle( 4 );
    cta_required->Draw( "c" );
    TGraph *cta_goal = new TGraph( "Instrument/CTA_GoalAngRes.dat" );
    z = 0; 
    for( int i = 0; i < cta_goal->GetN(); i++ )
    {
      cta_goal->GetPoint( i, x, y );
      if( x > 0. )
      {
        cta_goal->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    cta_goal->SetLineColor( 3 );
    cta_goal->SetLineStyle( 4 );
    cta_goal->Draw( "c" );
  }

  if( i_bPrint ) printCanvas( c, "AngularResolution" );
}


void PlotCTAPerformancePlots::plot_and_print_EffectiveArea( string iIRF_root_file, bool bEffAreaWithoutDirectionCut,
    bool bNewCanvas, int iLineStyle,
    double x_marker, double y_marker,
    string iTextMarker, bool i_bPrint, string iASCIIFile )
{
  char hname[200];
  char htitle[200];

  TCanvas *c = 0;
  if( bNewCanvas )
  {
    if( bEffAreaWithoutDirectionCut )
    {
      sprintf( hname, "EffectiveAreaNoDirectionCut" );
    }
    else
    {
      sprintf( hname, "EffectiveArea" );
    }
    sprintf( htitle, "effective area (m^{2})" );
    c = plot_canvas( hname, htitle, true, 1, 5.e7, 700, 600, true );
    c->SetLogy( true );
  }
  else
  {
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "EffectiveArea" );
  }

  TGraphErrors *EffArea = 0;
  if( bEffAreaWithoutDirectionCut )
  {
    EffArea = getHistogramFromFile( iIRF_root_file, "EffectiveAreaEtrueNoTheta2cut", "effective area (m^2)" );
  }
  else
  {
    EffArea = getHistogramFromFile( iIRF_root_file, "EffectiveAreaEtrue", "effective area (m^2)" );
  }
  if( !EffArea )
  {
    cout << "error retrieving histogram EffArea from " << iIRF_root_file << endl;
    return;
  }

  EffArea->SetLineStyle( iLineStyle );
  EffArea->SetLineColor( fCTAColor );
  EffArea->Draw( "c" );

  if( iASCIIFile.size() > 0 )
  {
    if( bEffAreaWithoutDirectionCut )
    {
      writeASCII_File( EffArea, "Effective area (after gamma/hadron separation cuts)", iASCIIFile, 0.05 );
    }
    else
    {
      writeASCII_File( EffArea, "Effective area (after gamma/hadron separation and direction cuts)", iASCIIFile, 0.05 );
    }
  }

  TLegend *iL = new TLegend( x_marker, y_marker, x_marker + 0.50, y_marker + 0.28 );
  iL->SetLineColor( 0 );
  iL->SetTextFont( 42 );
  iL->SetFillStyle( 0 );
  iL->SetBorderSize( 0. );
  iL->SetTextSize( 0.045 );
  iL->AddEntry( EffArea, iTextMarker.c_str(), "l" );
  iL->Draw();

  // effecitve area text
  TText *iDF = 0;
  if( bEffAreaWithoutDirectionCut )
  {
    iDF =  new TText( -2.0, 1.e7, "Effective area (after gamma/hadron separation cuts)" );
  }
  else
  {
    iDF =  new TText( -2.0, 1.e7, "Effective area (after gamma/hadron separation and direction cuts)" );
  }
  iDF->SetTextFont( 42 );
  iDF->SetTextSize( iDF->GetTextSize()*0.6);
  iDF->Draw();

  if( bPlotOtherInstruments )
  {
    // MAGIC
    TGraph *magic = new TGraph( "Instrument/MAGIC_EffectiveArea.dat" );
    magic->SetLineStyle( 5 );
    magic->SetLineColor( fMAGICColor );
    magic->Draw( "c" );
    TText *iLMagic = new TText( 1., 196918., "MAGIC" );
    iLMagic->SetTextColor( 1 );
    iLMagic->SetTextSize( iLMagic->GetTextSize()*0.5 );
    iLMagic->SetTextFont( 42 );
    iLMagic->SetTextColor( fMAGICColor );
    iLMagic->Draw();
  }

  ////////////////////////////////////////////////////////////
  // requirements
  // from MAN-PO/121004, version 2.5, June 14, 2013
  // section 6.1.5
  if( bPlotRequirements )
  {
    TGraph *goal_effectiveArea_m2 = new TGraph( 4 );
    // South
    if( iIRF_root_file.find( "Aar" ) != string::npos )
    {
      goal_effectiveArea_m2->SetPoint( 0, log10(   0.03 ), 1.e4 );
      goal_effectiveArea_m2->SetPoint( 1, log10(   0.10 ), 6.e4 );
      goal_effectiveArea_m2->SetPoint( 2, log10(   1.00 ), 1.e6 );
      goal_effectiveArea_m2->SetPoint( 3, log10(  10.03 ), 4.e6 );
      goal_effectiveArea_m2->SetPoint( 4, log10( 100.00 ), 7.e6 );
    }
    else
    {
      goal_effectiveArea_m2->SetPoint( 0, log10(   0.03 ), 1.e4 );
      goal_effectiveArea_m2->SetPoint( 1, log10(   0.10 ), 6.e4 );
      goal_effectiveArea_m2->SetPoint( 2, log10(   1.00 ), 0.5e6 );
      goal_effectiveArea_m2->SetPoint( 3, log10(  10.03 ), 1.e6 );
    }
    goal_effectiveArea_m2->SetMarkerStyle( 23 );
    goal_effectiveArea_m2->SetMarkerColor( 3 );
    goal_effectiveArea_m2->Draw( "p" );
  }


  if( i_bPrint ) printCanvas( c, c->GetName() );
}


void PlotCTAPerformancePlots::plot_and_print_Backgroundrate( string iIRF_root_file, bool bNewCanvas,
    int iLineStyle, double x_marker, double y_marker, 
    string iTextMarker, int iColor,
    bool i_bPrint, bool i_bSquDegree, string iASCIIFILE )
{
  TCanvas *c = 0;
  if( bNewCanvas )
  {
    if( i_bSquDegree )
    {
      c = plot_canvas( "BGRateSquDeg", "background rate (Hz/deg^{2})", true, 1.e-6, 12., 700, 600, false );
    }
    else
    {
      c = plot_canvas( "BGRate", "background rate (Hz)", true, 5.e-8, 8.e-1, 700, 600, false );
    }
    c->SetLogy( true );
  }
  else
  {
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "BGRate" );
  }

  TGraphErrors *BGRate = 0;
  if( i_bSquDegree )
  {
    BGRate = getHistogramFromFile( iIRF_root_file, "BGRatePerSqDeg", "background rate (Hz/deg^{2})", true, true );
  }
  else
  {
    BGRate = getHistogramFromFile( iIRF_root_file, "BGRate", "background rate (Hz)", true, true );
  }
  if( !BGRate )
  {
    cout << "error retrieving histogram BGRate from " << iIRF_root_file << endl;
    return;
  }

  BGRate->SetLineWidth( 2 );
  BGRate->SetMarkerStyle( iLineStyle );
  BGRate->SetLineColor( iColor );
  BGRate->SetMarkerColor( iColor );
  BGRate->Draw( "p" );

  if( iASCIIFILE.size() > 0 )
  {
    if( i_bSquDegree )
    {
      writeASCII_File( BGRate, "background rates vs energy (Hz/deg^2)", iASCIIFILE );
    }
    else
    {
      writeASCII_File( BGRate, "background rates vs energy (Hz)", iASCIIFILE );
    }
  }

  TLegend *iL = new TLegend( x_marker, y_marker, x_marker + 0.50, y_marker + 0.28 );
  iL->SetLineColor( 0 );
  iL->SetTextFont( 42 );
  iL->SetFillStyle( 0 );
  iL->SetBorderSize( 0. );
  iL->SetTextSize( 0.040 );
  iL->AddEntry( BGRate, iTextMarker.c_str(), "pl" );
  iL->Draw();

  if( i_bPrint )
  {
    if( i_bSquDegree ) printCanvas( c, "BackgroundRateSquDeg" );
    else               printCanvas( c, "BackgroundRate" );
  }

}

void PlotCTAPerformancePlots::plot_and_print_EnergyResolution( string iIRF_root_file, bool bNewCanvas, 
    int iLineStyle, bool bUseMedian,
    double x_marker, double y_marker,
    string iTextMarker, bool i_bPrint, bool bXaxisIsEtrue,
    string iASCII_File )
{
  // read migration matrix from file
  TFile *f = new TFile( iIRF_root_file.c_str() );
  if( f->IsZombie() )
  {
    return;
  }

  TH2F *migmatrix = (TH2F*)f->Get( "MigMatrix" );
  if( !migmatrix )
  {
    return;
  }
  // plot migration matrix
  TCanvas *c_mig = plot_canvas( "migrationMatrix", "energy_{MC} (TeV)", false, log10( 0.008 ), log10( 200. ), 700, 600, false );
  migmatrix->SetStats( 0 );
  migmatrix->Draw( "same col" );

  //////////////////////////////////////////////////////////////////////////////////
  // calculate energy resolution as the 68% interval around the most probably value 
  // for the energy reconstruction

  // canvas
  TCanvas *c = 0;
  if( bNewCanvas )
  {
    c = plot_canvas( "energyResolution", "#Delta E/E (68\% containment)", false, 0, 0.35, 700, 600, bXaxisIsEtrue );
  }
  else
  {
    c = ( TCanvas* )gROOT->GetListOfCanvases()->FindObject( "energyResolution" );
  }
  char hname[600];

  // canvas to plot some selected distributions
  TCanvas *cD = new TCanvas( "cD", "energy dispersion", 200, 200, 700, 700 );
  cD->Divide( 3, 3 );
  TCanvas *cC = new TCanvas( "cC", "cum energy dispersion", 250, 200, 700, 700 );
  cC->Divide( 3, 3 );

  vector< int > iBinsToPlot;
  iBinsToPlot.push_back( 30 );
  iBinsToPlot.push_back( 40 );
  iBinsToPlot.push_back( 50 );
  iBinsToPlot.push_back( 60 );
  iBinsToPlot.push_back( 70 ); 
  iBinsToPlot.push_back( 80 );
  iBinsToPlot.push_back( 90 );
  iBinsToPlot.push_back( 100 );
  iBinsToPlot.push_back( 250 ); 
  //    iBinsToPlot.push_back( 311 ); 
  //    iBinsToPlot.push_back( 400 ); 

  TGraph *Eres = new TGraph( 1 );
  Eres->SetMarkerStyle( 20 );
  int z = 0;
  int i_z = 0;
  int zC = 1;
  float eres = 0.;
  float eres_temp = 0.;
  double energy = 0.;
  double energy_mop_med = 0.;

  TF1 *fGauss = new TF1( "fGaus", "[0]*exp(-0.5*((x-[1])/[2])**2)", -2., 2. );
  TH1F *hDispersion = new TH1F( "hDispersion", "dispersion", 100, -2., 2. );

  int nbins = 0;
  if( bXaxisIsEtrue )
  {
    nbins = migmatrix->GetNbinsY();
    cout << "Number of bins in true energy: " << nbins << endl;
  }
  else
  {
    nbins = migmatrix->GetNbinsX();
    cout << "Number of bins in reconstructed energy: " << nbins << endl;
  }

  ///////////////////////////////////////////////////////////////////////
  // loop over axis in energy of migration matrix
  for( int i = 1; i < nbins; i++ )
  {
    TH1F *hDist = 0;
    // fixed E_true, energy resolution from E_rec distribution
    if( bXaxisIsEtrue )
    {
      hDist = (TH1F*)migmatrix->ProjectionX( "hx", i, i );
      energy = migmatrix->GetYaxis()->GetBinCenter( i );
    }
    // fixed E_rec, energy resolution from E_true distribution
    else
    {
      hDist = (TH1F*)migmatrix->ProjectionY( "hy", i, i );
      energy = migmatrix->GetXaxis()->GetBinCenter( i );
    }
    if( fEnergy_TeV_Eres_min > 0. && energy < log10( fEnergy_TeV_Eres_min ) )
    {
      continue;
    }

    // require at least a few entries
    if( hDist->GetEntries() < 1 ) 
    {
      continue;
    }

    ///////////////////////////////////
    // Fermi LAT method
    // half-width of +-34% intervall around most probably energy

    // most probable energy
    double energy_mprob = hDist->GetXaxis()->GetBinCenter(  hDist->GetMaximumBin() );

    // cummulative distribution
    sprintf( hname, "%s_CUMU", hDist->GetName() );
    TH1F *hcum = ( TH1F* )hDist->Clone( hname );
    hcum->Reset();

    hcum->SetBinContent( 1, hDist->GetBinContent( 1 ) );
    for( int j = 2; j <= hDist->GetNbinsX(); j++ )
    {
      hcum->SetBinContent( j, hDist->GetBinContent( j ) + hcum->GetBinContent( j - 1 ) );
    }
    if( hcum->GetMaximum() > 0. )
    {
      hcum->Scale( 1. / hcum->GetMaximum() );
    }
    double energy_median = hcum->GetXaxis()->GetBinCenter( hcum->FindFirstBinAbove( 0.5 ) );

    // use median or mprob?
    if( bUseMedian )
    {
      energy_mop_med = energy_median;
    }
    else
    {
      energy_mop_med = energy_mprob;
    }

    // determine energy resolution:
    // find 34% interval around energy value
    float prob_mprob = hcum->GetBinContent( hcum->GetXaxis()->FindBin( energy_mop_med ) );
    float prob_mprob_minus = prob_mprob - 0.34;
    if( prob_mprob_minus < 0. ) prob_mprob_minus = 0.;
    float prob_mprob_plus  = prob_mprob + 0.34;
    if( prob_mprob_plus > 1. ) prob_mprob_plus = 1.-1.e-3;
    int prob_minus = hcum->FindFirstBinAbove( prob_mprob_minus );
    int prob_plus  = hcum->FindFirstBinAbove( prob_mprob_plus );
    eres = TMath::Power( 10., hcum->GetXaxis()->GetBinCenter( prob_plus ) ) - TMath::Power( 10., hcum->GetXaxis()->GetBinCenter( prob_minus ) );
    eres /= TMath::Power( 10., energy_mop_med );
    eres *= 0.5;  // half width
    // use one sided interval in the threshold area
    if( fEnergy_TeV_Erec_oneSided > 0. && energy < log10( fEnergy_TeV_Erec_oneSided ) )
    {
      cout << "Double sided Erec at " << energy << " TeV: " << eres * 100. << "%" << endl;
      eres = TMath::Power( 10., hcum->GetXaxis()->GetBinCenter( prob_plus ) ) - TMath::Power( 10., energy_mop_med );
      eres /= TMath::Power( 10., energy_mop_med );
      cout << "One sided Erec at " << energy << " TeV: " << eres * 100. << "%" << endl;
    }

    eres_temp += eres;
    i_z++;

    // average of 20 bins
    if( i_z == 20 )
    {
      Eres->SetPoint( z, energy, eres_temp / (float)i_z );
      z++;
      eres_temp = 0.;
      i_z = 0;
    }

    // check if we want to plot this distribution
    for( unsigned p = 0; p < iBinsToPlot.size(); p++ )
    {
      if( iBinsToPlot[p] == i )
      {
        cout << "Plotting " << i << endl;
        cout << "\t energy: " << energy << endl;
        cout << "\t mprob energy: " << energy_mprob << " (prob " << prob_mprob << ")" << endl;
        cout << "\t median energy: " << energy_median << endl;
        cout << "\t +34% " << prob_plus << "\t energy " << hcum->GetXaxis()->GetBinCenter( prob_plus ) << endl;
        cout << "\t -34% " << prob_minus << "\t energy " << hcum->GetXaxis()->GetBinCenter( prob_minus ) << endl;

        cD->cd( zC );
        sprintf( hname, "energy = %.3f TeV: res. %.1f%%", TMath::Power( 10., energy ), eres*100. );
        hDist->SetStats( 0 );
        hDist->SetTitle( hname );
        // get min and max filled values
        double b_xmin = 0.;
        double b_xmax = 0.;
        for( int b = 1; b <= hDist->GetNbinsX(); b++ )
        {
          if( hDist->GetBinContent( b ) > 0 )
          {
            b_xmin = hDist->GetXaxis()->GetBinCenter( b );
            break;
          }
        }
        for( int b = hDist->GetNbinsX(); b > 0; b-- )
        {
          if( hDist->GetBinContent( b ) > 0 )
          {
            b_xmax = hDist->GetXaxis()->GetBinCenter( b );
            break;
          }
        }
        hDist->SetAxisRange( b_xmin, b_xmax );
        hDist->DrawCopy();
        TBox *iB = new TBox( hcum->GetXaxis()->GetBinCenter( prob_minus ), 0., hcum->GetXaxis()->GetBinCenter( prob_plus ), hDist->GetMaximum() );
        iB->SetFillColor( 18 );
        iB->Draw();
        hDist->DrawCopy( "same" );
        TLine *iLL = new TLine( energy_mop_med, 0., energy_mop_med, hDist->GetMaximum() );
        iLL->SetLineStyle( 2 );
        iLL->Draw();
        cD->Update();

        cC->cd( zC );
        hcum->SetTitle( hname );
        hcum->DrawCopy();
        iB->Draw();
        iLL->Draw();
        hcum->DrawCopy( "same" );
        cC->Update();

        zC++;
      }
    }
    c->cd();

  }
  c->cd();
  Eres->SetLineStyle( iLineStyle );
  Eres->SetLineColor( fCTAColor );
  Eres->SetLineWidth( 3. );
  //   Eres->Draw( "l" );

  // smooth this graph using a Kernel of width 0.4
  TGraphSmooth *Eres_Smooth = new TGraphSmooth( "s" );
  TGraph *Eres_Smooth_graph = Eres_Smooth->SmoothKern( Eres, "normal", 0.4, 200 );

  Eres_Smooth_graph->SetLineStyle( iLineStyle );
  Eres_Smooth_graph->SetLineColor( fCTAColor );
  Eres_Smooth_graph->SetLineWidth( 3. );
  //   Eres_Smooth_graph->SetLineColor(2);
  Eres_Smooth_graph->Draw( "l" );

  if( iASCII_File.size() > 0 )
  {
    writeASCII_File( (TGraphErrors*)Eres_Smooth_graph, "", iASCII_File, 0.2 );
  }

  TText *iL = new TText( x_marker, y_marker, iTextMarker.c_str() );
  iL->SetTextFont( 42 );
  iL->SetTextColor( fCTAColor );
  iL->SetTextSize( iL->GetTextSize()*0.9 );
  iL->SetTextAngle( 0. );
  iL->Draw();

  if( bPlotOtherInstruments )
  {
    double x = 0.;
    double y = 0.;
    // MAGIC
    TGraph *magic = new TGraph( "Instrument/MAGIC_Eres.dat" );
    z = 0;
    for( int i = 0; i < magic->GetN(); i++ )
    {
      magic->GetPoint( i, x, y );
      if( x > 0. )
      {
        magic->SetPoint( z, log10( x/1.e3 ), y );
        z++;
      }
    }
    magic->SetLineStyle( 5 );
    magic->SetLineColor( fMAGICColor );
    magic->Draw( "c" );
    TText *iLMagic = new TText( -0.96, 0.20, "MAGIC" );
    iLMagic->SetTextSize( iLMagic->GetTextSize()*0.5 );
    iLMagic->SetTextColor( fMAGICColor );
    iLMagic->SetTextFont( 42 );
    iLMagic->Draw();
    // latp8
    TGraph *latp8 = new TGraph( "Instrument/LAT_pass8_Eres.dat" );
    z = 0;
    for( int i = 0; i < latp8->GetN(); i++ )
    {
      latp8->GetPoint( i, x, y );
      if( x > 0. )
      {
        latp8->SetPoint( z, log10( x/1.e6 ), y );
        z++;
      }
    }
    latp8->SetLineStyle( 5 );
    latp8->SetLineColor( fOtherInstrumentColor );
    latp8->Draw( "c" );
    TText *iLlatp8 = new TText( -1.9, 0.060, "Fermi LAT Pass 8" );
    iLlatp8->SetTextSize( iLlatp8->GetTextSize()*0.5 );
    iLlatp8->SetTextColor( fOtherInstrumentColor );
    iLlatp8->SetTextFont( 42 );
    iLlatp8->Draw();
  }
  if( bPlotRequirements )
  {
    ////////////////////////////////////////////////////////////
    // requirements
    // from MAN-PO/121004, version 2.5, June 14, 2013

    // from Figure 5
    TGraph *cta_required = new TGraph( "Instrument/CTA_RequiredERes.dat" );
    double x = 0.;
    double y = 0.;
    int z = 0; 
    for( int i = 0; i < cta_required->GetN(); i++ )
    {
      cta_required->GetPoint( i, x, y );
      if( x > 0. )
      {
        cta_required->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    cta_required->SetLineColor( 2 );
    cta_required->SetLineStyle( 4 );
    cta_required->Draw( "l" );
    TGraph *cta_goal = new TGraph( "Instrument/CTA_GoalERes.dat" );
    z = 0; 
    for( int i = 0; i < cta_goal->GetN(); i++ )
    {
      cta_goal->GetPoint( i, x, y );
      if( x > 0. )
      {
        cta_goal->SetPoint( z, log10( x ), y );
        z++;
      }
    }
    cta_goal->SetLineColor( 3 );
    cta_goal->SetLineStyle( 4 );
    cta_goal->Draw( "l" );
  }

  if( i_bPrint ) printCanvas( c, "EnergyResolution" );

}

void PlotCTAPerformancePlots::plotPerformanceNorth( int iFileType, bool iPlotEnergyResolution, bool bXaxisIsEtrue )
{
  string i_50h = "North/TDR_Performance_2Q_Aar_50h.root";
  string i_05h = "North/TDR_Performance_2Q_Aar_5h.root";
  string i_30m = "North/TDR_Performance_2Q_Aar_0.5h.root";
  // IFAE files
  if( iFileType == 1 )
  {
    i_50h = "North/Performance_2N_TEN_50h_20150511.root";
    i_05h = "North/Performance_2N_TEN_5h_20150511.root";
    i_30m = "North/Performance_2N_TEN_0.5h_20150511.root";
  }

  // differentical sensitivity
  plot_and_print_DifferentialFluxSensitivity(i_50h, true, 1, 1.10, 9.5e-14, 40., "CTA North 50 h", false, "CTA-Performance-North-50h-DiffSens.txt" );
  plot_and_print_DifferentialFluxSensitivity(i_05h, false, 2, 1.10, 1.0e-12, 40., "CTA North 5 h", false, "CTA-Performance-North-5h-DiffSens.txt"  );
  plot_and_print_DifferentialFluxSensitivity(i_30m, false, 3, 1.10, 1.0e-11, 40., "CTA North 0.5 h", true, "CTA-Performance-North-0.5-DiffSens.txt" );

  // angular resolution
  plot_and_print_AngularResolution( i_50h, true, 1, -1.36, 0.22, "CTA North", true, "CTA-Performance-North-Angres.txt" );

  // background rates
  plot_and_print_Backgroundrate( i_50h, true, 24, 0.2, 0.17, "CTA North 50 h", fCTAColor, false, false, "CTA-Performance-North-50h-Background.txt");
  plot_and_print_Backgroundrate(i_05h, false, 25, 0.2, 0.11,  "CTA North 5 h", 9, false, false, "CTA-Performance-North-5h-Background.txt");
  plot_and_print_Backgroundrate(i_30m, false, 26, 0.2, 0.05, "CTA North 0.5 h", 8, true, false, "CTA-Performance-North-0.5h-Background.txt"  );

  // background rates / deg^2
  plot_and_print_Backgroundrate( i_50h, true, 24, 0.2, 0.17, "CTA North 50 h", fCTAColor, true, true, "CTA-Performance-North-50h-BackgroundSqdeg.txt" );
  //TMP
  plot_and_print_Backgroundrate( i_05h, true, 24, 0.2, 0.17, "CTA North 5 h", fCTAColor, true, true, "CTA-Performance-North-5h-BackgroundSqdeg.txt" );
  plot_and_print_Backgroundrate( i_30m, true, 24, 0.2, 0.17, "CTA North 0.5 h", fCTAColor, true, true, "CTA-Performance-North-0.5-BackgroundSqdeg.txt" );
  //TMP

  // effective area (after direction cuts)
  setEnergyRange( 0.005, 300. );
  plot_and_print_EffectiveArea( i_50h, false, true, 1, 0.4, 0.34,   "CTA North 50 h ", true, "CTA-Performance-North-50h-EffArea.txt"  );
  plot_and_print_EffectiveArea( i_05h, false, false, 2, 0.4, 0.27,   "CTA North 5 h  ", true, "CTA-Performance-North-5h-EffArea.txt" );
  plot_and_print_EffectiveArea( i_30m, false, false, 3, 0.4, 0.20, "CTA North 0.5 h", true, "CTA-Performance-North-0.5h-EffArea.txt" );

  // effective area (no direction cuts)
  plot_and_print_EffectiveArea( i_50h, true, true, 1, 0.4, 0.34,   "CTA North 50 h ", true, "CTA-Performance-North-50h-EffAreaNoDirectionCut.txt" );
  //TMP
  plot_and_print_EffectiveArea( i_05h, true, true, 1, 0.4, 0.34,   "CTA North 5 h ", true, "CTA-Performance-North-5h-EffAreaNoDirectionCut.txt" );
  plot_and_print_EffectiveArea( i_30m, true, true, 1, 0.4, 0.34,   "CTA North 0.5 h ", true, "CTA-Performance-North-0.5h-EffAreaNoDirectionCut.txt" );
  //TMP

  setEnergyRange();

  // energy resolution
  if( iPlotEnergyResolution )
  {
    plot_and_print_EnergyResolution( i_50h, true, 1, true, -1.2, 0.27, "CTA North", true, bXaxisIsEtrue, "CTA-Performance-North-Eres.txt" );
  } 


}

void PlotCTAPerformancePlots::plotPerformanceComparison()
{
  setPlotOtherInstruments(true);
  setHAWCColor( 41, 45  );
  setMAGICColor( 2 );
  setVERITASColor( 8 );
  setHESSColor( 9 );
  string i_50h_South = "South/Performance_2Q_Aar_50h_20150505.root";
  string i_50h_North = "North/Performance_2N_TEN_50h_20150511.root"; 
  string i_5h_South = "South/Performance_2Q_Aar_5h_20150505.root";
  string i_5h_North = "North/Performance_2N_TEN_5h_20150511.root"; 

  // differential sensitivity
  plot_and_print_DifferentialFluxSensitivity( i_50h_South, true, 1, 1.24, 4.3e-14, 42., "CTA South 50 h", false );
  setCTAPlottingStyle( 4 );
  plot_and_print_DifferentialFluxSensitivity( i_50h_North, false, 2, 1.24, 1.8e-13, 42., "CTA North 50 h", false );

  // angular resolution
  setCTAPlottingStyle( 1 );
  plot_and_print_AngularResolution( i_50h_South, true, 1, -1.36, 0.22, "CTA South", true );
  setCTAPlottingStyle( 4 );
  //    plot_and_print_AngularResolution( i_50h_North, false, 1, -1.25, 0.19, "CTA North", true );
}    


/*

   plot all performance plots for CTA South

 */
void PlotCTAPerformancePlots::plotPerformanceSouth( int iFileType, bool iPlotEnergyResolution, bool bXaxisIsEtrue )
{
  string i_50h = "South/TDR_Performance_2Q_Aar_50h.root";
  string i_05h = "South/TDR_Performance_2Q_Aar_5h.root";
  string i_30m = "South/TDR_Performance_2Q_Aar_0.5h.root";
  // IFAE files
  if( iFileType == 1 )
  {
    i_50h = "South/Performance_2Q_Aar_50h_20150505.root";
    i_05h = "South/Performance_2Q_Aar_5h_20150505.root";
    i_30m = "South/Performance_2Q_Aar_0.5h_20150505.root";
  }

  // differentical sensitivity
  plot_and_print_DifferentialFluxSensitivity(i_50h, true, 1, 1.10, 9.5e-14, 40., "CTA South 50 h", false, "CTA-Performance-South-50h-DiffSens.txt" );
  plot_and_print_DifferentialFluxSensitivity(i_05h, false, 2, 1.10, 1.0e-12, 40., "CTA South 5 h", false, "CTA-Performance-South-5h-DiffSens.txt" );
  plot_and_print_DifferentialFluxSensitivity(i_30m, false, 3, 1.10, 1.0e-11, 40., "CTA South 0.5 h", true, "CTA-Performance-South-0.5h-DiffSens.txt" );

  // angular resolution
  plot_and_print_AngularResolution( i_50h, true, 1, -1.36, 0.22, "CTA South", true, "CTA-Performance-South-Angres.txt" );

  // background rates
  plot_and_print_Backgroundrate( i_50h, true, 24, 0.2, 0.17, "CTA South 50 h", fCTAColor, false, false, "CTA-Performance-South-50h-Background.txt" );
  plot_and_print_Backgroundrate(i_05h, false, 25, 0.2, 0.11,  "CTA South 5 h", 9, false, false, "CTA-Performance-South-5h-Background.txt" );
  plot_and_print_Backgroundrate(i_30m, false, 26, 0.2, 0.05, "CTA South 0.5 h", 8, true, false, "CTA-Performance-South-0.5h-Background.txt" );

  // background rates / deg^2
  plot_and_print_Backgroundrate( i_50h, true, 24, 0.2, 0.17, "CTA South 50 h", fCTAColor, true, true, "CTA-Performance-South-50h-BackgroundSqdeg.txt" );
  // TMP
  plot_and_print_Backgroundrate( i_05h, true, 24, 0.2, 0.17, "CTA South 5 h", fCTAColor, true, true, "CTA-Performance-South-5h-BackgroundSqdeg.txt" );
  plot_and_print_Backgroundrate( i_30m, true, 24, 0.2, 0.17, "CTA South 0.5 h", fCTAColor, true, true, "CTA-Performance-South-0.5h-BackgroundSqdeg.txt" );
  // TMP

  // effective area (after direction cuts)
  setEnergyRange( 0.005, 300. );
  plot_and_print_EffectiveArea( i_50h, false, true, 1, 0.4, 0.34,   "CTA South 50 h ", true, "CTA-Performance-South-50h-EffArea.txt" );
  plot_and_print_EffectiveArea( i_05h, false, false, 2, 0.4, 0.27,   "CTA South 5 h  ", true, "CTA-Performance-South-5h-EffArea.txt" );
  plot_and_print_EffectiveArea( i_30m, false, false, 3, 0.4, 0.20, "CTA South 0.5 h", true, "CTA-Performance-South-0.5h-EffArea.txt" );

  // effective area (no direction cuts)
  plot_and_print_EffectiveArea( i_50h, true, true, 1, 0.4, 0.34,   "CTA South 50 h ", true, "CTA-Performance-South-50h-EffAreaNoDirectionCut.txt" );
  // TMP
  plot_and_print_EffectiveArea( i_05h, true, false, 1, 0.4, 0.34,   "CTA South 5 h ", true, "CTA-Performance-South-5h-EffAreaNoDirectionCut.txt" );
  plot_and_print_EffectiveArea( i_30m, true, false, 1, 0.4, 0.34,   "CTA South 0.5 h ", true, "CTA-Performance-South-0.5h-EffAreaNoDirectionCut.txt" );
  // TMP

  setEnergyRange();

  // energy resolution
  if( iPlotEnergyResolution )
  {
    plot_and_print_EnergyResolution( i_50h, true, 1, false, -1.2, 0.27, "CTA South", true, bXaxisIsEtrue, "CTA-Performance-South-Eres.txt"  );
  } 
}

void PlotCTAPerformancePlots::plot_offAxisSensitivity( bool iSouth, string iASCII_File )
{
  string iFile = "South/OffAxis_Aar_2Q.root";
  if( !iSouth )
  {
    iFile = "North/OffAxis_TEN_2N.root";
  }
  TFile *f = new TFile( iFile.c_str() );
  if( f->IsZombie() )
  {
    cout << "plotOffAxisSensitivity error: file not found: " << iFile << endl;
    return;
  }
  TCanvas *c1 = (TCanvas*)f->Get( "c1" );
  if( !c1 )
  {
    cout << "plotOffAxisSensitivity error: canvas not found " << endl;
    f->ls();
    return;
  }

  c1->SetGridx( 0 );
  c1->SetGridy( 0 );
  c1->SetLeftMargin( 0.125 );

  c1->Draw();
  c1->cd();

  string iText = "CTA South";
  if( !iSouth )
  {
    iText = "CTA North";
  }

  TText *iL = new TText( 0.5, 8.9, iText.c_str() );
  iL->SetTextColor( 1 );
  iL->SetTextFont( 62 );
  iL->SetTextSize( 0.05 );
  iL->Draw();

  TLine *iL1 = new TLine( 0., 1., 4.5, 1. );
  iL1->SetLineStyle( 2 );
  iL1->Draw();

  if( bPlotRequirements )
  {
    TLine *iR2 = new TLine( 0., 2., 3., 2. );
    iR2->SetLineStyle( 3 );
    iR2->Draw();
    TLine *iR3 = new TLine( 1., 1., 1., 2. );
    iR3->SetLineStyle( 3 );
    iR3->Draw();
    TLine *iR4 = new TLine( 3., 1., 3., 2. );
    iR4->SetLineStyle( 3 );
    iR4->Draw();
  }

  if( fVersionNumber.size() > 0 )
  {
    TText *iV = new TText();
    iV->SetNDC( true );
    iV->SetTextFont( 42 );
    iV->SetTextAngle( 90. );
    iV->SetTextSize( iV->GetTextSize() * 0.5 );
    iV->DrawTextNDC( 0.98, 0.2, fVersionNumber.c_str() );
  }

  if( iASCII_File.size() > 0 )
  {
    // get graphs from root file
    vector< TGraph* > i_offsetGraph;
    char hname[200];
    for( int i = 0; i < 4; i++ )
    {
      sprintf( hname, "gr_0%d", i+1 );
      i_offsetGraph.push_back( (TGraph*)c1->GetListOfPrimitives()->FindObject( hname ) );
    }
    // angles to print
    vector< double > i_ang;
    i_ang.push_back( 0.5 );
    i_ang.push_back( 1.0 );
    i_ang.push_back( 1.5 );
    i_ang.push_back( 2.5 );
    i_ang.push_back( 3.5 );
    i_ang.push_back( 4.5 );

    // write ascii file
    ofstream ofile;
    ofile.open( iASCII_File.c_str() );

    double x = 0.;
    double y = 0.;

    ofile << endl;
    ofile << "CTA Performance - Northern site" << endl;
    ofile << "for further details, see " << fVersionNumber << endl;
    ofile << "=======================================================" << endl;
    ofile << "Off-axis sensitivity" << endl;
    ofile << "(point-source sensitivity relative to FoV center)" << endl;
    ofile << "=======================================================" << endl;
    ofile << endl;
    ofile << "Angle w.r. to       50-80 GeV    0.5 - 0.8 TeV     5-8 TeV      50 - 80 TeV" << endl; 
    ofile << "the FoV center   " << endl;
    ofile << "(deg)                      " << endl;               
    ofile << "---------------   -------------  -------------  -------------  -------------" << endl;
    TGraph *g = 0;
    if( i_offsetGraph.size() > 0 && i_offsetGraph[i_offsetGraph.size()-1] )
    {
      g = i_offsetGraph[i_offsetGraph.size()-1];
    }

    for( unsigned int i = 0; i < i_ang.size(); i++ )
    {
      ofile << i_ang[i] << "                      ";
      ofile << setprecision( 3 );
      for( unsigned int j = 0; j < i_offsetGraph.size(); j++ )
      {
        if( i_offsetGraph[j] )
        {
          ofile << i_offsetGraph[j]->Eval( i_ang[i] ) << "          ";
        }
      }

      ofile << endl;
    }
    ofile.close();



  }
}

