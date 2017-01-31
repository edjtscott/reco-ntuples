// First attempt to open and analyse Clemens' HGC ntuples
// compile using, for example: g++ -o make_hists histplotter.cxx  `root-config --cflags --libs`
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLegend.h"

#include "interface/AEvent.h"
#include "interface/AObData.h"

#include "usefulClasses.h"
#include "HistSlicer.h"
#include "HistsByCategory.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::abs;
using std::vector;
using std::pair;
using std::make_pair;


//nice colour setup from Shane
/*void SetColor(TH1 * obj, int pos, int max)
{
    double modifier(0.15), colorIndex;
    int colour(1);
    double fraction = (double)(pos)/(double)(max-1);

    if( pos > max-1 || pos < 0 || max < 0 ) colour = 1;
    else
    {
        colorIndex = (fraction * (1.0-2.0*modifier) + modifier) * gStyle->GetNumberOfColors();
        colour = gStyle->GetColorPalette(colorIndex);
    }
    obj->SetLineColor(colour);
    obj->SetMarkerColor(colour);
}*/


template<class T>
void SetColor(T * obj, int pos, int max)
{
    double modifier(0.20), colorIndex;
    int colour(1);
    double fraction = (double)(pos)/(double)(max-1);

    if( pos > max-1 || pos < 0 || max < 0 ) colour = 1;
    else
    {
        colorIndex = (fraction * (1.0-2.0*modifier) + modifier) * gStyle->GetNumberOfColors();
        colour = gStyle->GetColorPalette(colorIndex);
    }
    obj->SetLineColor(colour);
    obj->SetMarkerColor(colour);
}

void MultiDrawerThree( TCanvas *c, string name, string title, string conversion,  int legendX=1, int legendY=1, int isNum=0, string outdir="HGCPlots/" ) {
  vector<string> ptValues;
  ptValues.push_back("35");
  ptValues.push_back("10");
  ptValues.push_back("5");
  //ptValues.push_back("2");
  vector< TGraph* > vecGraphs;
  TMultiGraph *multiGraph = new TMultiGraph();
  for( uint i=0; i<ptValues.size(); i++ ) {
    //string tempFileName = "output_Photon_Pt" + ptValues[i] + "_" + conversion + ".root";
    //string tempFileName = "Output/output_Photon_Pt" + ptValues[i] + "_" + conversion + "_Old.root";
    string tempFileName = "Output/output_Photon_Pt" + ptValues[i] + "_" + conversion + "_SensorDependent.root";
    //string tempFileName = "Output/output_Pion_Pt" + ptValues[i] + "_" + conversion + "_SensorDependent.root";
    TFile *tempFile = new TFile( tempFileName.c_str() );
    vecGraphs.push_back( (TGraph*)tempFile->Get( name.c_str() ) );
    tempFile->Close();
    delete tempFile;
    SetColor( vecGraphs[i], i, ptValues.size() );
    vecGraphs[i]->SetLineWidth(1);
    multiGraph->Add( vecGraphs[i] );
  }
  multiGraph->Draw("AC*");
  multiGraph->SetTitle( title.c_str() );
  if( name.find("VsRadius") != string::npos ) multiGraph->GetXaxis()->SetTitle("r / cm");
  else if( name.find("VsCut") != string::npos ) multiGraph->GetXaxis()->SetTitle("Number of required 2D clusters");
  if( isNum==2 ) {
    multiGraph->SetMinimum(0.);
    multiGraph->SetMaximum(3.);
    multiGraph->GetYaxis()->SetTitle("Mean number of selected multiclusters");
    multiGraph->GetYaxis()->SetTitleOffset(1.2);
  }
  else if( isNum==0 ) {
    multiGraph->SetMinimum(0.75);
    multiGraph->SetMaximum(1.25);
    multiGraph->GetYaxis()->SetTitle("Mean fraction of gen energy");
    multiGraph->GetYaxis()->SetTitleOffset(1.2);
  }
  else if( isNum==1 ) {
    multiGraph->SetMinimum(0.);
    multiGraph->SetMaximum(1.2);
    multiGraph->GetYaxis()->SetTitle("Mean fraction of gen energy");
    multiGraph->GetYaxis()->SetTitleOffset(1.2);
  }
  else if( isNum==3 ) {
    multiGraph->SetMinimum(0.);
    multiGraph->SetMaximum(0.15);
    multiGraph->GetYaxis()->SetTitle("Resolution of gen energy fraction");
    multiGraph->GetYaxis()->SetTitleOffset(1.2);
  }
  float legXpos = 0.2 + legendX*0.4;
  float legYpos = 0.2 + legendY*0.4;
  TLegend* l = new TLegend(legXpos,legYpos,legXpos+0.2,legYpos+0.2);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  for( uint i=0; i<vecGraphs.size(); i++ ) {
    string legendText = "p_{T} = " + ptValues[i] + " GeV";
    l->AddEntry( vecGraphs[i], legendText.c_str(), "l" );
  }
  l->Draw();
  string printName = outdir + name + conversion + ".pdf";
  c->Print(printName.c_str());
  printName = outdir + name + conversion + ".png";
  c->Print(printName.c_str());
  delete l;
}


int main( int argc, char *argv[] )
{
  TCanvas *c = new TCanvas("c","c",600,500);
  c->cd();

  gStyle->SetPalette(kBird);

  //void MultiDrawerThree( TCanvas *c, string name, string title, string conversion,  int legendX=1, int legendY=1, int isNum=0, string outdir="HGCPlots/" ) {
  
  /*
  MultiDrawerThree( c, "gNumSelectedVsCut", "Mean number of selected multiclusters vs required number of 2D clusters", "Unconverted", 1, 1, 2 );
  MultiDrawerThree( c, "gSelectedEnergyVsCut", "Mean fraction of gen energy of selected multiclusters vs required number of 2D clusters", "Unconverted", 0, 0, 1 );
  MultiDrawerThree( c, "gRecHitEFracMeanVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the gen photon", "Unconverted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracMeanBestVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the best multicluster", "Unconverted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionBestVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of the best multicluster", "Unconverted", 1, 1, 3 );
  MultiDrawerThree( c, "gRecHitEFracMeanSelectedVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the selected multiclusters", "Unconverted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionSelectedVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of any selected multicluster", "Unconverted", 1, 1, 3 );
  */

  MultiDrawerThree( c, "gNumSelectedVsCut", "Mean number of selected multiclusters vs required number of 2D clusters", "Converted", 1, 1, 2 );
  MultiDrawerThree( c, "gSelectedEnergyVsCut", "Mean fraction of gen energy of selected multiclusters vs required number of 2D clusters", "Converted", 0, 0, 1 );
  MultiDrawerThree( c, "gRecHitEFracMeanVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the gen photon", "Converted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracMeanBestVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the best multicluster", "Converted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionBestVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of the best multicluster", "Converted", 1, 1, 3 );
  MultiDrawerThree( c, "gRecHitEFracMeanSelectedVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the selected multiclusters", "Converted", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionSelectedVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of any selected multicluster", "Converted", 1, 1, 3 );

  /*
  MultiDrawerThree( c, "gNumSelectedVsCut", "Mean number of selected multiclusters vs required number of 2D clusters", "All", 1, 1, 2 );
  MultiDrawerThree( c, "gSelectedEnergyVsCut", "Mean fraction of gen energy of selected multiclusters vs required number of 2D clusters", "All", 0, 0, 1 );
  MultiDrawerThree( c, "gRecHitEFracMeanVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the gen photon", "All", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracMeanBestVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the best multicluster", "All", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionBestVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of the best multicluster", "All", 1, 1, 3 );
  MultiDrawerThree( c, "gRecHitEFracMeanSelectedVsRadius", "Mean fraction of gen energy contained by recHits within r cm of the selected multiclusters", "All", 1, 0 );
  MultiDrawerThree( c, "gRecHitEFracResolutionSelectedVsRadius", "Resolution of gen energy fraction contained by recHits within r cm of any selected multicluster", "All", 1, 1, 3 );
  */
  
  return 0;
}
