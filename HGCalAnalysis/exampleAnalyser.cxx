// Analyser that can be used for any pdg id
// compile using, for example: g++ -o exampleAnalyser exampleAnalyser.cxx  `root-config --cflags --libs`

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

using std::string;
using std::vector;

int main( int argc, char *argv[] )
{
  // get file and tree
  string fileName = "infile.root"; // specified sample
  TFile *inputFile = new TFile( fileName.c_str() ); 
  TTree *theTree = (TTree*)inputFile->Get("ana/hgc");

  // define quantities in ntuple
  vector<float> *multiclus_pt = 0;
  theTree->SetBranchAddress("multiclus_pt", &multiclus_pt);

  // define histograms
  TH1* hMultiPt = new TH1F( "hMultiPt", "multicluster pt", 50, 0., 50. ); 

  // loop over events
  uint nEntries = theTree->GetEntries();
  for( uint evtIndex = 0; evtIndex < nEntries; evtIndex++ )
  {
    // get values
    theTree->GetEntry( evtIndex );
     
    // loop over multiclusters
    for( uint multiIndex = 0; multiIndex < multiclus_pt->size(); multiIndex++ ) {
      hMultiPt->Fill( multiclus_pt->at(multiIndex) );
    }
  }

  // draw histograms
  TCanvas *c0 = new TCanvas();
  hMultiPt->Draw("hist");
  c0->Print( "hMultiPt.pdf" );
  c0->Print( "hMultiPt.png" );

  return 0;
}
