// Analyser that can be used for any pdg id
// compile using, for example: g++ -o flexibleAnalyser flexibleAnalyser.cxx  `root-config --cflags --libs`
//
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TStopwatch.h"

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
using std::sort;


int main( int argc, char *argv[] )
{
  if( argc < 4 ) {
    cout << "ERROR - need at least three arguments" << endl;
    cout << "1. Sample, e.g. Photon_Pt35" << endl;
    cout << "2. Conversion, e.g. Unconverted" << endl;
    cout << "3. Extension, e.g. Old" << endl;
    cout << "4. (OPTIONAL) Boolean for doing selected multicluster plots" << endl;
    cout << "5. (OPTIONAL) Layer beyond which to ignore objects, e.g. 28 for photons" << endl;
    exit(1);
  }

  const bool debug       = false;
  const bool drawAndSave = true;
  const bool write       = true;

  const float etaMin     = 1.7;
  const float etaMax     = 2.7;
  const float dvzCut     = 300.;

  const int numLayers    = 52;
  
  cout << "I'm in main, hello" << endl;

  for( unsigned argIndex=0; argIndex < argc; argIndex++ ) cout << "Argument " << argIndex << " is " << argv[argIndex] << endl;
  string sampleType = string( argv[1] );
  string convType   = string( argv[2] );
  string extType    = string( argv[3] );
  string doSelectedMultisStr = "";
  if( argc >= 5 ) doSelectedMultisStr = string( argv[4] );
  int layerCut = 9999;
  if( argc >= 6 ) layerCut = stoi( string(argv[5]) );
  if( layerCut != 9999 ) cout << "layerCut has been set to " << layerCut << endl;
  bool testConversion = true;
  bool converted = false;
  if( convType == "Converted" ) converted = true;
  else if( convType != "Unconverted" ) testConversion = false;
  bool doSelectedMultis = false;
  if( doSelectedMultisStr == "True" || doSelectedMultisStr == "true" ) doSelectedMultis = true;
  if( doSelectedMultis ) cout << "doSelectedMultis has been set to true" << endl;

  string fileName = "/afs/cern.ch/work/e/escott/public/HGCStudies/Ntuples/partGun_" + sampleType + "_" + extType + ".root"; // specified sample
  TFile *inputFile = new TFile( fileName.c_str() ); 
  if( !inputFile )
  {
    cout << "Wrong file address - no .root file opened" << endl;
    return 0;
  } 
  else cout << "Got file" << endl; // this step is not actully working as designed - inputFile still returns true when I don't find a file. Maybe set equal to 0/NULL first?

  TTree *theTree = (TTree*)inputFile->Get("ana/hgc");
  if( !theTree )
  {
    cout << "Wrong tree address - no vaild tree in existence" << endl;
    return 0;
  }
  else cout << "Got tree" << endl;

  // define collections
  //AEvent                    *myEvent; \\ not working, don't know why
  AGenPartCollection        *genParticles  = 0;
  ARecHitCollection         *recHits       = 0;
  ACluster2dCollection      *clusters2D    = 0;
  AMultiClusterCollection   *multiClusters = 0;
  ASimClusterCollection     *simClusters   = 0;
  APFClusterCollection      *pfClusters    = 0;
  ACaloParticleCollection   *caloParticles = 0;

  // set addresses of collections appropriately
  //theTree->SetBranchAddress("event", myEvent); // not working, don't know why
  theTree->SetBranchAddress("particles",     &genParticles);
  theTree->SetBranchAddress("rechits",       &recHits);
  //theTree->SetBranchAddress("rechits_raw",       &recHits);
  theTree->SetBranchAddress("cluster2d",     &clusters2D);
  theTree->SetBranchAddress("multicluster",  &multiClusters);
  theTree->SetBranchAddress("simcluster",    &simClusters);
  theTree->SetBranchAddress("pfcluster",     &pfClusters);
  theTree->SetBranchAddress("caloparticles", &caloParticles);

  // define histograms
  TH1* hGenEta          = new TH1F( "hGenEta", "gen particle eta", 50, -5., 5. ); 
  TH1* hGenPhi          = new TH1F( "hGenPhi", "gen particle phi", 50, -3.14159, 3.14159 ); 
  TH1* hGenEnergy       = new TH1F( "hGenEnergy", "gen particle energy", 50, 0., 500. ); 
  TH1* hGenPt           = new TH1F( "hGenPt", "gen particle pt", 50, 0., 100. ); 
  TH1* hGenUnconv       = new TH1F( "hGenUnconv", "number of unconverted photons vs eta", 10, 1.7, 2.7 ); 
  TH1* hGenConv         = new TH1F( "hGenConv", "number of converted photons vs eta", 10, 1.7, 2.7 );

  TH1* hNumWithExtraVsEta = new TH1F( "hNumWithExtraVsEta", "number of events with more than one selected multicluster vs eta", 5, 1.7, 2.7 ); 
  TH1* hNumWithExtraContigVsEta = new TH1F( "hNumWithExtraContigVsEta", "number of events with more than one multicluster passing contiguity cut vs eta", 5, 1.7, 2.7 ); 
  TH1* hNumVsEta          = new TH1F( "hNumVsEta", "number of events vs eta", 5, 1.7, 2.7 );
  TH1* hDistanceBetweenSelected = new TH1F( "hDistanceBetweenSelected", "distance between consecutive 2D clusters in selected", 40, 0, 8. );
  TH1* hDistanceContigBest = new TH1F( "hDistanceContigBest", "distance between multis passing contiguity and best", 30, 0, 15. );
  TH1* hSelectedBestDeta = new TH1F( "hSelectedBestDeta", "deta between selected multis and \"best\" multi", 50, -0.2, 0.2 );
  TH1* hSelectedBestDphi = new TH1F( "hSelectedBestDphi", "dphi between selected multis and \"best\" multi", 50, -0.2, 0.2 );
  TH1* hSelectedBestDr = new TH1F( "hSelectedBestDr", "dr between selected multis and \"best\" multi", 50, 0, 1.0 );
  TH1* hMultisPassingContigPtFrac = new TH1F( "hMultisPassingContigPtFrac", "fraction of gen pt contained by the multiclusters passing contiguity (but not best)", 100, 0., 2. );

  //TH1* hSingleMCEFrac   = new TH1F( "hSingleMCEFrac",  "fraction of gen energy contained by the \"best\" multicluster", 100, 0., 2. );
  TH1* hSingleMCEFrac   = new TH1F( "hSingleMCEFrac",  "fraction of gen energy contained by the \"best\" multicluster", 40, 0.9, 1.3 );
  TH1* hSingleMCPtFrac  = new TH1F( "hSingleMCPtFrac", "fraction of gen pt contained by the \"best\" multicluster", 100, 0., 2. );
  TH2* hSingleMCPtFracVsEta  = new TH2F( "hSingleMCPtFracVsEta", "fraction of gen pt contained by the \"best\" multicluster vs Eta", 10, 1.7, 2.7, 100, 0., 2. );
  TH1* hSingleMCGenDr   = new TH1F( "hSingleMCGenDr", "dr between \"best\" multicluster and gen", 50, 0, 0.1 );
  TH1* hSingleMCGenDeta = new TH1F( "hSingleMCGenDeta", "delta eta of best multicluster relative to gen", 50, -0.02, 0.02 );
  TH1* hSingleMCGenDphi = new TH1F( "hSingleMCGenDphi", "delta phi of best multicluster relative to gen", 50, -0.02, 0.02 );
  TH2* hSingleMCGenDetaDphi = new TH2F( "hSingleMCGenDetaDphi", "delta eta delta phi plane of best multicluster relative to gen", 50, -0.02, 0.02, 50, -0.02, 0.02 );
  TH2* hSingleMCGenDetaDphiWeighted = new TH2F( "hSingleMCGenDetaDphiWeighted", "delta eta delta phi plane of best multicluster relative to gen, weighted by energy", 
                                             50, -0.02, 0.02, 50, -0.02, 0.02 );

  TH1* hSingleMCLimitedEFrac   = new TH1F( "hSingleMCLimitedEFrac",  "fraction of gen energy contained by 2D clusters in the EE and in the \"best\" multicluster", 100, 0., 2. );
  TH1* hTotalClusterEFrac = new TH1F( "hTotalClusterEFrac",  "fraction of gen energy contained by all 2D clusters", 100, 0., 2. );
  TH1* hPseudoSuperClusterEFrac = new TH1F( "hPseudoSuperClusterEFrac",  "fraction of gen energy contained by multiclusters within deltaR of 0.1 of best and with more than 1\% of best energy", 100, 0., 2. );
  TH1* hBHclusterEta = new TH1F( "hBHclusterEta",  "eta of clusters in BH", 50, -5., 5. );
  TH1* hBHclusterPhi = new TH1F( "hBHclusterPhi",  "phi of clusters in BH", 50, -5., 5. );
  TH1* hBHrecHitEta = new TH1F( "hBHrecHitEta",  "eta of recHits in BH", 50, -5., 5. );
  TH1* hBHrecHitPhi = new TH1F( "hBHrecHitPhi",  "phi of recHits in BH", 50, -5., 5. );
  TH1* hBestPfEFrac = new TH1F( "hBestPfEFrac",  "fraction of gen energy contained by best pf cluster", 100, 0., 2. );

  TH1* hDoubleMCEFrac   = new TH1F( "hDoubleMCEFrac",  "fraction of gen energy contained by the two \"best\" multiclusters", 100, 0., 2. );
  TH1* hDoubleMCPtFrac  = new TH1F( "hDoubleMCPtFrac", "fraction of gen pt contained by the two \"best\" multiclusters", 100, 0., 2. );
  TH1* hSecondMCEFrac   = new TH1F( "hSecondMCEFrac",  "fraction of gen energy contained by the two \"best\" multiclusters", 100, 0., 2. );
  TH1* hSecondMCPtFrac  = new TH1F( "hSecondMCPtFrac", "fraction of gen pt contained by the two \"best\" multiclusters", 100, 0., 2. );
  TH1* hSecondRelEta    = new TH1F( "hSecondRelEta", "eta second MC relative to best", 50, -0.2, 0.2 );
  TH1* hSecondRelPhi    = new TH1F( "hSecondRelPhi", "phi second MC relative to best", 50, -0.2, 0.2 );
  TH1* hSecondRelDr     = new TH1F( "hSecondRelDr", "dr second MC relative to best",  50, 0., 0.5 );
  TH1* hNumMultiClustersHighPt = new TH1F( "hNumMultiClustersHighPt", "number of multiclusters with pt > 0.8 times nominal", 5, -0.5, 4.5 );
  TH1* hNumSelectedMultis = new TH1F( "hNumSelectedMultis", "number of selected multiclusters", 11, -0.5, 10.5 );
  TH1* hNumSelectedMultisContiguous = new TH1F( "hNumSelectedMultisContiguous", "number of selected multiclusters", 11, -0.5, 10.5 );
  TH1* hEFracContiguous = new TH1F( "hEFracContiguous",  "fraction of gen energy contained by the multiclusters passing contiguity selection", 100, 0., 2. );
  TH1* hLayerOfClustersInSelected = new TH1F( "hLayerOfClustersInSelected", "layer of 2D clusters in selected but not best multicluster", numLayers+1, -0.5, numLayers+0.5 );
  TH1* hLayerOfClustersInSelectedWeighted = new TH1F( "hLayerOfClustersInSelectedWeighted", "layer of 2D clusters in selected but not best multicluster, weighted by 2D energy", numLayers+1, -0.5, numLayers+0.5 );
  TH1* hLayerOfClustersInBest = new TH1F( "hLayerOfClustersInBest", "layer of 2D clusters in best multicluster", numLayers+1, -0.5, numLayers+0.5 );
  TH1* hLayerOfClustersInBestWeighted = new TH1F( "hLayerOfClustersInBestWeighted", "layer of 2D clusters in best multicluster, weighted by 2D energy", numLayers+1, -0.5, numLayers+0.5 );
  TH2* hMultiClusterDeltaEtaDeltaPhi = new TH2F( "hMultiClusterDeltaEtaDeltaPhi", "delta eta delta phi plane of multiclusters relative to \"best\" one", 50, -3.14, 3.14, 50, -1.5, 1.5 );
  TH2* hMultiClusterDeltaEtaDeltaPhiWeighted = new TH2F( "hMultiClusterDeltaEtaDeltaPhiWeighted", "delta eta delta phi plane of multiclusters relative to \"best\" one, weighted by energy", 
                                                         50, -3.14, 3.14, 50, -1.5, 1.5 );
  TH2* hMultiClusterDeltaEtaDeltaPhiZoomed = new TH2F( "hMultiClusterDeltaEtaDeltaPhiZoomed", "delta eta delta phi plane of multiclusters relative to \"best\" one",
                                                       50, -0.5, 0.5, 50, -0.5, 0.5 );
  TH2* hMultiClusterDeltaEtaDeltaPhiWeightedZoomed = new TH2F( "hMultiClusterDeltaEtaDeltaPhiWeightedZoomed", "delta eta delta phi plane of multiclusters relative to \"best\" one," 
                                                               "weighted by energy", 50, -0.5, 0.5, 50, -0.5, 0.5 );
  TH2* hMultiClusterDeltaEtaDeltaPhiZoomZoom = new TH2F( "hMultiClusterDeltaEtaDeltaPhiZoomZoom", "delta eta delta phi plane of multiclusters relative to \"best\" one",
                                                       50, -0.2, 0.2, 50, -0.2, 0.2 );
  TH2* hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom = new TH2F( "hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom", "delta eta delta phi plane of multiclusters relative to \"best\" one," 
                                                               "weighted by energy", 50, -0.2, 0.2, 50, -0.2, 0.2 );
  TH2* hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom = new TH2F( "hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom", "delta eta delta phi plane of multiclusters relative to \"best\" one", 
                                                                 50, -0.2, 0.2, 50, -0.2, 0.2 );
  TH2* hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide = new TH2F( "hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide", "delta eta delta phi plane of multiclusters relative to \"best\" one", 
                                                                 50, -0.5, 0.5, 50, -0.2, 0.2 );
  TH2* hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted = new TH2F( "hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted", "delta eta delta phi plane of multiclusters relative to \"best\" one", 
                                                                 50, -0.2, 0.2, 50, -0.2, 0.2 );
  TH2* hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted = new TH2F( "hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted", "delta eta delta phi plane of multiclusters relative to \"best\" one", 
                                                                 50, -0.5, 0.5, 50, -0.2, 0.2 );
  TH2* hNumClusters2DPerMultiVsPt = new TH2F( "hNumClusters2DPerMultiVsPt", "no of 2D clusters in each multicluster vs pt", 20, 0., 40., numLayers+1, -0.5, numLayers+0.5 );
  TH2* hNumClusters2DPerBestMultiVsPt = new TH2F( "hNumClusters2DPerBestMultiVsPt", "no of 2D clusters in the best multicluster vs pt", 20, 0., 40., numLayers-9, 9.5, numLayers+0.5 );
  TH2* hNumClusters2DPerSecondBestMultiVsPt = new TH2F( "hNumClusters2DPerSecondBestMultiVsPt", "no of 2D clusters in the best multicluster vs pt", 20, 0., 40., numLayers, 0.5, numLayers+0.5 );
  TH2* hNumLayersPerMultiVsPt     = new TH2F( "hNumLayersPerMultiVsPt", "no of layers with a 2d cluster in each multicluster vs pt", numLayers+1, -0.5, numLayers+0.5, 50, 0, 50 );
  TH2* hSingleMCEFracVsNumSelected = new TH2F( "hSingleMCEFracVsNumSelected", "Number of selected multiclusters vs the fraction of gen energy in best multicluster", 100, 0., 2., 5, 0.5, 5.5 );

  TH1* hMultiDeltaRtoBest = new TH1F( "hMultiDeltaRtoBest", "deltaR of multiclusters relative to \"best\" one, with 1\% of best energy cut", 50, 0., 1. );

  TH1* hTotRecHitEFrac = new TH1F( "hTotRecHitEFrac",  "fraction of gen energy contained by all recHits", 100, 0., 2. );
  TH2* hRecHitEtaPhi = new TH2F( "hRecHitEtaPhi", " eta phi plane of recHits", 50, -3.14, 3.14, 50, -3.0, 3.0 );
  TH2* hRecHitEtaPhiZoom = new TH2F( "hRecHitEtaPhiZoom", " eta phi plane of recHits", 50, -0.2, 0.2, 50, -0.2, 0.2 );
  TH1* hRecHitPhiMean = new TH1F( "hRecHitPhiMean", "desc", 40, -0.01, 0.01 );
  TH1* hRecHitPhiRms = new TH1F( "hRecHitPhiRms", "desc", 40, 0., 0.04 );
  TH1* hRecHitEtaMean = new TH1F( "hRecHitEtaMean", "desc", 40, -0.01, 0.01 );
  TH1* hRecHitEtaRms = new TH1F( "hRecHitEtaRms", "desc", 40, 0., 0.04 );
  TH2* hRecHitLayerVsClusterZ = new TH2F( "hRecHitLayerVsClusterZ", "desc", 53, -0.5, 52.5, 50, 300., 500. );

  const int numRadii = 29;
  /*vector<TH1*> vecRecHitEnergyRadiusHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithin" + NumberToString( (i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( (i+2) );
    vecRecHitEnergyRadiusHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 40, 0.7, 1.1 ) );
  }
  vector<TH1*> vecRecHitEnergyBestHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithinBest" + NumberToString( (i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( (i+2) ) + " cm of best multicluster";
    vecRecHitEnergyBestHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 40, 0.7, 1.1 ) );
  }
  vector<TH1*> vecRecHitEnergySelectedHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithinSelected" + NumberToString( (i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( (i+2) ) + " cm from any selected multicluster";
    vecRecHitEnergySelectedHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 40, 0.7, 1.1 ) );
  }*/
  vector<TH1*> vecRecHitEnergyRadiusHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithin" + NumberToString( 0.5*(i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( 0.5*(i+2) );
    vecRecHitEnergyRadiusHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 100, 0., 2. ) );
  }
  vector<TH1*> vecRecHitEnergyBestHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithinBest" + NumberToString( 0.5*(i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( 0.5*(i+2) ) + " cm of best multicluster";
    //vecRecHitEnergyBestHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 100, 0., 2. ) );
    vecRecHitEnergyBestHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 40, 0.9, 1.3 ) );
  }
  vector<TH1*> vecRecHitEnergySelectedHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithinSelected" + NumberToString( 0.5*(i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( 0.5*(i+2) ) + " cm from any selected multicluster";
    //vecRecHitEnergySelectedHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 100, 0., 2. ) );
    vecRecHitEnergySelectedHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 40, 0.9, 1.3 ) );
  }
  vector<TH1*> vecRecHitEnergyContiguousHists;
  for( uint i=0; i<numRadii; i++ ) { 
    string histName = "hRecHitEFracWithinContiguous" + NumberToString( 0.5*(i+2) );
    string histDesc = "fraction of gen energy contained in recHits within radius of " + NumberToString( 0.5*(i+2) ) + " cm from any selected multicluster passing contiguity cut";
    vecRecHitEnergyContiguousHists.push_back( new TH1F( histName.c_str(),  histDesc.c_str(), 100, 0., 2. ) );
  }

  //const int numLayers = 28;
  vector<TH1*> vecNumSelectedMultiHists;
  for( uint i=0; i<numLayers; i++ ) {
    string histName = "hNumSelectedMultisHist" + NumberToString( i );
    string histDesc = "number of multiclusters with more than " + NumberToString( i ) + " 2D clusters";
    vecNumSelectedMultiHists.push_back( new TH1F( histName.c_str(), histDesc.c_str(), 11, -0.5, 10.5 ) );
  }
  vector<TH1*> vecSelectedMultiEnergyHists;
  for( uint i=0; i<numLayers; i++ ) {
    string histName = "hSelectedMultisEnergyHist" + NumberToString( i );
    string histDesc = "total energy in multiclusters with more than " + NumberToString( i ) + " 2D clusters";
    //vecSelectedMultiEnergyHists.push_back( new TH1F( histName.c_str(), histDesc.c_str(), 100, 0., 2. ) );
    vecSelectedMultiEnergyHists.push_back( new TH1F( histName.c_str(), histDesc.c_str(), 40, 0.9, 1.3 ) );
  }

  //const int numClusterCut = 3; // NB this is implemented as > than, NOT >= than.
  const int numClusterCut = 2; // NB this is implemented as > than, NOT >= than.

  // loop over events
  uint nEntries = theTree->GetEntries();
  if(debug) cout << "nEntries " << nEntries << endl;
  for( uint evtIndex = 0; evtIndex < nEntries; evtIndex++ )
  {
    if( evtIndex==0 ) cout << "Processing 0th event" << endl;
    else if( evtIndex%100 == 0 ) cout << "Processing " << evtIndex << "th event of " << nEntries << " total" << endl;
    theTree->GetEntry( evtIndex );
    if(debug) { if( evtIndex==0 ) cout << "Got entry of 0th event" << endl; }

    if( debug ) {
      cout << "GEN PHOTON ONE ETA,PHI = " << (*genParticles)[0].eta << ",    " << (*genParticles)[0].phi << endl;
      cout << "GEN PHOTON TWO ETA,PHI = " << (*genParticles)[1].eta << ",    " << (*genParticles)[1].phi << endl;
    }

    // loop over gen particles
    uint nGenParticles = genParticles->size();
    if(debug) cout << "nGenParticles " << nGenParticles << endl;
    for( uint gpIndex = 0; gpIndex < nGenParticles; gpIndex++ )
    {
      if(debug) cout << "gen particle " << gpIndex << " has energy " << (*genParticles)[gpIndex].energy << endl;
      if(debug) cout << "gen particle " << gpIndex << " has eta "    << (*genParticles)[gpIndex].eta    << endl;

      float genEta    = (*genParticles)[gpIndex].eta;
      float genPhi    = (*genParticles)[gpIndex].phi;
      float genEnergy = (*genParticles)[gpIndex].energy;
      float genPt     = (*genParticles)[gpIndex].pt;
      float genDz     = (*genParticles)[gpIndex].dvz;
      if(debug) cout << "genDz = " << genDz << endl;

      // require good containment
      if( abs(genEta) < etaMin || abs(genEta) > etaMax ) continue; 
      // fill prerequisites for conversion ratio histogram
      if( abs(genDz) < dvzCut ) hGenConv->Fill( abs(genEta) );
      else if ( abs(genDz) > dvzCut ) hGenUnconv->Fill( abs(genEta) );
      if( testConversion ) {
	// skip (un)converted events as necessary
	if( converted    && abs(genDz) > dvzCut ) continue;
	if( (!converted) && abs(genDz) < dvzCut ) continue;
      }

      hGenEta->Fill( genEta );
      hGenPhi->Fill( genPhi );
      hGenEnergy->Fill( genEnergy );
      hGenPt->Fill( genPt );

      //loop over pfClusters, count them and find best 
      int numPfClusters = 0;
      float bestPfEnergy = 0.;
      int bestPfIndex= -1.;
      for( uint pfIndex=0; pfIndex<pfClusters->size(); pfIndex++ ) {
        if( pfClusters->at(pfIndex).eta*genEta < 0 ) continue;
        numPfClusters++;
        if( pfClusters->at(pfIndex).energy > bestPfEnergy ) {
          bestPfEnergy = pfClusters->at(pfIndex).energy;
          bestPfIndex= pfIndex;
        }
      }
      if(debug) cout << "numPfClusters = " << numPfClusters << endl;
      if( convType == "SinglePfCluster" && numPfClusters != 1 ) continue;
      hBestPfEFrac->Fill( pfClusters->at(bestPfIndex).energy / genEnergy );

      // loop over multiclusters and match to gen particles using dR
      float matchedMultiClusterEnergy = 0.;
      float matchedMultiClusterPt     = 0.;
      float matchedMultiClusterEta    = 0.;
      int   matchedMultiClusterIndex  = -1; // used later to only include 2d clusters from "correct" multicluster
      int   numHighPtMultis = 0;
      vector<int> vecNumSelectedMultis(numLayers,0);
      vector< vector<int> > vecVecSelectedMultiIndices(numLayers);
      if(debug) cout << "Entering first loop over multiclusters" << endl;
      for( uint mcIndex = 0; mcIndex < multiClusters->size(); mcIndex++ )
      {
        float mcEta = multiClusters->at( mcIndex ).eta;
        if( mcEta*genEta < 0. ) continue;
        //if( multiClusters->at(mcIndex).nclus != 1 ) continue; //characterising the single-layer multiclusters
        float mcPhi = multiClusters->at( mcIndex ).phi;
        float tempDr = deltaR( mcEta, genEta, mcPhi, genPhi );
        if( multiClusters->at( mcIndex ).energy > matchedMultiClusterEnergy ) {
          matchedMultiClusterEnergy = multiClusters->at( mcIndex ).energy;
          matchedMultiClusterPt     = multiClusters->at( mcIndex ).pt;
          matchedMultiClusterEta     = multiClusters->at( mcIndex ).eta;
          matchedMultiClusterIndex  = mcIndex;
        }
        if( multiClusters->at( mcIndex ).pt > 0.8*35 ) numHighPtMultis++;
        for( uint selIndex=0; selIndex<vecNumSelectedMultis.size(); selIndex++ )
        {
          if( multiClusters->at( mcIndex ).nclus > selIndex ) {
            vecNumSelectedMultis[selIndex] += 1;
            vecVecSelectedMultiIndices[selIndex].push_back( mcIndex );
          }
        }
        hNumClusters2DPerMultiVsPt->Fill( multiClusters->at(mcIndex).pt, multiClusters->at(mcIndex).nclus );
      } // end of first loop over multiclusters
      if( genEnergy > 222 && genEnergy < 223 && debug ) {
	cout << "event number = " << evtIndex << endl;
	cout << "gen Energy = " << genEnergy << endl;
	cout << "number of 2Ds in best multicluster = " << multiClusters->at( matchedMultiClusterIndex ).nclus << endl;
	cout << "number of multiclusters = " << multiClusters->size() << endl << endl;
      }
      if(debug) cout << "Exiting first loop over multiclusters" << endl;
      hSingleMCEFrac->Fill(  matchedMultiClusterEnergy / genEnergy );
      hSingleMCPtFrac->Fill( matchedMultiClusterPt     / genPt     );
      hSingleMCPtFracVsEta->Fill( matchedMultiClusterEta, matchedMultiClusterPt / genPt     );
      hSingleMCGenDr->Fill( deltaR( matchedMultiClusterEta, genEta, multiClusters->at( matchedMultiClusterIndex ).phi, genPhi ) );
      hSingleMCGenDeta->Fill(  matchedMultiClusterEta - genEta );
      hSingleMCGenDphi->Fill(  multiClusters->at( matchedMultiClusterIndex ).phi - genPhi );
      hSingleMCGenDetaDphi->Fill(  multiClusters->at( matchedMultiClusterIndex ).phi - genPhi, matchedMultiClusterEta - genEta );
      hSingleMCGenDetaDphiWeighted->Fill(  multiClusters->at( matchedMultiClusterIndex ).phi - genPhi, matchedMultiClusterEta - genEta, matchedMultiClusterEnergy );
      hNumMultiClustersHighPt->Fill( numHighPtMultis );
      vector<float> vecTotalSelectedEnergy(numLayers,0.);
      vector< vector<int> > vecVecSelectedMulti2DLayers( vecNumSelectedMultis[numClusterCut] ); // for selecting multis with >=numCluscterCut CONTIGUOUS 2D clusters
      for( uint i=0; i<vecNumSelectedMultis.size(); i++ ) {
        vecNumSelectedMultiHists[i]->Fill( vecNumSelectedMultis[i] );
        if( i==numClusterCut ) hSingleMCEFracVsNumSelected->Fill( multiClusters->at( matchedMultiClusterIndex ).energy / genEnergy, vecNumSelectedMultis[i] );
        if( i==numClusterCut ) hNumVsEta->Fill( multiClusters->at( matchedMultiClusterIndex ).eta );
        if( i==numClusterCut && vecNumSelectedMultis[i]>1 ) hNumWithExtraVsEta->Fill( multiClusters->at( matchedMultiClusterIndex ).eta );
        for( uint j=0; j<vecVecSelectedMultiIndices[i].size(); j++ ) {
          vecTotalSelectedEnergy[i] += multiClusters->at( vecVecSelectedMultiIndices[i][j] ).energy;
          if( (i == numClusterCut) && (vecVecSelectedMultiIndices[i][j] != matchedMultiClusterIndex) ) {
	    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->Fill( deltaPhi( multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).phi ), 
								 deltaEta( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).eta ) );
	    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->Fill( deltaPhi( multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).phi ), 
								 deltaEta( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).eta ) );
	    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->Fill( deltaPhi( multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).phi ), 
								 deltaEta( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).eta ), 
                                                                 multiClusters->at( vecVecSelectedMultiIndices[i][j] ).energy );
	    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->Fill( deltaPhi( multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).phi ), 
								 deltaEta( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( vecVecSelectedMultiIndices[i][j] ).eta ),
                                                                 multiClusters->at( vecVecSelectedMultiIndices[i][j] ).energy );
            /*cout << "evtIndex = " << evtIndex << endl;
            cout << "gpIndex  = " << gpIndex << endl;
	    cout << "pt    = " << multiClusters->at( vecVecSelectedMultiIndices[i][j] ).pt << endl;
	    cout << "eta   = " << multiClusters->at( vecVecSelectedMultiIndices[i][j] ).eta << endl;
	    cout << "nclus = " << multiClusters->at( vecVecSelectedMultiIndices[i][j] ).nclus << endl;
	    cout << "seed value = " << multiClusters->at( vecVecSelectedMultiIndices[i][j] ).cl2dSeed << endl;
	    cout << "seed layer = " << clusters2D->at( multiClusters->at( vecVecSelectedMultiIndices[i][j] ).cl2dSeed ).layer << endl;*/
            float prevX = 0.;
            float prevY = 0.;
            for( uint tempIndex=0; tempIndex<clusters2D->size(); tempIndex++ ) {
              if( clusters2D->at( tempIndex ).multicluster != vecVecSelectedMultiIndices[i][j] ) continue;
              vecVecSelectedMulti2DLayers[j].push_back( clusters2D->at( tempIndex ).layer );
              float tempX = clusters2D->at( tempIndex ).x;
              float tempY = clusters2D->at( tempIndex ).y;
              if( abs(prevX) > 0.00001 ) hDistanceBetweenSelected->Fill( deltaX( prevX, tempX, prevY, tempY ) );
              prevX = tempX;
              prevY = tempY;
              /*cout << "consitutent 2D cluster number " << tempIndex << endl;
              cout << "layer = " << clusters2D->at( tempIndex ).layer << endl;
              cout << "eta   = " << clusters2D->at( tempIndex ).eta << endl;
              cout << "phi   = " << clusters2D->at( tempIndex ).phi << endl;
              cout << "pt    = " << clusters2D->at( tempIndex ).pt << endl;*/
              hLayerOfClustersInSelected->Fill( clusters2D->at( tempIndex ).layer );
              hLayerOfClustersInSelectedWeighted->Fill( clusters2D->at( tempIndex ).layer, clusters2D->at( tempIndex ).energy );
            } // end of loop over 2D clusters
            //cout << endl;
          }
          else if( (i == numClusterCut) && (vecVecSelectedMultiIndices[i][j] == matchedMultiClusterIndex) ) {
            for( uint tempIndex=0; tempIndex<clusters2D->size(); tempIndex++ ) {
              if( clusters2D->at( tempIndex ).multicluster != matchedMultiClusterIndex ) continue;
              vecVecSelectedMulti2DLayers[j].push_back( clusters2D->at( tempIndex ).layer );
              /*cout << "consitutent 2D cluster number " << tempIndex << endl;
              cout << "layer = " << clusters2D->at( tempIndex ).layer << endl;
              cout << "eta   = " << clusters2D->at( tempIndex ).eta << endl;
              cout << "phi   = " << clusters2D->at( tempIndex ).phi << endl;
              cout << "pt    = " << clusters2D->at( tempIndex ).pt << endl;*/
              hLayerOfClustersInBest->Fill( clusters2D->at( tempIndex ).layer );
              hLayerOfClustersInBestWeighted->Fill( clusters2D->at( tempIndex ).layer, clusters2D->at( tempIndex ).energy );
            }
          }
        }
        vecSelectedMultiEnergyHists[i]->Fill( vecTotalSelectedEnergy[i] / genEnergy );
      }
      hNumClusters2DPerBestMultiVsPt->Fill( multiClusters->at(matchedMultiClusterIndex).pt, multiClusters->at(matchedMultiClusterIndex).nclus );

      //testing deltaphi correction 
      /*cout << "genEta       = " << genEta << endl;
      float genPz = genEnergy * cos( 2*atan( exp( -1*genEta ) ) );
      cout << "genPz        = " << genPz << endl;
      float phiCorrection = -320.*0.001*2.998*3.8 / genPz;
      cout << "genPhi       = " << genPhi << endl;
      cout << "dPhi         = " << phiCorrection << endl;
      cout << "correctedPhi = " << genPhi + phiCorrection << endl;
      cout << "mcPhi        = " << multiClusters->at(matchedMultiClusterIndex).phi << endl;
      cout << "old diff = " << multiClusters->at(matchedMultiClusterIndex).phi - genPhi << endl;
      cout << "new diff = " << multiClusters->at(matchedMultiClusterIndex).phi - (genPhi + phiCorrection) << endl << endl;*/

      // do selection on contiguous clusters
      int numSelectedContiguous = 0;
      vector< int > vecMultisPassingContigIndices;
      for( uint contigIndex=0; contigIndex<vecVecSelectedMulti2DLayers.size(); contigIndex++ ) {
        if(debug) {
	  cout << "j = " << contigIndex << endl;
	  cout << "vecVecSelectedMulti2DLayers[j].size() = " << vecVecSelectedMulti2DLayers[contigIndex].size() << endl;
	  cout << "vecVecSelectedMultiIndices[numClusterCut][j].nclus = " << multiClusters->at( vecVecSelectedMultiIndices[numClusterCut][contigIndex] ).nclus << endl;
        }
        vector<int> vecLayers = vecVecSelectedMulti2DLayers[contigIndex];
        sort( vecLayers.begin(), vecLayers.end() );
        int prevLayer = -999;
        int contigCount = 1; // got this wrong first time - should obviously start from one
        for( uint layerIndex=0; layerIndex<vecLayers.size(); layerIndex++ ) {
          if( vecLayers[layerIndex] == prevLayer + 1 ) contigCount++;
          else contigCount = 1;
          if( contigCount == 3 ) break;
          prevLayer = vecLayers[layerIndex];
        }
        if( contigCount == 3 ) { 
          numSelectedContiguous++;
          vecVecSelectedMulti2DLayers[contigIndex][0] = 1000;
          vecMultisPassingContigIndices.push_back( vecVecSelectedMultiIndices[numClusterCut][contigIndex] );
        }
      }
      if(debug) {
        cout << "numSelectedMultis     = " << vecNumSelectedMultis[numClusterCut] << endl;
        cout << "numSelectedContiguous = " << numSelectedContiguous << endl;
      }
      hNumSelectedMultis->Fill( vecNumSelectedMultis[numClusterCut] );
      hNumSelectedMultisContiguous->Fill( numSelectedContiguous );
      if( numSelectedContiguous > 1 ) hNumWithExtraContigVsEta->Fill( multiClusters->at( matchedMultiClusterIndex ).eta );
      for( uint i=0; i<vecMultisPassingContigIndices.size(); i++ ) {
        if( vecMultisPassingContigIndices[i] != matchedMultiClusterIndex ) { 
          int contigIndex = vecMultisPassingContigIndices[i];
          /*float matchedX = etaPhiZtoX( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( matchedMultiClusterIndex ).z );
          float matchedY = etaPhiZtoY( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( matchedMultiClusterIndex ).z );
          float contigX = etaPhiZtoX( multiClusters->at( contigIndex ).eta, multiClusters->at( contigIndex ).phi, multiClusters->at( contigIndex ).z );
          float contigY = etaPhiZtoY( multiClusters->at( contigIndex ).eta, multiClusters->at( contigIndex ).phi, multiClusters->at( contigIndex ).z );*/
          float matchedX = etaPhiZtoX( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( matchedMultiClusterIndex ).phi, 320. );
          float matchedY = etaPhiZtoY( multiClusters->at( matchedMultiClusterIndex ).eta, multiClusters->at( matchedMultiClusterIndex ).phi, 320. );
          float contigX = etaPhiZtoX( multiClusters->at( contigIndex ).eta, multiClusters->at( contigIndex ).phi, 320. );
          float contigY = etaPhiZtoY( multiClusters->at( contigIndex ).eta, multiClusters->at( contigIndex ).phi, 320. ); // changed to constant val to allow better comparison
	  hDistanceContigBest->Fill( deltaX( matchedX, contigX, matchedY, contigY ) );
	  hMultisPassingContigPtFrac->Fill( multiClusters->at( vecMultisPassingContigIndices[i] ).pt / genPt );
        }
      }

      // energy in contiguous
      float contigEnergy = 0;
      for( uint contigIndex=0; contigIndex<vecNumSelectedMultis[numClusterCut]; contigIndex++ ) {
        if( vecVecSelectedMulti2DLayers[contigIndex][0] == 1000 ) contigEnergy += multiClusters->at( vecVecSelectedMultiIndices[numClusterCut][contigIndex] ).energy;
      }
      if(debug) cout << "contigEnergy / genEnergy = " << contigEnergy / genEnergy << endl;
      hEFracContiguous->Fill( contigEnergy / genEnergy );

      //energy in best MC EE only
      float singleMCLimitedEnergy = 0.;
      float totalClusterEnergy = 0.;
      for( uint limitIndex=0; limitIndex<clusters2D->size(); limitIndex++ ) {
        if( clusters2D->at( limitIndex ).eta*genEta<0. ) continue;
        if( clusters2D->at( limitIndex ).multicluster == matchedMultiClusterIndex && clusters2D->at( limitIndex ).layer<=28 ) singleMCLimitedEnergy += clusters2D->at( limitIndex ).energy;
        totalClusterEnergy += clusters2D->at( limitIndex ).energy;
        if( clusters2D->at( limitIndex ).layer > 40 ) hBHclusterEta->Fill( clusters2D->at( limitIndex ).eta );
        if( clusters2D->at( limitIndex ).layer > 40 ) hBHclusterPhi->Fill( clusters2D->at( limitIndex ).phi );
      }
      hSingleMCLimitedEFrac->Fill( singleMCLimitedEnergy / genEnergy );
      hTotalClusterEFrac->Fill( totalClusterEnergy / genEnergy );

      // now loop again for deltaEta deltaPhi plot and second highest multicluster
      float secondMultiClusterEnergy = 0.;
      float secondMultiClusterPt     = 0.;
      float secondMultiClusterEta    = 0.;
      float secondMinDr              = 1000.;
      int   secondMultiClusterIndex  = -1; // used later to only include 2d clusters from "correct" multicluster
      float pseudoSuperClusterEnergy = matchedMultiClusterEnergy;
      if(debug) cout << "Entering second loop over multiclusters" << endl;
      for( uint mcIndex(0); mcIndex<multiClusters->size(); mcIndex++ ) {
        if( multiClusters->at( mcIndex ).eta*genEta < 0 ) continue;
        if( mcIndex == matchedMultiClusterIndex ) continue;
        if( multiClusters->at(mcIndex).nclus != 1 ) continue; //characterising the single-layer multiclusters
        float bestEta = multiClusters->at( matchedMultiClusterIndex ).eta;
        float bestPhi = multiClusters->at( matchedMultiClusterIndex ).phi;
        float tempEta = multiClusters->at( mcIndex ).eta;
        float tempPhi = multiClusters->at( mcIndex ).phi;
        hMultiClusterDeltaEtaDeltaPhi->Fill(          deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta) );
        hMultiClusterDeltaEtaDeltaPhiWeighted->Fill(  deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta), multiClusters->at( mcIndex ).energy );
        hMultiClusterDeltaEtaDeltaPhiZoomed->Fill(          deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta) );
        hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->Fill(  deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta), multiClusters->at( mcIndex ).energy );
        hMultiClusterDeltaEtaDeltaPhiZoomZoom->Fill(          deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta) );
        hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->Fill(  deltaPhi(bestPhi, tempPhi), deltaEta(bestEta, tempEta), multiClusters->at( mcIndex ).energy );
        
        float tempEnergy = multiClusters->at( mcIndex ).energy;
        if( tempEnergy > secondMultiClusterEnergy ) {
          secondMultiClusterEnergy = tempEnergy;
          secondMultiClusterPt     = multiClusters->at( mcIndex ).pt;
          secondMultiClusterEta    = tempEta;
          secondMultiClusterIndex  = mcIndex;
        }

        float tempDr = deltaR( tempEta, bestEta, tempEta, bestPhi );
        if( tempEnergy > 0.01*matchedMultiClusterEnergy ) hMultiDeltaRtoBest->Fill( tempDr );
        // want to add a loop over possible depth superclustering dR values here
        if( tempDr < 0.1 && tempEnergy > 0.01*matchedMultiClusterEnergy ) pseudoSuperClusterEnergy += tempEnergy;
      } // end of second loop over multiclusters
      if(debug) cout << "Exiting second loop over multiclusters" << endl;
      hDoubleMCEFrac->Fill(  (matchedMultiClusterEnergy + secondMultiClusterEnergy) / genEnergy );
      hDoubleMCPtFrac->Fill( (matchedMultiClusterPt + secondMultiClusterPt)   / genPt     );
      hSecondMCEFrac->Fill(  (secondMultiClusterEnergy) / genEnergy );
      hSecondMCPtFrac->Fill( (secondMultiClusterPt)   / genPt     );
      hSecondRelEta->Fill( multiClusters->at( secondMultiClusterIndex ).eta - multiClusters->at( matchedMultiClusterIndex ).eta );
      hSecondRelPhi->Fill( multiClusters->at( secondMultiClusterIndex ).phi - multiClusters->at( matchedMultiClusterIndex ).phi);
      hSecondRelDr->Fill( deltaR( multiClusters->at( secondMultiClusterIndex ).eta, multiClusters->at( matchedMultiClusterIndex ).eta, 
                                   multiClusters->at( matchedMultiClusterIndex ).phi, multiClusters->at( matchedMultiClusterIndex ).phi ) );
      hNumClusters2DPerSecondBestMultiVsPt->Fill( multiClusters->at( secondMultiClusterIndex ).pt, multiClusters->at( secondMultiClusterIndex ).nclus );
      hPseudoSuperClusterEFrac->Fill( pseudoSuperClusterEnergy / genEnergy );

      if( debug ) {
        cout << "best MC eta = " << multiClusters->at( matchedMultiClusterIndex ).eta << endl;
        cout << "gen eta     = " << genEta << endl;
        cout << "best MC phi = " << multiClusters->at( matchedMultiClusterIndex ).phi << endl;
        cout << "gen phi     = " << genPhi << endl << endl;
      }

      // same for recHits. Have three energy containment plots: within r cm of gen, within r cm of best multi, and within r com of any selected multi
      vector< pair<float,float> > vecRecHitEnergyWithinRadius;
      vecRecHitEnergyWithinRadius.resize( numRadii );
      //for( uint i=0; i<vecRecHitEnergyWithinRadius.size(); i++ ) vecRecHitEnergyWithinRadius[i] = make_pair( 1.*(i+2), 0. );
      for( uint i=0; i<vecRecHitEnergyWithinRadius.size(); i++ ) vecRecHitEnergyWithinRadius[i] = make_pair( 0.5*(i+2), 0. );
      vector< pair<float,float> > vecRecHitEnergyWithinBest;
      vecRecHitEnergyWithinBest.resize( numRadii );
      //for( uint i=0; i<vecRecHitEnergyWithinBest.size(); i++ ) vecRecHitEnergyWithinBest[i] = make_pair( 1.*(i+2), 0. );
      for( uint i=0; i<vecRecHitEnergyWithinBest.size(); i++ ) vecRecHitEnergyWithinBest[i] = make_pair( 0.5*(i+2), 0. );
      vector< pair<float,float> > vecRecHitEnergyWithinSelected;
      vecRecHitEnergyWithinSelected.resize( numRadii );
      //for( uint i=0; i<vecRecHitEnergyWithinSelected.size(); i++ ) vecRecHitEnergyWithinSelected[i] = make_pair( 1.*(i+2), 0. );
      for( uint i=0; i<vecRecHitEnergyWithinSelected.size(); i++ ) vecRecHitEnergyWithinSelected[i] = make_pair( 0.5*(i+2), 0. );
      vector< pair<float,float> > vecRecHitEnergyWithinContiguous;
      vecRecHitEnergyWithinContiguous.resize( numRadii );
      //for( uint i=0; i<vecRecHitEnergyWithinContiguous.size(); i++ ) vecRecHitEnergyWithinContiguous[i] = make_pair( 1.*(i+2), 0. );
      for( uint i=0; i<vecRecHitEnergyWithinContiguous.size(); i++ ) vecRecHitEnergyWithinContiguous[i] = make_pair( 0.5*(i+2), 0. );
      if(debug) cout << "Entering loop over recHits" << endl;
      float rhTotEnergy = 0.;
      TH2* hTemp = new TH2F( "hTemp", "e weighted eta phi plane of recHits", 50, -0.2, 0.2, 50, -0.2, 0.2 );
      for( uint rhIndex = 0; rhIndex < recHits->size(); rhIndex++ )
      {

        float recHitEta      = recHits->at( rhIndex ).eta;
        if( recHitEta*genEta < 0 ) continue;
        float recHitPhi      = recHits->at( rhIndex ).phi;
        float recHitX        = recHits->at( rhIndex ).x;
        float recHitY        = recHits->at( rhIndex ).y;
        float recHitEnergy   = recHits->at( rhIndex ).energy;
        float recHitLayer    = recHits->at( rhIndex ).layer;
        float recHitCluster  = recHits->at( rhIndex ).cluster2d;
        hRecHitLayerVsClusterZ->Fill( recHitLayer, clusters2D->at(recHitCluster).z );
        if(recHitLayer > 40 ) hBHrecHitEta->Fill( recHitEta );
        if(recHitLayer > 40 ) hBHrecHitPhi->Fill( recHitPhi );
        if( recHitLayer > layerCut ) continue;
        if( recHitEnergy > 0.001*genEnergy ) hRecHitEtaPhi->Fill( deltaPhi(genPhi,recHitPhi) , deltaEta(genEta,recHitEta) );
        if( recHitEnergy > 0.001*genEnergy ) hRecHitEtaPhiZoom->Fill( deltaPhi(genPhi,recHitPhi) , deltaEta(genEta,recHitEta) );
        hTemp->Fill( deltaPhi(genPhi,recHitPhi) , deltaEta(genEta,recHitEta), recHitEnergy );
        rhTotEnergy += recHitEnergy; 
        if(debug) {
          cout << "rh actual z    = " << recHits->at( rhIndex ).z << endl;
          cout << "rh layer to z  = " << layerToZ( recHitLayer, recHitEta ) << endl;
          cout << "rh actual x    = " << recHits->at( rhIndex ).x << endl;
          cout << "rh etaphi to x = " << etaPhiZtoX( recHitEta, recHitPhi, recHits->at( rhIndex ).z ) << endl;
        }
        float genRecHitDs    = dsGenRecHit( genEta, genPhi, recHitLayer, recHitX, recHitY ); // now giving good results on old ntuple
        float mcBestEta      = multiClusters->at( matchedMultiClusterIndex ).eta;
        float mcBestPhi      = multiClusters->at( matchedMultiClusterIndex ).phi;
        float mcBestRecHitDs = dsGenRecHit( mcBestEta, mcBestPhi, recHitLayer, recHitX, recHitY );
        for( uint i=0; i<vecRecHitEnergyWithinRadius.size(); i++ ) {
          if( genRecHitDs    < vecRecHitEnergyWithinRadius[i].first ) vecRecHitEnergyWithinRadius[i].second += recHitEnergy;
          if( mcBestRecHitDs < vecRecHitEnergyWithinBest[i].first   ) vecRecHitEnergyWithinBest[i].second   += recHitEnergy;

	  if( doSelectedMultis ) {
	    for( uint j=0; j<vecVecSelectedMultiIndices[numClusterCut].size(); j++ ) {
	      float mcTempEta = multiClusters->at( vecVecSelectedMultiIndices[numClusterCut][j] ).eta;
	      float mcTempPhi = multiClusters->at( vecVecSelectedMultiIndices[numClusterCut][j] ).phi;
	      float mcTempRecHitDs = dsGenRecHit( mcTempEta, mcTempPhi, recHitLayer, recHitX, recHitY );
	      if( mcTempRecHitDs  < vecRecHitEnergyWithinSelected[i].first ) {
		vecRecHitEnergyWithinSelected[i].second += recHitEnergy;
		for( uint k=0; k<vecMultisPassingContigIndices.size(); k++ ) {
		  if( vecVecSelectedMultiIndices[numClusterCut][j] == vecMultisPassingContigIndices[k] ) vecRecHitEnergyWithinContiguous[i].second += recHitEnergy;
		}
		break;
	      }
	    }
          }
        }
      }  // end of loop over recHits
      double tempStats[7];
      hTemp->GetStats(tempStats);
      delete hTemp;
      hRecHitPhiMean->Fill( tempStats[2]/tempStats[0]  );
      hRecHitPhiRms->Fill( sqrt( tempStats[3]/tempStats[0] ) );
      hRecHitEtaMean->Fill( tempStats[4]/tempStats[0] );
      hRecHitEtaRms->Fill( sqrt( tempStats[5]/tempStats[0] ) );
      hTotRecHitEFrac->Fill( rhTotEnergy / genEnergy);
      if(debug) cout << "Exiting loop over recHits" << endl;
      for( uint i=0; i<vecRecHitEnergyRadiusHists.size(); i++ ) vecRecHitEnergyRadiusHists[i]->Fill( vecRecHitEnergyWithinRadius[i].second / genEnergy );
      for( uint i=0; i<vecRecHitEnergyBestHists.size(); i++ ) vecRecHitEnergyBestHists[i]->Fill( vecRecHitEnergyWithinBest[i].second / genEnergy );
      for( uint i=0; i<vecRecHitEnergySelectedHists.size(); i++ ) vecRecHitEnergySelectedHists[i]->Fill( vecRecHitEnergyWithinSelected[i].second / genEnergy );
      for( uint i=0; i<vecRecHitEnergyContiguousHists.size(); i++ ) vecRecHitEnergyContiguousHists[i]->Fill( vecRecHitEnergyWithinContiguous[i].second / genEnergy );
  
    } // end of loop over gen particles

  } // end of loop over events
  cout << "Finished processing events" << endl;

  hGenUnconv->Add( hGenConv ); // now total photons so can be used for ratio
  TGraphAsymmErrors* gConvRatio = new TGraphAsymmErrors( hGenConv, hGenUnconv );
  gConvRatio->SetName("gConvRatio");

  TGraphAsymmErrors* gExtraRatio = new TGraphAsymmErrors( hNumWithExtraVsEta, hNumVsEta );
  gExtraRatio->SetName("gExtraRatio");

  TGraphAsymmErrors* gContigRatio = new TGraphAsymmErrors( hNumWithExtraContigVsEta, hNumVsEta );
  gContigRatio->SetName("gContigRatio");
  
  TGraph* gNumSelectedVsCut = new TGraph( 20 );
  gNumSelectedVsCut->SetName("gNumSelectedVsCut"); 
  TGraph* gSelectedEnergyVsCut = new TGraph( 20 );
  gSelectedEnergyVsCut->SetName("gSelectedEnergyVsCut");

  TGraph* gRecHitEFracMeanVsRadius  = new TGraph( numRadii );
  gRecHitEFracMeanVsRadius->SetName("gRecHitEFracMeanVsRadius");
  TGraph* gRecHitEFracSigmaVsRadius = new TGraph( numRadii );
  gRecHitEFracSigmaVsRadius->SetName("gRecHitEFracSigmaVsRadius");
  TGraph* gRecHitEFracMeanBestVsRadius  = new TGraph( numRadii );
  gRecHitEFracMeanBestVsRadius->SetName("gRecHitEFracMeanBestVsRadius");
  TGraph* gRecHitEFracSigmaBestVsRadius = new TGraph( numRadii );
  gRecHitEFracSigmaBestVsRadius->SetName("gRecHitEFracSigmaBestVsRadius");
  TGraph* gRecHitEFracResolutionBestVsRadius = new TGraph( numRadii );
  gRecHitEFracResolutionBestVsRadius->SetName("gRecHitEFracResolutionBestVsRadius");
  TGraph* gRecHitEFracMeanSelectedVsRadius  = new TGraph( numRadii );
  gRecHitEFracMeanSelectedVsRadius->SetName("gRecHitEFracMeanSelectedVsRadius");
  TGraph* gRecHitEFracSigmaSelectedVsRadius = new TGraph( numRadii );
  gRecHitEFracSigmaSelectedVsRadius->SetName("gRecHitEFracSigmaSelectedVsRadius");
  TGraph* gRecHitEFracResolutionSelectedVsRadius = new TGraph( numRadii );
  gRecHitEFracResolutionSelectedVsRadius->SetName("gRecHitEFracResolutionSelectedVsRadius");
  TGraph* gRecHitEFracMeanContiguousVsRadius  = new TGraph( numRadii );
  gRecHitEFracMeanContiguousVsRadius->SetName("gRecHitEFracMeanContiguousVsRadius");
  TGraph* gRecHitEFracSigmaContiguousVsRadius = new TGraph( numRadii );
  gRecHitEFracSigmaContiguousVsRadius->SetName("gRecHitEFracSigmaContiguousVsRadius");
  TGraph* gRecHitEFracResolutionContiguousVsRadius = new TGraph( numRadii );
  gRecHitEFracResolutionContiguousVsRadius->SetName("gRecHitEFracResolutionContiguousVsRadius");

  // draw and save histograms
  if( drawAndSave ) {
    TCanvas *c0 = new TCanvas();
    hGenEta->Draw("hist");
    hGenEta->GetXaxis()->SetTitleSize(0.05);
    hGenEta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenEta.pdf" );
    c0->Print( "HGCPlots/hGenEta.png" );
    hGenPhi->Draw("hist");
    hGenPhi->GetXaxis()->SetTitleSize(0.05);
    hGenPhi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenPhi.pdf" );
    c0->Print( "HGCPlots/hGenPhi.png" );
    hGenEnergy->Draw("hist");
    hGenEnergy->GetXaxis()->SetTitleSize(0.05);
    hGenEnergy->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenEnergy.pdf" );
    c0->Print( "HGCPlots/hGenEnergy.png" );
    hGenPt->Draw("hist");
    hGenPt->GetXaxis()->SetTitleSize(0.05);
    hGenPt->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenPt.pdf" );
    c0->Print( "HGCPlots/hGenPt.png" );
    hGenConv->Draw("hist");
    hGenConv->GetXaxis()->SetTitleSize(0.05);
    hGenConv->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenConv.pdf" );
    c0->Print( "HGCPlots/hGenConv.png" );
    hGenUnconv->Add( hGenConv, -1. );
    hGenUnconv->Draw("hist");
    hGenUnconv->GetXaxis()->SetTitleSize(0.05);
    hGenUnconv->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hGenUnconv.pdf" );
    c0->Print( "HGCPlots/hGenUnconv.png" );
    hSingleMCEFrac->Draw("hist");
    hSingleMCEFrac->GetXaxis()->SetTitleSize(0.05);
    hSingleMCEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCEFrac.pdf" );
    c0->Print( "HGCPlots/hSingleMCEFrac.png" );
    hSingleMCLimitedEFrac->Draw("hist");
    hSingleMCLimitedEFrac->GetXaxis()->SetTitleSize(0.05);
    hSingleMCLimitedEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCLimitedEFrac.pdf" );
    c0->Print( "HGCPlots/hSingleMCLimitedEFrac.png" );
    hTotalClusterEFrac->Draw("hist");
    hTotalClusterEFrac->GetXaxis()->SetTitleSize(0.05);
    hTotalClusterEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hTotalClusterEFrac.pdf" );
    c0->Print( "HGCPlots/hTotalClusterEFrac.png" );
    hPseudoSuperClusterEFrac->Draw("hist");
    hPseudoSuperClusterEFrac->GetXaxis()->SetTitleSize(0.05);
    hPseudoSuperClusterEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hPseudoSuperClusterEFrac.pdf" );
    c0->Print( "HGCPlots/hPseudoSuperClusterEFrac.png" );
    hBHclusterEta->Draw("hist");
    hBHclusterEta->GetXaxis()->SetTitleSize(0.05);
    hBHclusterEta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hBHclusterEta.pdf" );
    c0->Print( "HGCPlots/hBHclusterEta.png" );
    hBHclusterPhi->Draw("hist");
    hBHclusterPhi->GetXaxis()->SetTitleSize(0.05);
    hBHclusterPhi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hBHclusterPhi.pdf" );
    c0->Print( "HGCPlots/hBHclusterPhi.png" );
    hBHrecHitEta->Draw("hist");
    hBHrecHitEta->GetXaxis()->SetTitleSize(0.05);
    hBHrecHitEta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hBHrecHitEta.pdf" );
    c0->Print( "HGCPlots/hBHrecHitEta.png" );
    hBHrecHitPhi->Draw("hist");
    hBHrecHitPhi->GetXaxis()->SetTitleSize(0.05);
    hBHrecHitPhi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hBHrecHitPhi.pdf" );
    c0->Print( "HGCPlots/hBHrecHitPhi.png" );
    hBestPfEFrac->Draw("hist");
    hBestPfEFrac->GetXaxis()->SetTitleSize(0.05);
    hBestPfEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hBestPfEFrac.pdf" );
    c0->Print( "HGCPlots/hBestPfEFrac.png" );
    hSingleMCPtFrac->Draw("hist");
    hSingleMCPtFrac->GetXaxis()->SetTitle("p_{T} multicluster / p_{T} gen");
    hSingleMCPtFrac->GetXaxis()->SetTitleSize(0.05);
    hSingleMCPtFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCPtFrac.pdf" );
    c0->Print( "HGCPlots/hSingleMCPtFrac.png" );
    hSingleMCGenDr->Draw("hist");
    hSingleMCGenDr->GetXaxis()->SetTitleSize(0.05);
    hSingleMCGenDr->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCGenDr.pdf" );
    c0->Print( "HGCPlots/hSingleMCGenDr.png" );
    hSingleMCGenDeta->Draw("hist");
    hSingleMCGenDeta->GetXaxis()->SetTitleSize(0.05);
    hSingleMCGenDeta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCGenDeta.pdf" );
    c0->Print( "HGCPlots/hSingleMCGenDeta.png" );
    hSingleMCGenDphi->Draw("hist");
    hSingleMCGenDphi->GetXaxis()->SetTitleSize(0.05);
    hSingleMCGenDphi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCGenDphi.pdf" );
    c0->Print( "HGCPlots/hSingleMCGenDphi.png" );
    hMultiDeltaRtoBest->Draw("hist");
    hMultiDeltaRtoBest->GetXaxis()->SetTitleSize(0.05);
    hMultiDeltaRtoBest->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hMultiDeltaRtoBest.pdf" );
    c0->Print( "HGCPlots/hMultiDeltaRtoBest.png" );
    hDoubleMCEFrac->Draw("hist");
    hDoubleMCEFrac->GetXaxis()->SetTitleSize(0.05);
    hDoubleMCEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hDoubleMCEFrac.pdf" );
    c0->Print( "HGCPlots/hDoubleMCEFrac.png" );
    hDoubleMCPtFrac->Draw("hist");
    hDoubleMCPtFrac->GetXaxis()->SetTitleSize(0.05);
    hDoubleMCPtFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hDoubleMCPtFrac.pdf" );
    c0->Print( "HGCPlots/hDoubleMCPtFrac.png" );
    hSecondMCEFrac->Draw("hist");
    hSecondMCEFrac->GetXaxis()->SetTitleSize(0.05);
    hSecondMCEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSecondMCEFrac.pdf" );
    c0->Print( "HGCPlots/hSecondMCEFrac.png" );
    hSecondMCPtFrac->Draw("hist");
    hSecondMCPtFrac->GetXaxis()->SetTitleSize(0.05);
    hSecondMCPtFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSecondMCPtFrac.pdf" );
    c0->Print( "HGCPlots/hSecondMCPtFrac.png" );
    hSecondRelEta->Draw("hist");
    hSecondRelEta->GetXaxis()->SetTitleSize(0.05);
    hSecondRelEta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSecondRelEta.pdf" );
    c0->Print( "HGCPlots/hSecondRelEta.png" );
    hSecondRelPhi->Draw("hist");
    hSecondRelPhi->GetXaxis()->SetTitleSize(0.05);
    hSecondRelPhi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSecondRelPhi.pdf" );
    c0->Print( "HGCPlots/hSecondRelPhi.png" );
    hSecondRelDr->Draw("hist");
    hSecondRelDr->GetXaxis()->SetTitleSize(0.05);
    hSecondRelDr->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSecondRelDr.pdf" );
    c0->Print( "HGCPlots/hSecondRelDr.png" );
    hNumMultiClustersHighPt->Draw("hist");
    hNumMultiClustersHighPt->Draw("hist");
    hNumMultiClustersHighPt->GetXaxis()->SetTitleSize(0.05);
    hNumMultiClustersHighPt->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hNumMultiClustersHighPt.pdf" );
    c0->Print( "HGCPlots/hNumMultiClustersHighPt.png" );
    hNumSelectedMultis->Draw("hist");
    hNumSelectedMultis->GetXaxis()->SetTitleSize(0.05);
    hNumSelectedMultis->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hNumSelectedMultis.pdf" );
    c0->Print( "HGCPlots/hNumSelectedMultis.png" );
    hNumSelectedMultisContiguous->Draw("hist");
    hNumSelectedMultisContiguous->GetXaxis()->SetTitleSize(0.05);
    hNumSelectedMultisContiguous->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hNumSelectedMultisContiguous.pdf" );
    c0->Print( "HGCPlots/hNumSelectedMultisContiguous.png" );
    hEFracContiguous->Draw("hist");
    hEFracContiguous->GetXaxis()->SetTitleSize(0.05);
    hEFracContiguous->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hEFracContiguous.pdf" );
    c0->Print( "HGCPlots/hEFracContiguous.png" );
    hLayerOfClustersInSelected->Draw("hist");
    hLayerOfClustersInSelected->GetXaxis()->SetTitleSize(0.05);
    hLayerOfClustersInSelected->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hLayerOfClustersInSelected.pdf" );
    c0->Print( "HGCPlots/hLayerOfClustersInSelected.png" );
    hLayerOfClustersInSelectedWeighted->Draw("hist");
    hLayerOfClustersInSelectedWeighted->GetXaxis()->SetTitleSize(0.05);
    hLayerOfClustersInSelectedWeighted->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hLayerOfClustersInSelectedWeighted.pdf" );
    c0->Print( "HGCPlots/hLayerOfClustersInSelectedWeighted.png" );
    hLayerOfClustersInBest->Draw("hist");
    hLayerOfClustersInBest->GetXaxis()->SetTitleSize(0.05);
    hLayerOfClustersInBest->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hLayerOfClustersInBest.pdf" );
    c0->Print( "HGCPlots/hLayerOfClustersInBest.png" );
    hLayerOfClustersInBestWeighted->Draw("hist");
    hLayerOfClustersInBestWeighted->GetXaxis()->SetTitleSize(0.05);
    hLayerOfClustersInBestWeighted->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hLayerOfClustersInBestWeighted.pdf" );
    c0->Print( "HGCPlots/hLayerOfClustersInBestWeighted.png" );
    hDistanceBetweenSelected->Draw("hist");
    hDistanceBetweenSelected->GetXaxis()->SetTitle("distance / cm");
    hDistanceBetweenSelected->GetXaxis()->SetTitleSize(0.05);
    hDistanceBetweenSelected->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hDistanceBetweenSelected.pdf" );
    c0->Print( "HGCPlots/hDistanceBetweenSelected.png" );
    hDistanceContigBest->Draw("hist");
    hDistanceContigBest->GetXaxis()->SetTitle("distance / cm");
    hDistanceContigBest->GetXaxis()->SetTitleSize(0.05);
    hDistanceContigBest->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hDistanceContigBest.pdf" );
    c0->Print( "HGCPlots/hDistanceContigBest.png" );
    hMultisPassingContigPtFrac->Draw("hist");
    hMultisPassingContigPtFrac->GetXaxis()->SetTitleSize(0.05);
    hMultisPassingContigPtFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hMultisPassingContigPtFrac.pdf" );
    c0->Print( "HGCPlots/hMultisPassingContigPtFrac.png" );
    hTotRecHitEFrac->Draw("hist");
    hTotRecHitEFrac->GetXaxis()->SetTitleSize(0.05);
    hTotRecHitEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hTotRecHitEFrac.pdf" );
    c0->Print( "HGCPlots/hTotRecHitEFrac.png" );

    HistSlicer slice = HistSlicer( hSingleMCPtFracVsEta, 5 );
    auto tempVec = slice.getOneDimHists();
    slice.saveOneDimHists( tempVec, "HGCPlots/" );

    // these are the energy radius hists - use them to make TGraphs by fitting gaussians as well
    drawAndSaveHists( vecRecHitEnergyRadiusHists, c0 );
    c0->SetLogy(0);
    for( uint i=0; i<vecRecHitEnergyRadiusHists.size(); i++ ) {
      vecRecHitEnergyRadiusHists[i]->Draw("hist");
      TF1 *fGaus  = new TF1("fGaus", "gaus");
      vecRecHitEnergyRadiusHists[i]->Fit( fGaus );
      double mean = fGaus->GetParameter(1);
      double sigma = fGaus->GetParameter(2);
      //gRecHitEFracMeanVsRadius->SetPoint(  i, 1.*(i+2), mean );
      //gRecHitEFracSigmaVsRadius->SetPoint( i, 1.*(i+2), sigma );
      gRecHitEFracMeanVsRadius->SetPoint(  i, 0.5*(i+2), mean );
      gRecHitEFracSigmaVsRadius->SetPoint( i, 0.5*(i+2), sigma );
    }
    gRecHitEFracMeanVsRadius->Draw("AC*");
    gRecHitEFracMeanVsRadius->SetTitle("Mean fraction of gen energy contained in recHits within r cm");
    gRecHitEFracMeanVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracMeanVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracMeanVsRadius.png");
    gRecHitEFracSigmaVsRadius->Draw("AC*");
    gRecHitEFracSigmaVsRadius->SetTitle("Sigma fraction of gen energy contained in recHits within r cm");
    gRecHitEFracSigmaVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracSigmaVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracSigmaVsRadius.png");

    drawAndSaveHists( vecRecHitEnergyBestHists, c0 );
    c0->SetLogy(0);
    /*for( uint i=0; i<vecRecHitEnergyBestHists.size(); i++ ) {
      vecRecHitEnergyBestHists[i]->Draw("hist");
      TF1 *fGaus  = new TF1("fGaus", "gaus");
      vecRecHitEnergyBestHists[i]->Fit( fGaus );
      double mean = fGaus->GetParameter(1);
      double sigma = fGaus->GetParameter(2);
      //gRecHitEFracMeanBestVsRadius->SetPoint(  i, 1.*(i+2), mean );
      //gRecHitEFracSigmaBestVsRadius->SetPoint( i, 1.*(i+2), sigma );
      //if( sigma > 0 ) gRecHitEFracResolutionBestVsRadius->SetPoint( i, 1.*(i+2), sigma/mean );
      gRecHitEFracMeanBestVsRadius->SetPoint(  i, 0.5*(i+2), mean );
      gRecHitEFracSigmaBestVsRadius->SetPoint( i, 0.5*(i+2), sigma );
      if( mean != 0. ) gRecHitEFracResolutionBestVsRadius->SetPoint( i, 0.5*(i+2), sigma/mean );
    }*/
    gRecHitEFracMeanBestVsRadius->Draw("AC*");
    gRecHitEFracMeanBestVsRadius->SetTitle("Mean fraction of gen energy contained in recHits within r cm");
    gRecHitEFracMeanBestVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracMeanBestVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracMeanBestVsRadius.png");
    gRecHitEFracSigmaBestVsRadius->Draw("AC*");
    gRecHitEFracSigmaBestVsRadius->SetTitle("Sigma fraction of gen energy contained in recHits within r cm");
    gRecHitEFracSigmaBestVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracSigmaBestVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracSigmaBestVsRadius.png");
    gRecHitEFracResolutionBestVsRadius->Draw("AC*");
    gRecHitEFracResolutionBestVsRadius->SetTitle("Resolution of gen energy contained in recHits within r cm");
    gRecHitEFracResolutionBestVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracResolutionBestVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracResolutionBestVsRadius.png");

    drawAndSaveHists( vecRecHitEnergySelectedHists, c0 );
    c0->SetLogy(0);
    /*for( uint i=0; i<vecRecHitEnergySelectedHists.size(); i++ ) {
      vecRecHitEnergySelectedHists[i]->Draw("hist");
      TF1 *fGaus  = new TF1("fGaus", "gaus");
      vecRecHitEnergySelectedHists[i]->Fit( fGaus );
      double mean = fGaus->GetParameter(1);
      double sigma = fGaus->GetParameter(2);
      //gRecHitEFracMeanSelectedVsRadius->SetPoint(  i, 1.*(i+2), mean );
      //gRecHitEFracSigmaSelectedVsRadius->SetPoint( i, 1.*(i+2), sigma );
      //if( sigma > 0 ) gRecHitEFracResolutionSelectedVsRadius->SetPoint( i, 1.*(i+2), sigma/mean );
      gRecHitEFracMeanSelectedVsRadius->SetPoint(  i, 0.5*(i+2), mean );
      gRecHitEFracSigmaSelectedVsRadius->SetPoint( i, 0.5*(i+2), sigma );
      if( mean != 0 ) gRecHitEFracResolutionSelectedVsRadius->SetPoint( i, 0.5*(i+2), sigma/mean );
    }*/
    gRecHitEFracMeanSelectedVsRadius->Draw("AC*");
    gRecHitEFracMeanSelectedVsRadius->SetTitle("Mean fraction of gen energy contained in recHits within r cm");
    gRecHitEFracMeanSelectedVsRadius->GetXaxis()->SetTitle("r / cm");
    string selectName = "HGCPlots/gRecHitEFracMeanSelected" + NumberToString(numClusterCut) + "VsRadius.pdf";
    c0->Print( selectName.c_str() );
    selectName = "HGCPlots/gRecHitEFracMeanSelected" + NumberToString(numClusterCut) + "VsRadius.png";
    c0->Print( selectName.c_str() );
    gRecHitEFracSigmaSelectedVsRadius->Draw("AC*");
    gRecHitEFracSigmaSelectedVsRadius->SetTitle("Sigma fraction of gen energy contained in recHits within r cm");
    gRecHitEFracSigmaSelectedVsRadius->GetXaxis()->SetTitle("r / cm");
    selectName = "HGCPlots/gRecHitEFracSigmaSelected" + NumberToString(numClusterCut) + "VsRadius.pdf";
    c0->Print( selectName.c_str() );
    selectName = "HGCPlots/gRecHitEFracSigmaSelected" + NumberToString(numClusterCut) + "VsRadius.png";
    c0->Print( selectName.c_str() );
    gRecHitEFracResolutionSelectedVsRadius->Draw("AC*");
    gRecHitEFracResolutionSelectedVsRadius->SetTitle("Resolution of gen energy contained in recHits within r cm");
    gRecHitEFracResolutionSelectedVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracResolutionSelectedVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracResolutionSelectedVsRadius.png");

    drawAndSaveHists( vecRecHitEnergyContiguousHists, c0 );
    c0->SetLogy(0);
    for( uint i=0; i<vecRecHitEnergyContiguousHists.size(); i++ ) {
      vecRecHitEnergyContiguousHists[i]->Draw("hist");
      TF1 *fGaus  = new TF1("fGaus", "gaus");
      vecRecHitEnergyContiguousHists[i]->Fit( fGaus );
      double mean = fGaus->GetParameter(1);
      double sigma = fGaus->GetParameter(2);
      //gRecHitEFracMeanContiguousVsRadius->SetPoint(  i, 1.*(i+2), mean );
      //gRecHitEFracSigmaContiguousVsRadius->SetPoint( i, 1.*(i+2), sigma );
      //if( sigma > 0 ) gRecHitEFracResolutionContiguousVsRadius->SetPoint( i, 1.*(i+2), sigma/mean );
      gRecHitEFracMeanContiguousVsRadius->SetPoint(  i, 0.5*(i+2), mean );
      gRecHitEFracSigmaContiguousVsRadius->SetPoint( i, 0.5*(i+2), sigma );
      if( mean != 0 ) gRecHitEFracResolutionContiguousVsRadius->SetPoint( i, 0.5*(i+2), sigma/mean );
    }
    gRecHitEFracMeanContiguousVsRadius->Draw("AC*");
    gRecHitEFracMeanContiguousVsRadius->SetTitle("Mean fraction of gen energy contained in recHits within r cm");
    gRecHitEFracMeanContiguousVsRadius->GetXaxis()->SetTitle("r / cm");
    selectName = "HGCPlots/gRecHitEFracMeanContiguous" + NumberToString(numClusterCut) + "VsRadius.pdf";
    c0->Print( selectName.c_str() );
    selectName = "HGCPlots/gRecHitEFracMeanContiguous" + NumberToString(numClusterCut) + "VsRadius.png";
    c0->Print( selectName.c_str() );
    gRecHitEFracSigmaContiguousVsRadius->Draw("AC*");
    gRecHitEFracSigmaContiguousVsRadius->SetTitle("Sigma fraction of gen energy contained in recHits within r cm");
    gRecHitEFracSigmaContiguousVsRadius->GetXaxis()->SetTitle("r / cm");
    selectName = "HGCPlots/gRecHitEFracSigmaContiguous" + NumberToString(numClusterCut) + "VsRadius.pdf";
    c0->Print( selectName.c_str() );
    selectName = "HGCPlots/gRecHitEFracSigmaContiguous" + NumberToString(numClusterCut) + "VsRadius.png";
    c0->Print( selectName.c_str() );
    gRecHitEFracResolutionContiguousVsRadius->Draw("AC*");
    gRecHitEFracResolutionContiguousVsRadius->SetTitle("Resolution of gen energy contained in recHits within r cm");
    gRecHitEFracResolutionContiguousVsRadius->GetXaxis()->SetTitle("r / cm");
    c0->Print("HGCPlots/gRecHitEFracResolutionContiguousVsRadius.pdf");
    c0->Print("HGCPlots/gRecHitEFracResolutionContiguousVsRadius.png");

    gConvRatio->Draw("AP");
    gConvRatio->SetTitle("Conversion ratio vs eta");
    gConvRatio->GetXaxis()->SetTitle("#eta");
    c0->Print("HGCPlots/gConvRatio.pdf");
    c0->Print("HGCPlots/gConvRatio.png");

    gExtraRatio->Draw("AP");
    gExtraRatio->SetTitle("Fraction of events with >1 selected multicluster vs eta");
    gExtraRatio->GetXaxis()->SetTitle("#eta");
    c0->Print("HGCPlots/gExtraRatio.pdf");
    c0->Print("HGCPlots/gExtraRatio.png");

    gContigRatio->Draw("AP");
    gContigRatio->SetTitle("Fraction of events with >1 multicluster passing contiguity cut vs eta");
    gContigRatio->GetXaxis()->SetTitle("#eta");
    c0->Print("HGCPlots/gContigRatio.pdf");
    c0->Print("HGCPlots/gContigRatio.png");

    drawAndSaveHists( vecNumSelectedMultiHists, c0 );
    drawAndSaveHists( vecSelectedMultiEnergyHists, c0 );
    for( uint i=2; i<22; i++) {
      //gNumSelectedVsCut->SetPoint( i-2, i, vecNumSelectedMultiHists[i]->GetMean() );
      //gSelectedEnergyVsCut->SetPoint( i-2, i, vecSelectedMultiEnergyHists[i]->GetMean() );
      gNumSelectedVsCut->SetPoint( i-2, i+1, vecNumSelectedMultiHists[i]->GetMean() );
      gSelectedEnergyVsCut->SetPoint( i-2, i+1, vecSelectedMultiEnergyHists[i]->GetMean() );
    }
    c0->SetLogy(0);
    gNumSelectedVsCut->Draw("AC*");
    gNumSelectedVsCut->SetTitle("Mean number of selected multiclusters vs required num 2Ds");
    gNumSelectedVsCut->GetXaxis()->SetTitle("2D clusters required");
    gNumSelectedVsCut->GetYaxis()->SetTitle("Mean selected multis");
    c0->Print("HGCPlots/gNumSelectedVsCut.pdf");
    c0->Print("HGCPlots/gNumSelectedVsCut.png");
    gSelectedEnergyVsCut->Draw("AC*");
    gSelectedEnergyVsCut->SetTitle("Mean total energy of selected multiclusters vs required num 2Ds");
    gSelectedEnergyVsCut->GetXaxis()->SetTitle("2D clusters required");
    gSelectedEnergyVsCut->GetYaxis()->SetTitle("Mean total energy");
    c0->Print("HGCPlots/gSelectedEnergyVsCut.pdf");
    c0->Print("HGCPlots/gSelectedEnergyVsCut.png");
    hRecHitPhiMean->Draw("hist");
    hRecHitPhiMean->SetTitle("Mean rechit gen dphi, per event");
    hRecHitPhiMean->GetXaxis()->SetTitle("mean dphi");
    c0->Print("HGCPlots/hRecHitPhiMean.pdf");
    c0->Print("HGCPlots/hRecHitPhiMean.png");
    hRecHitPhiRms->Draw("hist");
    hRecHitPhiRms->SetTitle("rms rechit gen dphi, per event");
    hRecHitPhiRms->GetXaxis()->SetTitle("rms dphi");
    c0->Print("HGCPlots/hRecHitPhiRms.pdf");
    c0->Print("HGCPlots/hRecHitPhiRms.png");
    hRecHitEtaMean->Draw("hist");
    hRecHitEtaMean->SetTitle("Mean rechit gen deta, per event");
    hRecHitEtaMean->GetXaxis()->SetTitle("mean deta");
    c0->Print("HGCPlots/hRecHitEtaMean.pdf");
    c0->Print("HGCPlots/hRecHitEtaMean.png");
    hRecHitEtaRms->Draw("hist");
    hRecHitEtaRms->SetTitle("rms rechit gen deta, per event");
    hRecHitEtaRms->GetXaxis()->SetTitle("rms deta");
    c0->Print("HGCPlots/hRecHitEtaRms.pdf");
    c0->Print("HGCPlots/hRecHitEtaRms.png");

    gStyle->SetPalette( 57 );
    gStyle->SetNumberContours( 255 );
    gStyle->SetOptStat( 0 );
    c0->SetLogy(0);
    hMultiClusterDeltaEtaDeltaPhi->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhi->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhi->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhi->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhi.pdf" );
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhi.png" );
    hMultiClusterDeltaEtaDeltaPhiWeighted->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhiWeighted->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhiWeighted->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhiWeighted->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhiWeighted->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeighted.pdf");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeighted.png");
    hMultiClusterDeltaEtaDeltaPhiZoomed->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhiZoomed->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhiZoomed->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhiZoomed->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhiZoomed->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhiZoomed.pdf" );
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhiZoomed.png" );
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeightedZoomed.pdf");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeightedZoomed.png");
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhiZoomZoom.pdf" );
    c0->Print( "HGCPlots/hMultiClusterDeltaEtaDeltaPhiZoomZoom.png" );
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->Draw("colz");
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->GetYaxis()->SetTitle("#Delta#eta");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom.pdf");
    c0->Print("HGCPlots/hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom.png");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->Draw("colz");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->GetYaxis()->SetTitle("#Delta#eta");
    string selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoom.pdf";
    c0->Print(selectedDeltaName.c_str());
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoom.png";
    c0->Print(selectedDeltaName.c_str());
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->Draw("colz");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->GetYaxis()->SetTitle("#Delta#eta");
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWide.pdf";
    c0->Print(selectedDeltaName.c_str());
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWide.png";
    c0->Print(selectedDeltaName.c_str());
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->Draw("colz");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->GetYaxis()->SetTitle("#Delta#eta");
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWeighted.pdf";
    c0->Print(selectedDeltaName.c_str());
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWeighted.png";
    c0->Print(selectedDeltaName.c_str());
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->Draw("colz");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->GetXaxis()->SetTitle("#Delta#phi");
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->GetXaxis()->SetTitleSize(0.05);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->GetXaxis()->SetTitleOffset(0.8);
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->GetYaxis()->SetTitle("#Delta#eta");
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWideWeighted.pdf";
    c0->Print(selectedDeltaName.c_str());
    selectedDeltaName = "HGCPlots/hMultiClusterSelected" + NumberToString(numClusterCut) + "DeltaEtaDeltaPhiZoomZoomWideWeighted.png";
    c0->Print(selectedDeltaName.c_str());
    hNumClusters2DPerMultiVsPt->Draw("colz");
    hNumClusters2DPerMultiVsPt->GetXaxis()->SetTitle("Multicluster p_{T}");
    hNumClusters2DPerMultiVsPt->GetXaxis()->SetTitleSize(0.05);
    hNumClusters2DPerMultiVsPt->GetXaxis()->SetTitleOffset(0.8);
    hNumClusters2DPerMultiVsPt->GetYaxis()->SetTitle("Number of constituent 2D clusters");
    c0->Print("HGCPlots/hNumClusters2DPerMultiVsPt.pdf");
    c0->Print("HGCPlots/hNumClusters2DPerMultiVsPt.png");
    hNumClusters2DPerBestMultiVsPt->Draw("colz");
    hNumClusters2DPerBestMultiVsPt->GetXaxis()->SetTitle("Multicluster p_{T}");
    hNumClusters2DPerBestMultiVsPt->GetXaxis()->SetTitleSize(0.05);
    hNumClusters2DPerBestMultiVsPt->GetXaxis()->SetTitleOffset(0.8);
    hNumClusters2DPerBestMultiVsPt->GetYaxis()->SetTitle("Number of constituent 2D clusters");
    c0->Print("HGCPlots/hNumClusters2DPerBestMultiVsPt.pdf");
    c0->Print("HGCPlots/hNumClusters2DPerBestMultiVsPt.png");
    hNumClusters2DPerSecondBestMultiVsPt->Draw("colz");
    hNumClusters2DPerSecondBestMultiVsPt->GetXaxis()->SetTitle("Multicluster p_{T}");
    hNumClusters2DPerSecondBestMultiVsPt->GetXaxis()->SetTitleSize(0.05);
    hNumClusters2DPerSecondBestMultiVsPt->GetXaxis()->SetTitleOffset(0.8);
    hNumClusters2DPerSecondBestMultiVsPt->GetYaxis()->SetTitle("Number of constituent 2D clusters");
    c0->Print("HGCPlots/hNumClusters2DPerSecondBestMultiVsPt.pdf");
    c0->Print("HGCPlots/hNumClusters2DPerSecondBestMultiVsPt.png");
    hSingleMCPtFracVsEta->Draw("colz");
    hSingleMCPtFracVsEta->GetXaxis()->SetTitleSize(0.05);
    hSingleMCPtFracVsEta->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCPtFracVsEta.pdf" );
    c0->Print( "HGCPlots/hSingleMCPtFracVsEta.png" );
    hSingleMCGenDetaDphi->Draw("colz");
    hSingleMCGenDetaDphi->GetXaxis()->SetTitleSize(0.05);
    hSingleMCGenDetaDphi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCGenDetaDphi.pdf" );
    c0->Print( "HGCPlots/hSingleMCGenDetaDphi.png" );
    hSingleMCGenDetaDphiWeighted->Draw("colz");
    hSingleMCGenDetaDphiWeighted->GetXaxis()->SetTitleSize(0.05);
    hSingleMCGenDetaDphiWeighted->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCGenDetaDphiWeighted.pdf" );
    c0->Print( "HGCPlots/hSingleMCGenDetaDphiWeighted.png" );
    hSingleMCEFracVsNumSelected->Draw("colz");
    hSingleMCEFracVsNumSelected->GetXaxis()->SetTitleSize(0.05);
    hSingleMCEFracVsNumSelected->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hSingleMCEFracVsNumSelected.pdf" );
    c0->Print( "HGCPlots/hSingleMCEFracVsNumSelected.png" );
    hRecHitEtaPhi->Draw("colz");
    hRecHitEtaPhi->GetXaxis()->SetTitleSize(0.05);
    hRecHitEtaPhi->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hRecHitEtaPhi.pdf" );
    c0->Print( "HGCPlots/hRecHitEtaPhi.png" );
    hRecHitEtaPhiZoom->Draw("colz");
    hRecHitEtaPhiZoom->GetXaxis()->SetTitleSize(0.05);
    hRecHitEtaPhiZoom->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hRecHitEtaPhiZoom.pdf" );
    c0->Print( "HGCPlots/hRecHitEtaPhiZoom.png" );
    hRecHitEtaPhiZoom->Draw("colz");
    hRecHitEtaPhiZoom->GetXaxis()->SetTitleSize(0.05);
    hRecHitEtaPhiZoom->GetXaxis()->SetTitleOffset(0.8);
    c0->Print( "HGCPlots/hRecHitEtaPhiZoom.pdf" );
    c0->Print( "HGCPlots/hRecHitEtaPhiZoom.png" );
    hRecHitLayerVsClusterZ->Draw("colz");
    hRecHitLayerVsClusterZ->SetTitle("Rec hit layer vs cluster z");
    hRecHitLayerVsClusterZ->GetXaxis()->SetTitle("rec hit layer");
    hRecHitLayerVsClusterZ->GetYaxis()->SetTitle("cluster z");
    c0->Print("HGCPlots/hRecHitLayerVsClusterZ.pdf");
    c0->Print("HGCPlots/hRecHitLayerVsClusterZ.png");

    cout << "entering pretty bit" << endl << endl;
    hSingleMCEFrac->Draw();
    gStyle->SetOptStat( 0 );
    hSingleMCEFrac->SetTitle("Fraction of gen energy in \"best\" multicluster");
    hSingleMCEFrac->GetXaxis()->SetTitle("E_{multi} / E_{gen}");
    hSingleMCEFrac->GetXaxis()->SetTitleSize(0.05);
    hSingleMCEFrac->GetXaxis()->SetTitleOffset(0.8);
    c0->Print("HGCPlots/PrettyPlot1.pdf");
    c0->Print("HGCPlots/PrettyPlot1.png");
    hSingleMCEFrac->Fit("gaus");
    gStyle->SetOptStat( 0 );
    c0->Print("HGCPlots/PrettyPlot1a.pdf");
    c0->Print("HGCPlots/PrettyPlot1a.png");
    vecRecHitEnergyBestHists[2]->Draw();
    gStyle->SetOptStat( 0 );
    vecRecHitEnergyBestHists[2]->SetTitle("Fraction of gen energy in RecHits within 2cm of \"best\" multicluster");
    vecRecHitEnergyBestHists[2]->GetXaxis()->SetTitle("E_{#Sigma RecHits} / E_{gen}");
    vecRecHitEnergyBestHists[2]->GetXaxis()->SetTitleSize(0.05);
    vecRecHitEnergyBestHists[2]->GetXaxis()->SetTitleOffset(0.8);
    c0->Print("HGCPlots/PrettyPlot2.pdf");
    c0->Print("HGCPlots/PrettyPlot2.png");
    vecRecHitEnergyBestHists[2]->Fit("gaus");
    gStyle->SetOptStat( 0 );
    c0->Print("HGCPlots/PrettyPlot2a.pdf");
    c0->Print("HGCPlots/PrettyPlot2a.png");
    vecRecHitEnergySelectedHists[8]->Draw();
    gStyle->SetOptStat( 0 );
    vecRecHitEnergySelectedHists[8]->SetTitle("Energy in RecHits within 5cm of selected multiclusters");
    vecRecHitEnergySelectedHists[8]->GetXaxis()->SetTitle("E_{#Sigma RecHits} / E_{gen}");
    vecRecHitEnergySelectedHists[8]->GetXaxis()->SetTitleSize(0.05);
    vecRecHitEnergySelectedHists[8]->GetXaxis()->SetTitleOffset(0.8);
    c0->Print("HGCPlots/PrettyPlot3.pdf");
    c0->Print("HGCPlots/PrettyPlot3.png");
    vecRecHitEnergySelectedHists[8]->Fit("gaus");
    gStyle->SetOptStat( 0 );
    c0->Print("HGCPlots/PrettyPlot3a.pdf");
    c0->Print("HGCPlots/PrettyPlot3a.png");
    vecSelectedMultiEnergyHists[2]->Draw();
    gStyle->SetOptStat( 0 );
    vecSelectedMultiEnergyHists[2]->SetTitle("Fraction of gen energy in selected (nClus>2) multiclusters");
    vecSelectedMultiEnergyHists[2]->GetXaxis()->SetTitle("E_{#Sigma Multis} / E_{gen}");
    vecSelectedMultiEnergyHists[2]->GetXaxis()->SetTitleSize(0.05);
    vecSelectedMultiEnergyHists[2]->GetXaxis()->SetTitleOffset(0.8);
    c0->Print("HGCPlots/PrettyPlot4.pdf");
    c0->Print("HGCPlots/PrettyPlot4.png");
    vecSelectedMultiEnergyHists[2]->Fit("gaus");
    gStyle->SetOptStat( 0 );
    c0->Print("HGCPlots/PrettyPlot4a.pdf");
    c0->Print("HGCPlots/PrettyPlot4a.png");
  }


  // write histograms
  if( write ) {
    string outputName = "Output/output_Photon_Pt35.root";
    if( sampleType != "" ) outputName = "Output/output_" + sampleType + "_" + convType + "_" + extType + ".root";
    TFile* outputFile = new TFile( outputName.c_str(), "recreate" );

    hGenEta->Write();
    hGenPhi->Write();
    hGenEnergy->Write();
    hGenPt->Write();
    hSingleMCEFrac->Write();
    hSingleMCLimitedEFrac->Write();
    hTotalClusterEFrac->Write();
    hBHclusterEta->Write();
    hBHclusterPhi->Write();
    hBHrecHitEta->Write();
    hBHrecHitPhi->Write();
    hBestPfEFrac->Write();
    hSingleMCPtFrac->Write();
    hSingleMCGenDr->Write();
    hSingleMCGenDeta->Write();
    hSingleMCGenDphi->Write();
    hDoubleMCEFrac->Write();
    hDoubleMCPtFrac->Write();
    hSecondMCEFrac->Write();
    hSecondMCPtFrac->Write();
    hSecondRelEta->Write();
    hSecondRelPhi->Write();
    hSecondRelDr->Write();
    hNumMultiClustersHighPt->Write();
    hMultiClusterDeltaEtaDeltaPhi->Write();
    hMultiClusterDeltaEtaDeltaPhiWeighted->Write();
    hMultiClusterDeltaEtaDeltaPhiZoomed->Write();
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomed->Write();
    hMultiClusterDeltaEtaDeltaPhiZoomZoom->Write();
    hMultiClusterDeltaEtaDeltaPhiWeightedZoomZoom->Write();
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoom->Write();
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWeighted->Write();
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWide->Write();
    hMultiClusterSelectedDeltaEtaDeltaPhiZoomZoomWideWeighted->Write();
    hNumClusters2DPerMultiVsPt->Write();
    hNumClusters2DPerBestMultiVsPt->Write();
    hNumClusters2DPerSecondBestMultiVsPt->Write();
    hSingleMCPtFracVsEta->Write();
    hSingleMCEFracVsNumSelected->Write();
    hLayerOfClustersInSelected->Write();
    hLayerOfClustersInSelectedWeighted->Write();
    hLayerOfClustersInBest->Write();
    hLayerOfClustersInBestWeighted->Write();
    hDistanceBetweenSelected->Write();
    hDistanceContigBest->Write();
    hMultisPassingContigPtFrac->Write();
    hNumSelectedMultis->Write();
    hNumSelectedMultisContiguous->Write();
    hEFracContiguous->Write();
    hMultiDeltaRtoBest->Write();
    hPseudoSuperClusterEFrac->Write();

    gRecHitEFracMeanVsRadius->Write();
    gRecHitEFracSigmaVsRadius->Write();
    gRecHitEFracMeanBestVsRadius->Write();
    gRecHitEFracSigmaBestVsRadius->Write();
    gRecHitEFracResolutionBestVsRadius->Write();
    gRecHitEFracMeanSelectedVsRadius->Write();
    gRecHitEFracSigmaSelectedVsRadius->Write();
    gRecHitEFracResolutionSelectedVsRadius->Write();
    gConvRatio->Write();
    gExtraRatio->Write();
    gNumSelectedVsCut->Write();
    gSelectedEnergyVsCut->Write();

    writeHists( vecRecHitEnergyRadiusHists );
    writeHists( vecRecHitEnergyBestHists );
    writeHists( vecRecHitEnergySelectedHists );
    writeHists( vecRecHitEnergyContiguousHists );

    outputFile->Close();
  }


  return 0;
}