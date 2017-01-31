// Takes a 2D histogram and separates it into 1D histograms
// Ed Scott, 11.05.2016

#ifndef HistSlicer_h
#define HistSlicer_h

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
#include "TCanvas.h"
#include "NumberToString.h"
using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::vector;

class HistSlicer
{
    public:
        HistSlicer()  {}
        ~HistSlicer() {}

        HistSlicer( TH2* twoDimHist ) { 
            twoDimHist_ = twoDimHist;
            splitNum_ = 1;
        }

        HistSlicer( TH2* twoDimHist, int splitNum )
        {
            twoDimHist_ = twoDimHist;
            splitNum_ = splitNum;
        }
        
        void setSplitNum(   int splitNum     ) { splitNum_   = splitNum;   }
        void setTwoDimHist( TH2* twoDimHist ) { twoDimHist_ = twoDimHist; }

        int   getSplitNum  () const { return splitNum_;   }
        TH2* getTwoDimHist() const { return twoDimHist_; }

        vector<TH1D*> getOneDimHists() const
        {
            vector<TH1D*> vecOneDimHists;

            int nXBins      = twoDimHist_->GetNbinsX();
            //cout << "nXBins = " << nXBins << endl;
            int remainder   = nXBins % splitNum_;
            //cout << "remainder = " << remainder << endl;
            int newNumBins  = nXBins / splitNum_;
            //cout << "newNumBins = " << newNumBins << endl;
            float xMin   = twoDimHist_->GetXaxis()->GetXmin();
            //cout << "xMin = " << xMin << endl;
            float xMax   = twoDimHist_->GetXaxis()->GetXmax();
            //cout << "xMax = " << xMax << endl;
            float xWidth = (xMax - xMin) / nXBins;
            //cout << "xWidth = " << xWidth << endl;

            int lowXBin  = 0;
            int highXBin = 0;
            for( int histIndex = 0; histIndex < splitNum_; histIndex++ ) {
                string histNumber = NumberToString( histIndex + 1 );
                string histName(  twoDimHist_->GetName() );
                histName += "_oneDimHist" + NumberToString( histIndex + 1 );

                lowXBin  = highXBin + 1;
                highXBin = lowXBin  + newNumBins - 1;
                if( histIndex < remainder ) highXBin++;
                //cout << "lowXBin  = " << lowXBin  << endl;
                //cout << "highXBin = " << highXBin << endl;
                vecOneDimHists.push_back( twoDimHist_->ProjectionY( histName.c_str(), lowXBin, highXBin ) );

                string oldXTitle( twoDimHist_->GetXaxis()->GetTitle() );
                string oldYTitle( twoDimHist_->GetYaxis()->GetTitle() );
                string histTitle = "Projection " + histNumber + " / " + NumberToString( splitNum_ ) + " of " + oldYTitle + " vs " + oldXTitle;
                vecOneDimHists[histIndex]->SetTitle( histTitle.c_str() );
                float xLow  = xMin + ( lowXBin - 1 ) * xWidth;
                //cout << "xLow = " << xLow << endl;
                float xHigh = xMin + ( highXBin    ) * xWidth;
                //cout << "xHigh = " << xHigh << endl;
                string histYTitle = oldXTitle + " from " + NumberToString( xLow ) + " to " + NumberToString( xHigh );
                vecOneDimHists[histIndex]->GetYaxis()->SetTitle( histYTitle.c_str() );
            }
            return vecOneDimHists;
        }

        vector<TH1D*> getOneDimHists( vector<float> &edges ) const
        {
            vector<TH1D*> vecOneDimHists;

            string oldXTitle( twoDimHist_->GetXaxis()->GetTitle() );
            string oldYTitle( twoDimHist_->GetYaxis()->GetTitle() );

            int nXBins      = twoDimHist_->GetNbinsX();
            //cout << "nXBins = " << nXBins << endl;
            float xMin   = twoDimHist_->GetXaxis()->GetXmin();
            //cout << "xMin = " << xMin << endl;
            float xMax   = twoDimHist_->GetXaxis()->GetXmax();
            //cout << "xMax = " << xMax << endl;
            float xWidth = (xMax - xMin) / nXBins;
            //cout << "xWidth = " << xWidth << endl;

            for( int histIndex = 0; histIndex < edges.size() - 1; histIndex++ ) {
                string histNumber = NumberToString( histIndex + 1 );
                string histName(  twoDimHist_->GetName() );
                histName += "_oneDimHist" + NumberToString( histIndex + 1 );

                int lowXBin  = floor( ( edges[histIndex]     - xMin ) / xWidth ) + 1;
                int highXBin = floor( ( edges[histIndex + 1] - xMin ) / xWidth )    ;
                //cout << "lowXBin  = " << lowXBin  << endl;
                //cout << "highXBin = " << highXBin << endl;
                vecOneDimHists.push_back( twoDimHist_->ProjectionY( histName.c_str(), lowXBin, highXBin ) );

                string histTitle = "Projection " + histNumber + " / " + NumberToString( edges.size() - 1 ) + " of " + oldYTitle + " vs " + oldXTitle;
                vecOneDimHists[histIndex]->SetTitle( histTitle.c_str() );
                float xLow  = xMin + ( lowXBin  - 1 ) * xWidth;
                //cout << "xLow = " << xLow << endl;
                float xHigh = xMin + ( highXBin     ) * xWidth;
                //cout << "xHigh = " << xHigh << endl;
                string histYTitle = oldXTitle + " from " + NumberToString( xLow ) + " to " + NumberToString( xHigh );
                vecOneDimHists[histIndex]->GetYaxis()->SetTitle( histYTitle.c_str() );
            }
            return vecOneDimHists;
        }

        void saveOneDimHists( vector<TH1D*> &v ) const {
            saveOneDimHists( v, "" );
        }

        void saveOneDimHists( vector<TH1D*> &v, string loc ) const {
            TCanvas *c = new TCanvas();
            for( int i = 0; i < v.size(); i++ ) {
                //v[i]->Draw();
                v[i]->Draw("hist");
                string baseName( v[i]->GetName() );
                string printName = loc + baseName + ".png";
                c->Print( printName.c_str() );
                printName = loc + baseName + ".pdf";
                c->Print( printName.c_str() );
            }
            delete c;
        }

        void writeOneDimHists( vector<TH1D*> &v ) const {
            writeOneDimHists( v, "" );
        }

        void writeOneDimHists( vector<TH1D*> &v, string loc ) const {
            for( int i = 0; i < v.size(); i++ ) {
                v[i]->Write( loc.c_str() );
            }
        }

    private:
        TH2 * twoDimHist_;
        int    splitNum_; 
};

#endif

// Local Variables:
// // mode:c++
// // indent-tabs-mode:nil
// // tab-width:4
// // c-basic-offset:4
// // End:
// // vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
