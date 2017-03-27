// code that automatically makes, fills and writes hists when there are multiple categories
// currently configured for the eta and energy ratio categories - will generalise later
vector<TH1F*> MakeHistsByCategory( string name, string title, int nBins, float xLow, float xHigh, string xTitle = "", string yTitle = "" ) {
  vector<TH1F*> v;
  vector<string> cats;
  cats.push_back( "LowEB");
  cats.push_back( "HighEB");
  cats.push_back( "LowEE");
  cats.push_back( "HighEE");
  for( int i = 0; i < cats.size(); i++ ) {
    string fullname  = name  + cats[i];
    string fulltitle = title + ", " + cats[i];
    TH1F* h = new TH1F( fullname.c_str(), fulltitle.c_str(), nBins, xLow, xHigh );
    h->GetXaxis()->SetTitle( xTitle.c_str() );
    h->GetYaxis()->SetTitle( yTitle.c_str() );
    v.push_back( h );
    //v[i]->GetXaxis()->SetTitle( xTitle.c_str() );
    //v[i]->GetYaxis()->SetTitle( yTitle.c_str() );
    //delete h;
  }
  return v;
}

vector<TH2F*> MakeHistsByCategory( string name, string title, int nXBins, float xLow, float xHigh, int nYBins, float yLow, float yHigh, string xTitle = "", string yTitle = "" ) {
  vector<TH2F*> v;
  vector<string> cats;
  cats.push_back( "LowEB");
  cats.push_back( "HighEB");
  cats.push_back( "LowEE");
  cats.push_back( "HighEE");
  for( int i = 0; i < cats.size(); i++ ) {
    string fullname  = name  + cats[i];
    string fulltitle = title + ", " + cats[i];
    TH2F* h = new TH2F( fullname.c_str(), fulltitle.c_str(), nXBins, xLow, xHigh, nYBins, yLow, yHigh );
    h->GetXaxis()->SetTitle( xTitle.c_str() );
    h->GetYaxis()->SetTitle( yTitle.c_str() );
    v.push_back( h );
    //v[i]->GetXaxis()->SetTitle( xTitle.c_str() );
    //v[i]->GetYaxis()->SetTitle( yTitle.c_str() );
    //delete h;
  }
  return v;
}

template <typename T>
  void FillHistsByCategory( vector<TH1F*> &v, float eta, float ratio, T val, float wei = 1. ) {
    if(      abs( eta ) < 1.5 && ratio < 0.8) v[0]->Fill( val, wei );
    else if( abs( eta ) < 1.5 && ratio < 1.2) v[1]->Fill( val, wei );
    else if( abs( eta ) < 2.5 && ratio < 0.8) v[2]->Fill( val, wei );
    else if( abs( eta ) < 2.5 && ratio < 1.2) v[3]->Fill( val, wei );
  }

void FillHistsByCategory( vector<TH2F*> &v, float eta, float ratio, float xVal, float yVal, float wei = 1. ) {
  if(      abs( eta ) < 1.5 && ratio < 0.8) v[0]->Fill( xVal, yVal, wei );
  else if( abs( eta ) < 1.5 && ratio < 1.2) v[1]->Fill( xVal, yVal, wei );
  else if( abs( eta ) < 2.5 && ratio < 0.8) v[2]->Fill( xVal, yVal, wei );
  else if( abs( eta ) < 2.5 && ratio < 1.2) v[3]->Fill( xVal, yVal, wei );
}

void WriteHistsByCategory( vector<TH1F*> &v ) {
  for( int i = 0; i < v.size(); i++ ) {
    v[i]->Write();
  }
}

void WriteHistsByCategory( vector<TH2F*> &v ) {
  for( int i = 0; i < v.size(); i++ ) {
    v[i]->Write();
  }
}
