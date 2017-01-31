// Simple templated function to convert a number to a string

#ifndef NumberToString_h
#define NumberToString_h

#include <iostream>
#include <sstream>
#include <string>
using std::ostringstream;
using std::string;

template <typename T>
  string NumberToString( T Number )
  {
    ostringstream ss;
    ss << Number;
    return ss.str();
  }

#endif

// Local Variables:
// // mode:c++
// // indent-tabs-mode:nil
// // tab-width:4
// // c-basic-offset:4
// // End:
// // vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
