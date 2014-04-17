#include<iostream>
#include<TROOT.h>
using namespace std;

void setup(){
   // you can find the fastjet locations by doing:
   // $ scram setup fastjet
  
    
   
 
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/BaconAna/DataFormats/interface/");
   gSystem->AddIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/include");
   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
   gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/lib/libfastjet.so");
   gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/lib/libfastjettools.so");   
   gSystem->Load("libFWCoreFWLite.so");
   AutoLibraryLoader::enable(); 
   gSystem->Load("libDataFormatsFWLite.so");





  








}
