#include<iostream>
#include<TROOT.h>
using namespace std;

void setup(){
   // you can find the fastjet locations by doing:
   // $ scram setup fastjet
  
    
   
 
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/BaconAna/DataFormats/interface/");
   gSystem->AddIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/include");
   gSystem->AddIncludePath("-I$CMSSW_BASE/src/fjcontrib-1.009/include/fastjet/contrib");
   gSystem->AddIncludePath("-I/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/lib/slc5_amd64_gcc462/");
   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
   gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/lib/libfastjet.so");
   gSystem->Load("/afs/cern.ch/cms/slc5_amd64_gcc462/external/fastjet/3.0.1-cms3/lib/libfastjettools.so");
   gSystem->Load("libFWCoreFWLite.so");
   AutoLibraryLoader::enable();
   gSystem->Load("libDataFormatsFWLite.so");
    
   gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libCondFormatsJetMETObjects.so");



}
