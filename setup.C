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





   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconProdNtupler.so");
   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconProdUtils.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libCMGToolsExternal.so");   

gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libCommonToolsRecoUtils.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libCondFormatsEgammaObjects.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libCondFormatsJetMETObjects.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libDataFormatsJetReco.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libDataFormatsMETReco.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libEGammaEGammaAnalysisTools.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libJetMETCorrectionsType1MET.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libJetSubstructureSubstructureTools.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libJetToolsAnalyzerToolbox.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMMozerpowhegweight.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMuScleFitCalibration.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libPhysicsToolsKinFitter.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libPhysicsToolsPatAlgos.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libPhysicsToolsPatUtils.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libRecoJetsJetProducers.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libRecoMETMETPUSubtraction.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libSHarperHEEPAnalyzer.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginBaconProdNtuplerPlugins.so");

gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginCMGToolsExternalCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginCMGToolsExternalsplugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginCondFormatsEgammaObjectsCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginCondFormatsJetMETObjectsCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginDataFormatsJetRecoCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginDataFormatsMETRecoCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginJetMETCorrectionsType1MET_plugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginJetToolbox.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginJetToolsAnalyzerToolboxCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginMETFilters_plugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginMMozerpowhegweightPlugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginMuScleFitCalibrationPlugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginPFPUAssoMapPlugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginPhysicsToolsPatAlgos_plugins.so");


gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginPhysicsToolsPatAlgos_testAnalyzers.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginPhysicsToolsPatUtilsCapabilities.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginPhysicsToolsPatUtils_plugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoJetsJetProducers_plugins.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoMETMETAnalyzers.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoMETMETPUSubtraction_plugins.so");
 gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginRecoMETMETProducers.so");
gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/pluginSHarperHEEPAnalyzer_plugins.so");








}
