#include <vector>
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
//#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
#include "JetSubstructure/SubstructureTools/interface/Njettiness.hh"
#include "JetSubstructure/SubstructureTools/interface/Nsubjettiness.hh"
#include "JetSubstructure/SubstructureTools/src/QjetsPlugin.h"
#include "fastjet/tools/Subtractor.hh"

//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
//#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h"
#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#endif
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPFPart.hh"
#include "../interface/TGenEventInfo.hh"
//#include "../interface/TEventInfo.hh"
#include "../interface/TGenParticle.hh"
#include <TMath.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <string>
#include <fstream>
#include "TBranch.h"
#include "TLorentzVector.h"


using namespace std;
using namespace fastjet;

using namespace baconhep;

void readBacon(){


        TFile* fIn = new TFile("ntuple.root");
        TTree* tree = (TTree*) fIn->Get("Events");
        TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
        TClonesArray *fGenPart = new TClonesArray("baconhep::TGenParticle");
        TClonesArray *fVertex = new TClonesArray("baconhep::TVertex");
     
        tree->SetBranchAddress("PFPart", &fPFPart);
        tree->SetBranchAddress("GenParticle", &fGenPart);
        tree->SetBranchAddress("PV",&fVertex);

for(int i0 = 0; i0 <tree->GetEntriesFast(); i0++) { //1, event loop



tree->GetEntry(i0);

std::vector<fastjet::PseudoJet> particles; particles.clear();
        for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//2,entries loop,fill the vector particles with PF particles
            baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);

        double Px=pPartTmp->pt*cos(pPartTmp->phi);
        double Py= pPartTmp->pt*sin(pPartTmp->phi);
        double theta = 2*atan(exp(-pPartTmp->eta)); //eta = -ln(tan(theta/2))
        double Pz = pPartTmp->pt/tan(theta);
        double E = pPartTmp->e;
        double pdgId = pPartTmp->pfType;
        int charge = pPartTmp->q;
        fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
        tmp_psjet.set_user_info( new PseudoJetUserInfo(pdgId, charge) );
        particles.push_back(tmp_psjet);
  

        }//2,entries loop ,fill the vector particles with PFparticles


        double rhoEtaMax=4.4;
        fastjet::GridMedianBackgroundEstimator* mBgeGrid;
        fastjet::Subtractor* subtractor;
        subtractor=NULL;
        mBgeGrid= new fastjet::GridMedianBackgroundEstimator(rhoEtaMax, 0.55);
        mBgeGrid->set_particles(particles);
        if(subtractor) delete subtractor;
        subtractor= new fastjet::Subtractor(mBgeGrid);
        double rhoVal_grid=mBgeGrid->rho();
         fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
        fjActiveArea.set_fj2_placement(true);
        fastjet::AreaDefinition *mAreaDefinition;
        mAreaDefinition =new fastjet::AreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );
        fastjet::Selector selected_eta = fastjet::SelectorAbsEtaMax(2.4);
        
        double R=0.8;//define the jet cone radius to use 
        std::string JetAlgorithm = "AK"; //choose which jet algo to use , write AK for antikt, CA for Cambridge/Achen, KT for kt algo
        
        JetDefinition jet_def(antikt_algorithm,R);


      if (JetAlgorithm == "AK")jet_def.set_jet_algorithm( fastjet::antikt_algorithm );
      else if (JetAlgorithm == "CA")jet_def.set_jet_algorithm( fastjet::cambridge_algorithm );
      else if (JetAlgorithm == "KT")jet_def.set_jet_algorithm( fastjet::kt_algorithm );
       //else throw << " unknown jet algorithm " << std::endl;

        fastjet::ClusterSequenceArea thisClustering_area(particles, jet_def, *mAreaDefinition);
        fastjet::ClusterSequence thisClustering_basic(particles, jet_def);



        std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt( selected_eta(thisClustering_area.inclusive_jets(15.0)) );
        std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt( selected_eta(thisClustering_basic.inclusive_jets(15.0)) );
  


        fastjet::Filter trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm,0.2),fastjet::SelectorPtFractionMin(0.03)));

        fastjet::Filter filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm,0.3),fastjet::SelectorNHardest(3)));
  
        fastjet::Pruner pruner(fastjet::cambridge_algorithm, 0.1, 0.5);


//Now you can print the jet momentum , mass , pt etc like the following
        for(int j=0;j<)out_jets.size();j++{
        cout<<"jet pt = "<<out_jets.at(j).pt()<<endl;
        }






















}//1 event loop










}


