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
//#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "JetSubstructure/SubstructureTools/interface/PseudoJetUserInfo.h"
//#include "JetSubstructure/SubstructureTools/interface/JetSubstructureTools.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPFPart.hh"
#include "../interface/TGenEventInfo.hh"
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
/////Enter here the relevent informations about jets and algos////////
//////////User Entry Box start////////////////
//***Enter ConeRadius BELOW
double R[15] = {0};//define what cone radius jets you want to use
R[0]=0.3;
R[1]=0.35; 
R[2]=0.4;
R[3]=0.45;
//***Enter total No of variables  BELOW
int NoOfVariables=6;//Enter here the total no varibles you want to store ,like pt ,eta,phi.


//****write the variable types BELOW if you want to add a new variable
string VarType[NoOfVariables];
VarType[0]="pt";
VarType[1]="eta";
VarType[2]="phi";
VarType[3]="e";
VarType[4]="mass";
VarType[5]="area";

//****Enter the type of jet algorithm  BELOW
string JetAlgorithm = "AK";///write AK for antikt,KT for kt ,CA for cambridgeAchen ,one algo is possible at a time for now
//****Enter the following contants BELOW
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 4.4;
    static const int NUM_JET_MAX = 6;
//////////User Entry Box end /////////////////
/////////////////////////////////////////////////////////////////////
int NoOfdiffRJets=0;
for(int j=0;j<15;j++){
if(R[j] > 0.2){
NoOfdiffRJets=NoOfdiffRJets+1;
}
}
string Rs[NoOfdiffRJets];

std::ostringstream os[NoOfdiffRJets];
for(int ii=0;ii<NoOfdiffRJets;ii++){
os[ii] << R[ii];
}

for(int i=0;i<NoOfdiffRJets;i++){
std::string str = os[i].str();
Rs[i]=str;
cout<<"double to string=  "<<Rs[i]<<endl;
}

string BranchName[NoOfVariables];
ofstream write;
write.open("BranchNames.txt");
for(int jj=0;jj<NoOfdiffRJets;jj++){
//if you want to create a new branch create a new string as its name 

for(int n=0;n<NoOfVariables;n++){
BranchName[n]="GroomedJet_"+JetAlgorithm+Rs[jj]+"_"+VarType[n];
write<<BranchName[n]<<endl; 
}

}

int NoOfBranches = NoOfdiffRJets*NoOfVariables;

ifstream infile("BranchNames.txt");
TFile *f = new TFile("BaconGroomedJets.root","RECREATE");
TTree *tr = new TTree("GroomedJets","bacon groomed jets");
std::vector<double> b[NoOfBranches];   
TBranch *B[5*NoOfBranches];// = new TBranch("branch1",&b);
char s[NoOfBranches];
for(int i=0; i<NoOfBranches; i++){ //branch create loop
    infile>>s;
    cout<<"creating a branch named .........       "<<s<<endl;
    B[i] = tr->Branch(s,&b[i]);

}

//################

    TFile* fIn = new TFile("/afs/cern.ch/user/b/bmahakud/work/NewBaconFramework/CMSSW_5_3_13/src/BaconProd/Ntupler/python/ntuple.root");
    TTree* tree = (TTree*) fIn->Get("Events");
    TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
    TClonesArray *fGenPart = new TClonesArray("baconhep::TGenParticle");
    
    tree->SetBranchAddress("PFPart", &fPFPart);
    tree->SetBranchAddress("GenParticle", &fGenPart);
/////////////////////////////////////////////////////////////////

    for(int i0 = 0; i0 < tree->GetEntriesFast(); i0++) { //loop3 event loop

          for(int ij=0;ij<NoOfBranches;ij++){//clear the vectors
                                b[ij].clear();
                                             }//clear the vectors      
        tree->GetEntry(i0);

        
        // read in the event tree
    //From this PFpartcile stuff begins
        std::vector<fastjet::PseudoJet> particles;  particles.clear();       
        for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//entries loop
            baconhep::TPFPart *pPartTmp = (baconhep::TPFPart*)((*fPFPart)[i1]);
        //    std::cout << "fPFPart["<<"0"<<"].pt = " << pPartTmp->pt << std::endl;:
          //  if (i1 > 10) break;

        double Px=pPartTmp->pt*cos(pPartTmp->phi);
        double Py= pPartTmp->pt*sin(pPartTmp->phi);
        double theta = 2*atan(exp(-pPartTmp->eta));   //eta = -ln(tan(theta/2))
	double Pz = pPartTmp->pt/tan(theta);
        double E = pPartTmp->e;
        double pdgId = pPartTmp->pfType;
        int    charge = pPartTmp->q;
         fastjet::PseudoJet tmp_psjet(Px, Py, Pz, E);
         tmp_psjet.set_user_info( new PseudoJetUserInfo(pdgId, charge) ); 
 //tmp_psjet.set_user_info(new PseudoJetUserInfo(5,-1));
         particles.push_back(tmp_psjet);
  

        }//entries loop


    bool   mSaveConstituents;
 
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fjActiveArea.set_fj2_placement(true);
    fastjet::AreaDefinition *mAreaDefinition;
    mAreaDefinition =new fastjet::AreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );
    fastjet::Selector selected_eta = fastjet::SelectorAbsEtaMax(2.4);
///////reclustering starts from here//////////////////////////////////////////// 

for(int jjj=0;jjj<NoOfdiffRJets;jjj++){//loop over the diff cone radius  
JetDefinition jet_def(antikt_algorithm,R[jjj]);


//jet_def= new fastjet::JetDefinition(fastjet::antikt_algorithm, R[jjj]);


	if (JetAlgorithm == "AK")jet_def.set_jet_algorithm( fastjet::antikt_algorithm );
	else if (JetAlgorithm == "CA")jet_def.set_jet_algorithm( fastjet::cambridge_algorithm );
	else if (JetAlgorithm == "KT")jet_def.set_jet_algorithm( fastjet::kt_algorithm );
	//else throw  << " unknown jet algorithm " << std::endl;

//JetDefinition jet_def(antikt_algorithm, R[jjj]);
  // run the clustering, extract the jets
  fastjet::ClusterSequenceArea thisClustering_area(particles, jet_def, *mAreaDefinition);
  fastjet::ClusterSequence thisClustering_basic(particles, jet_def);

   std::vector<fastjet::PseudoJet> out_jets  = sorted_by_pt( selected_eta(thisClustering_area.inclusive_jets(15.0)) );
   std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt( selected_eta(thisClustering_basic.inclusive_jets(15.0)) );
   
fastjet::Filter trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm,0.2),fastjet::SelectorPtFractionMin(0.03)));

fastjet::Filter filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm,0.3),fastjet::SelectorNHardest(3)));
  
fastjet::Pruner pruner(fastjet::cambridge_algorithm, 0.1, 0.5);


std::vector<fastjet::Transformer const *> transformers;
transformers.push_back(&trimmer);
transformers.push_back(&filter);
transformers.push_back(&pruner);


int number_jet_central = out_jets.size();


for(int j = 0; j < number_jet_central && j<NUM_JET_MAX; j++) {//loop over the jets
b[5+(jjj*NoOfVariables)].push_back(out_jets.at(j).area());

//TLorentzVector jet_corr = getCorrectedJet(out_jets.at(j),0);


for ( std::vector<fastjet::Transformer const *>::const_iterator 
					itransf = transformers.begin(), itransfEnd = transformers.end(); 
					itransf != itransfEnd; ++itransf ) {




                        fastjet::PseudoJet transformedJet = out_jets.at(j);
		     	transformedJet = (**itransf)(transformedJet);

		    	fastjet::PseudoJet transformedJet_basic = out_jets_basic.at(j);
	          	transformedJet_basic = (**itransf)(transformedJet_basic);


}
}//loop over the jets



//cout << "Clustering with " << jet_def.description() << endl;
// print the jets
 // cout <<   "        pt y phi" << endl;
  for (unsigned i = 0; i < out_jets_basic.size(); i++) {//loop2

    
      b[0+(jjj*NoOfVariables)].push_back(out_jets_basic[i].perp());
      b[1+(jjj*NoOfVariables)].push_back(out_jets_basic[i].rap());
      b[2+(jjj*NoOfVariables)].push_back(out_jets_basic[i].phi());
      b[3+(jjj*NoOfVariables)].push_back(out_jets_basic[i].e());
      b[4+(jjj*NoOfVariables)].push_back(out_jets_basic[i].m());
//    vector<PseudoJet> constituents = out_jets_basic[i].constituents();
 //   for (unsigned j = 0; j < constituents.size(); j++) {//loop1

   // }//loop1

}//loop2



//std::cout<<"no of jets = "<<b[0].size()<<std::endl;
tr->Fill();

}//loop over diff cone radius 

}//loop3 event loop



f->Write();




}
