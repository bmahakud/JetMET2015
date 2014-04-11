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
#include "fjcontrib-1.009/include/fastjet/contrib/SoftDrop.hh"
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

void readBacon(char *bs=NULL,std::string start="0",std::string end="0"){//main programme
       /// char *bs=NULL; 


/////Enter here the relevent informations about jets and algos////////
//////////User Entry Box start////////////////
//***Enter ConeRadius BELOW
        double R[15] = {0};//define what cone radius jets you want to use
        R[0]=1.2;
      //  R[1]=1.2;
      //  R[2]=1.2;
        //R[3]=0.45;


        //***Enter total No of variables BELOW
        int NoOfVariables=53;//Enter here the total no varibles you want to store ,like pt ,eta,phi.


        //below here enter the GlobalTag
        std::string GlobalTag = "START53_V7G";

        //****write the variable types BELOW if you want to add a new variable
        string VarType[NoOfVariables];
        VarType[0]="pt_uncorr";
        VarType[1]="number_jet_central";
        VarType[2]="mass_uncorr";
        VarType[3]="mass_tr_uncorr";
        VarType[4]="pt_tr_uncorr";
        VarType[5]="mass_ft_uncorr";
        VarType[6]="pt_ft_uncorr";
        VarType[7]="mass_pr_uncorr";
        VarType[8]="pt_pr_uncorr";
        VarType[9]="jet_area";
        VarType[10]="jet_eta_uncorr";
        VarType[11]="jet_phi_uncorr";
        VarType[12]="pt_corr";
        VarType[13]="mass_corr";
        VarType[14]="eta_corr";
        VarType[15]="phi_corr";
        VarType[16]="e_corr";
        VarType[17]="mass_tr_corr";
        VarType[18]="pt_tr_corr";
        VarType[19]="eta_tr_corr";
        VarType[20]="phi_tr_corr";
        VarType[21]="e_tr_corr";
        VarType[22]="area_tr_corr";

        VarType[23]="mass_ft_corr";
        VarType[24]="pt_ft_corr";
        VarType[25]="eta_ft_corr";
        VarType[26]="phi_ft_corr";
        VarType[27]="e_ft_corr";
        VarType[28]="area_ft_corr";

        VarType[29]="mass_pr_corr";
        VarType[30]="pt_pr_corr";
        VarType[31]="eta_pr_corr";
        VarType[32]="phi_pr_corr";
        VarType[33]="e_pr_corr";
        VarType[34]="area_pr_corr";

        VarType[35]="prsubjet1_px";
        VarType[36]="prsubjet1_py";
        VarType[37]="prsubjet1_pz";
        VarType[38]="prsubjet1_e";
        VarType[39]="prsubjet2_px";
        VarType[40]="prsubjet2_py";
        VarType[41]="prsubjet2_pz";
        VarType[42]="prsubjet2_e";

        VarType[43]="massdrop_pr_uncorr";
        VarType[44]="massdrop_pr";
        VarType[45]="tau1";
        VarType[46]="tau2";
        VarType[47]="tau3";
        VarType[48]="tau4";
        VarType[49]="tau2tau1";

        VarType[50]="qjetmass";
        VarType[51]="qjet_massdrop";
        VarType[52]="n_PV"; 
        //****Enter the type of jet algorithm BELOW
        string JetAlgorithm = "AK";///write AK for antikt,KT for kt ,CA for cambridgeAchen ,one algo is possible at a time for now
        //****Enter the following contants BELOW
        int activeAreaRepeats = 1;
        double ghostArea = 0.01;
        double ghostEtaMax = 4.4;
        static const int NUM_JET_MAX = 6;
        //////////User Entry Box end /////////////////
        /////////////////////////////////////////////////////////////////////


        ///finds the no of different R jets to be reclustered
        int NoOfdiffRJets=0;
        for(int j=0;j<15;j++){//2
          if(R[j] > 0.2){//1
          NoOfdiffRJets=NoOfdiffRJets+1;
          }//1
        }//2

        string Rs[NoOfdiffRJets];//a string array for storing the R values
        ///finds the no of different R jets to be reclustered

        ///converts Jet R values to string
        std::ostringstream os[NoOfdiffRJets];
        for(int ii=0;ii<NoOfdiffRJets;ii++){//3
        os[ii] << R[ii];
        }//3

        for(int i=0;i<NoOfdiffRJets;i++){//4//stores the string values of R in a string array
        std::string str = os[i].str();
        Rs[i]=str;
        cout<<"double to string= "<<Rs[i]<<endl;
        }//4


        string BranchName[NoOfVariables];
        ofstream write;
        write.open("BranchNames.txt");
        for(int jj=0;jj<NoOfdiffRJets;jj++){//5
        //if you want to create a new branch create a new string as its name

          for(int n=0;n<NoOfVariables;n++){//6
          BranchName[n]="GroomedJet_"+JetAlgorithm+Rs[jj]+"_"+VarType[n];
          write<<BranchName[n]<<endl;
          } //6

        }//5

        ///calculate total no of branches to becreated
        int NoOfBranches = NoOfdiffRJets*NoOfVariables;


        ///reads the names of branches from the file and makes branches
        ifstream infile("BranchNames.txt");
        
        const char *startchar;
        const char *endchar; 
        startchar = start.c_str();
        endchar   =end.c_str();
        char outputfileName[100];
 
        sprintf(outputfileName, "QCDptEvent%sto%s.root", startchar, endchar);




        TFile *f = new TFile(outputfileName,"RECREATE");
        TTree *tr = new TTree("GroomedJets","bacon groomed jets");
        int n_PVertex;
        tr->Branch("n_PVertex",&n_PVertex,"n_PVertex/I");
        std::vector<double> b[NoOfBranches]; //create no of vectors equals to NoOfBranches
        TBranch *B[5*NoOfBranches];// = new TBranch("branch1",&b);
        char s[NoOfBranches];
        for(int i=0; i<NoOfBranches; i++){ //7, branch create loop
        infile>>s;
        cout<<"creating a branch named ......... "<<s<<endl;
        B[i] = tr->Branch(s,&b[i]);

        }//7,branch create loop
        ///reads the names of branches from the file and makes branches


///////STARTS READING THE EXPERT NTUPLE////////////
//////////////////////////////////////////////////
	 
	char buff[50]; // buff is large enough to hold the entire formatted string
	sprintf(buff, "%s", bs);
       int startNumber = atoi(start.c_str()); 
       int endNumber = atoi(end.c_str()); 



        TFile* fIn = new TFile(buff);
        TTree* tree = (TTree*) fIn->Get("Events");
        TClonesArray *fPFPart = new TClonesArray("baconhep::TPFPart");
        TClonesArray *fGenPart = new TClonesArray("baconhep::TGenParticle");
        TClonesArray *fVertex = new TClonesArray("baconhep::TVertex");
     
        tree->SetBranchAddress("PFPart", &fPFPart);
        tree->SetBranchAddress("GenParticle", &fGenPart);
        tree->SetBranchAddress("PV",&fVertex);
        /////////////////////////////////////////////////////////////////
        cout <<"total no of events = "<<tree->GetEntriesFast()<<endl;
        for(int i0 = startNumber; i0 < endNumber/*tree->GetEntriesFast()*/; i0++) { //8, event loop
          //     if(i0==start){
          // break;
          // }
        if(i0%10==0)cout<<"Evt = "<<i0<<endl;
          for(int ij=0;ij<NoOfBranches;ij++){//clear the vectors before each event loop
                                b[ij].clear();
                                             }//clear the vectors before each event loop
        tree->GetEntry(i0);
        n_PVertex = fVertex->GetEntries();
        //From this PFpartcile stuff begins
        //////////////////////////////////////////////////////////////////////////////////////////////&&&&&&&&
        std::vector<fastjet::PseudoJet> particles; particles.clear();
        for( int i1 = 0; i1 < fPFPart->GetEntriesFast(); i1++){//9,entries loop,fill the vector particles with PF particles
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
  

        }//9,entries loop ,fill the vector particles with PFparticles
        ////////////////////////////////////////////////////////////////////////////////////////////////////&&&&&&&&&

        ////////////////////////////////Read Jec txt files

        std::vector< JetCorrectorParameters > jecPars;
        std::vector< std::string > jecStr;
        jecStr.push_back(GlobalTag+"_L1FastJet_AK7PF.txt");
        jecStr.push_back(GlobalTag+"_L2Relative_AK7PF.txt");
        jecStr.push_back(GlobalTag+"_L3Absolute_AK7PF.txt");

        JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(GlobalTag+"_L1FastJet_AK7PF.txt");
        JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(GlobalTag+"_L2Relative_AK7PF.txt");
        JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(GlobalTag+"_L3Absolute_AK7PF.txt");

        jecPars.push_back(*L1JetPar);
        jecPars.push_back(*L2JetPar);
        jecPars.push_back(*L3JetPar);
        FactorizedJetCorrector* jec_ = new FactorizedJetCorrector(jecPars);
        JetCorrectionUncertainty* jecUnc_=new JetCorrectionUncertainty( GlobalTag+"_Uncertainty_AK7PF.txt" );

        ////////////////////////////////Read Jec txt files
      
        //calculate rhoValue_grid
        double rhoEtaMax=4.4;
        fastjet::GridMedianBackgroundEstimator* mBgeGrid;
        fastjet::Subtractor* subtractor;
        subtractor=NULL;
        mBgeGrid= new fastjet::GridMedianBackgroundEstimator(rhoEtaMax, 0.55);
        mBgeGrid->set_particles(particles);
        if(subtractor) delete subtractor;
        subtractor= new fastjet::Subtractor(mBgeGrid);
        double rhoVal_grid=mBgeGrid->rho();
       // h1->Fill(rhoVal_grid);



 
        fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
        fjActiveArea.set_fj2_placement(true);
        fastjet::AreaDefinition *mAreaDefinition;
        mAreaDefinition =new fastjet::AreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );
        fastjet::Selector selected_eta = fastjet::SelectorAbsEtaMax(2.4);
        ///////reclustering starts from here////////////////////////////////////////////

        for(int jjj=0;jjj<NoOfdiffRJets;jjj++){//9,loop over the diff cone radius
        JetDefinition jet_def(antikt_algorithm,R[jjj]);


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

        double beta_sd = 1.0;
        double zcut_sd = 0.1;
        double mu_sd = 1.0;
        contrib::SoftDropTagger soft_drop_mmdt(0.0, zcut_sd, mu_sd);
        contrib::SoftDropTagger soft_drop_sdb2(2.0, zcut_sd, mu_sd);
        contrib::SoftDropTagger soft_drop_sdm1(-1.0, zcut_sd, mu_sd);
        
        
        
        std::vector<fastjet::Transformer const *> transformers;
        transformers.push_back(&trimmer);
        transformers.push_back(&filter);
        transformers.push_back(&pruner);








        int nconstituents0;
        int number_jet_central = out_jets.size();
        b[1+(jjj*NoOfVariables)].push_back(number_jet_central);


        for(int j = 0; j < number_jet_central && j<NUM_JET_MAX; j++) {//loop over the jets

          if(j==0){//11
            if (out_jets_basic.at(j).constituents().size() >= 100) nconstituents0 = 100;
else nconstituents0 = (int) out_jets_basic.at(j).constituents().size();
std::vector<fastjet::PseudoJet> cur_constituents = sorted_by_pt(out_jets_basic.at(j).constituents());
for (int aa = 0; aa < nconstituents0; aa++){ //10
// constituents0_eta[aa] = cur_constituents.at(aa).eta();
// constituents0_phi[aa] = cur_constituents.at(aa).phi();
// constituents0_e[aa] = cur_constituents.at(aa).e();
} //10

                  }//11

        b[0+(jjj*NoOfVariables)].push_back(out_jets.at(j).pt());
        b[2+(jjj*NoOfVariables)].push_back(out_jets.at(j).m());
        b[9+(jjj*NoOfVariables)].push_back(out_jets.at(j).area());
        b[10+(jjj*NoOfVariables)].push_back(out_jets.at(j).eta());
        b[11+(jjj*NoOfVariables)].push_back(out_jets.at(j).phi());

        jec_->setJetEta( out_jets.at(j).eta() );
jec_->setJetPt ( out_jets.at(j).pt() );
jec_->setJetE ( out_jets.at(j).e() );
jec_->setJetA ( out_jets.at(j).area() );
jec_->setRho ( rhoVal_grid );
jec_->setNPV ( fVertex->GetEntriesFast() );
b[52+(jjj*NoOfVariables)].push_back(fVertex->GetEntries());
        double corr = jec_->getCorrection();
        TLorentzVector jet_corr(corr*(out_jets.at(j).px()),corr*(out_jets.at(j).py()),corr*(out_jets.at(j).pz()),corr*(out_jets.at(j).e()));

        b[12+(jjj*NoOfVariables)].push_back(jet_corr.Pt());
        b[13+(jjj*NoOfVariables)].push_back(jet_corr.M());
        b[14+(jjj*NoOfVariables)].push_back(jet_corr.Eta());
        b[15+(jjj*NoOfVariables)].push_back(jet_corr.Phi());
        b[16+(jjj*NoOfVariables)].push_back(jet_corr.Energy());

        //cout <<"size trans "<<transformers.size()<<endl;
        //cout <<"End trans "<<transformers.end()<<endl;
        int transctr = 0;
        for ( std::vector<fastjet::Transformer const *>::const_iterator
itransf = transformers.begin(), itransfEnd = transformers.end();
itransf != itransfEnd; ++itransf ) {//transformed jet

           fastjet::PseudoJet transformedJet = out_jets.at(j);
transformedJet = (**itransf)(transformedJet);

  fastjet::PseudoJet transformedJet_basic = out_jets_basic.at(j);
transformedJet_basic = (**itransf)(transformedJet_basic);
                 if (transctr == 0){//12 Trimmer
                    b[3+(jjj*NoOfVariables)].push_back(transformedJet.m());
                    b[4+(jjj*NoOfVariables)].push_back(transformedJet.pt());
                    TLorentzVector jet_tr_corr(corr*(transformedJet.px()),corr*(transformedJet.py()),corr*(transformedJet.pz()),corr*(transformedJet.e()));

                    b[17+(jjj*NoOfVariables)].push_back(jet_tr_corr.M());
                    b[18+(jjj*NoOfVariables)].push_back(jet_tr_corr.Pt());
                    b[19+(jjj*NoOfVariables)].push_back(jet_tr_corr.Eta());
                    b[20+(jjj*NoOfVariables)].push_back(jet_tr_corr.Phi());
                    b[21+(jjj*NoOfVariables)].push_back(jet_tr_corr.Energy());
                    b[22+(jjj*NoOfVariables)].push_back(transformedJet.area());


                    }//12

                 else if (transctr ==1 ){//Filter
                   b[5+(jjj*NoOfVariables)].push_back(transformedJet.m());
                   b[6+(jjj*NoOfVariables)].push_back(transformedJet.pt());
                   TLorentzVector jet_ft_corr(corr*(transformedJet.px()),corr*(transformedJet.py()),corr*(transformedJet.pz()),corr*(transformedJet.e()));
                                           
                   b[23+(jjj*NoOfVariables)].push_back(jet_ft_corr.M());
                   b[24+(jjj*NoOfVariables)].push_back(jet_ft_corr.Pt());
                   b[25+(jjj*NoOfVariables)].push_back(jet_ft_corr.Eta());
                   b[26+(jjj*NoOfVariables)].push_back(jet_ft_corr.Phi());
                   b[27+(jjj*NoOfVariables)].push_back(jet_ft_corr.Energy());
                   b[28+(jjj*NoOfVariables)].push_back(transformedJet.area());
                     }

                 else if(transctr ==2){//prunner
                   b[7+(jjj*NoOfVariables)].push_back(transformedJet.m());
                   b[8+(jjj*NoOfVariables)].push_back(transformedJet.pt());
                      
                   TLorentzVector jet_pr_corr(corr*(transformedJet.px()),corr*(transformedJet.py()),corr*(transformedJet.pz()),corr*(transformedJet.e()));
                   b[29+(jjj*NoOfVariables)].push_back(jet_pr_corr.M());
                   b[30+(jjj*NoOfVariables)].push_back(jet_pr_corr.Pt());
                   b[31+(jjj*NoOfVariables)].push_back(jet_pr_corr.Eta());
                   b[32+(jjj*NoOfVariables)].push_back(jet_pr_corr.Phi());
                   b[33+(jjj*NoOfVariables)].push_back(jet_pr_corr.Energy());
                   b[34+(jjj*NoOfVariables)].push_back(transformedJet.area());


//decompose into requested no of subjets
                   if (transformedJet_basic.constituents().size() > 1){//15 decompose

                   int nsubjetstokeep = 2;

                   std::vector<fastjet::PseudoJet> subjets = transformedJet_basic.associated_cluster_sequence()->exclusive_subjets(transformedJet_basic,nsubjetstokeep);

                   TLorentzVector sj1( subjets.at(0).px(),subjets.at(0).py(),subjets.at(0).pz(),subjets.at(0).e());
                   TLorentzVector sj2( subjets.at(1).px(),subjets.at(1).py(),subjets.at(1).pz(),subjets.at(1).e());

                   b[35+(jjj*NoOfVariables)].push_back(subjets.at(0).px());
                   b[36+(jjj*NoOfVariables)].push_back(subjets.at(0).py());
                   b[37+(jjj*NoOfVariables)].push_back(subjets.at(0).pz());
                   b[38+(jjj*NoOfVariables)].push_back(subjets.at(0).e());
                   b[39+(jjj*NoOfVariables)].push_back(subjets.at(1).px());
                   b[40+(jjj*NoOfVariables)].push_back(subjets.at(1).py());
                   b[41+(jjj*NoOfVariables)].push_back(subjets.at(1).pz());
                   b[42+(jjj*NoOfVariables)].push_back(subjets.at(1).e());

                   TLorentzVector fullj = sj1 + sj2;


           if (subjets.at(0).m() >= subjets.at(1).m()){//16
            b[43+(jjj*NoOfVariables)].push_back(subjets.at(0).m()/transformedJet.m());
            b[44+(jjj*NoOfVariables)].push_back(subjets.at(0).m()/jet_pr_corr.M());

           } //16
           else {//17
            b[43+(jjj*NoOfVariables)].push_back(subjets.at(0).m()/transformedJet.m());
            b[44+(jjj*NoOfVariables)].push_back(subjets.at(0).m()/jet_pr_corr.M());
           }//17

        }//15 decompose


        }
        else{ std::cout << "error in number of transformers" << std::endl;}
 
        transctr++;
        }//transformed jet

        double mNsubjettinessKappa=1;
        float tau1,tau2,tau3,tau4,tau2tau1;

        double beta = mNsubjettinessKappa; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadrat ic/classic k-means
        double R0 = R[jjj]; // Characteristic jet radius for normalization
        double Rcut = R[jjj]; // maximum R particles can be from axis to be included in jet

        fastjet::Nsubjettiness nSub1KT(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
fastjet::Nsubjettiness nSub2KT(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
fastjet::Nsubjettiness nSub3KT(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
fastjet::Nsubjettiness nSub4KT(4, Njettiness::onepass_kt_axes, beta, R0, Rcut);

        tau1 = nSub1KT(out_jets.at(j));
        tau2 = nSub2KT(out_jets.at(j));
        tau3 = nSub3KT(out_jets.at(j));
        tau4 = nSub4KT(out_jets.at(j));
        tau2tau1 = tau2/tau1;
        b[45+(jjj*NoOfVariables)].push_back(tau1);
        b[46+(jjj*NoOfVariables)].push_back(tau2);
        b[47+(jjj*NoOfVariables)].push_back(tau3);
        b[48+(jjj*NoOfVariables)].push_back(tau4);
        b[49+(jjj*NoOfVariables)].push_back(tau2tau1);
  


        // cores computation -------------
// Begining the core computation
std::vector<fastjet::PseudoJet> constits = thisClustering_area.constituents(out_jets.at(j));
for (int kk = 0; kk < 11; ++kk){
double coreCtr = (double) kk;
if (coreCtr < R[jjj]*10.){
float m_core = 0, pt_core = 0;


        fastjet::JetDefinition jetDef_rcore(fastjet::cambridge_algorithm, coreCtr/10.);

fastjet::ClusterSequence thisClustering(constits, jetDef_rcore);
        std::vector<fastjet::PseudoJet> out_jets_core = sorted_by_pt(thisClustering.inclusive_jets(0.0));
m_core = out_jets_core.at(0).m();
pt_core = out_jets_core.at(0).pt();	

// computeCore( constits, coreCtr/10., tmpm, tmppt );
// if (m_core > 0) rcores[kk][j] = m_core/out_jets_core.at(j).m();
// if (pt_core > 0) ptcores[kk][j] = pt_core/out_jets_core.at(j).pt();
}
}

//std::cout<< "Ending the core computation" << endl;



        //Begining the planarflow computation


        //Ending the planarflow computation


        //begin qjets computation -------------

        if(j==0){//do qjets only for the hardest jet in the event!
          double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1);

          QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity);
          fastjet::JetDefinition qjet_def(&qjet_plugin);
          vector<fastjet::PseudoJet> constitsq;
          unsigned int nqjetconstits = out_jets_basic.at(j).constituents().size();
          int mQJetsPreclustering=30;

         if (nqjetconstits < (unsigned int) mQJetsPreclustering) constitsq = out_jets_basic.at(j).constituents();
else constitsq = out_jets_basic.at(j).associated_cluster_sequence()->exclusive_subjets_up_to(out_jets_basic.at(j),mQJetsPreclustering);
        int mQJetsN=50;
        for(unsigned int ii = 0 ; ii < (unsigned int) mQJetsN ; ii++){//mQJetsN loop
        fastjet::ClusterSequence qjet_seq(constitsq, qjet_def);
vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq.inclusive_jets(20.0));
        if (inclusive_jets2.size()>0) {//20
b[50+(jjj*NoOfVariables)].push_back(inclusive_jets2[0].m());
if (inclusive_jets2[0].constituents().size() > 1){//21
vector<fastjet::PseudoJet> subjets_qjet = qjet_seq.exclusive_subjets(inclusive_jets2[0],2);
if (subjets_qjet.at(0).m() >= subjets_qjet.at(1).m()){//22
b[51+(jjj*NoOfVariables)].push_back( (subjets_qjet.at(0).m()/inclusive_jets2[0].m()));
}//22
else{//23
b[51+(jjj*NoOfVariables)].push_back( (subjets_qjet.at(1).m()/inclusive_jets2[0].m()));
}//23
}//21
else{//24
b[51+(jjj*NoOfVariables)].push_back(1.);
}//24
}//20


                       else{//19
b[51+(jjj*NoOfVariables)].push_back(1.);
}//19

        }//mQJetsN loop


        }//do qjets only for the hardest jet in the event!


//end qjets computation --------------



}//loop over the jets


tr->Fill();

}//9, loop over diff cone radius

}//8, event loop

//h1->Draw();
f->Write();


}//main programme
