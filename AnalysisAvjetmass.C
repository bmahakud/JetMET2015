#define AnalysisAvJetMass_cxx
#include "AnalysisAvJetMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>

void AnalysisAvJetMass::Loop()
{

//////////////////////////////////////////////////////////////////////////////////
//User options
int noEvents=1000;
string Jetalgo="AK";

double R[15] = {0};//define what cone radius jets you want to use for now only THREE are available,0.8,1.0,1.2
        R[0]=0.8;
        R[1]=1.0;
        R[2]=1.2;

//grooming techniques
string Gt[10];

Gt[0]="sd1";
Gt[1]="sd2";
Gt[2]="sd3";
Gt[3]="sd4";
Gt[4]="ft";
Gt[5]="tr";
Gt[6]="pr";
//User Options
////////////////////////////////////////////////////////////////////////////////////

          int NoOfdiffRJets=0;
          for(int j=0;j<15;j++){//2
          if(R[j] > 0.2){//1
          NoOfdiffRJets=NoOfdiffRJets+1;
          }//1
          }//2
         

          int NoGrooming=0;
          for(int j=0;j<10;j++){//2
          if(Gt[j] == "sd1" || Gt[j] == "sd2" || Gt[j] == "sd3" ||  Gt[j] == "sd4" || Gt[j] == "ft" || Gt[j] == "tr" || Gt[j] == "pr"){//1
          NoGrooming=NoGrooming+1;
          }//1
          }//2
          //cout<<"NoOfGrooming = "<<NoGrooming<<endl;
          static const int NumGrooming=NoGrooming; 

///////////////////////////////////////////


      char *a1="AK0_8";
      char *a2="pt";
      char *a3="sd1"; 
      char *a4="corr";
      char *a5="AK8_softdrop_sd1";
        
      char name[100]; 
      char label[500];
   
      sprintf(name, "GroomedJet_%s_%s_%s_%s->size() \n", a1,a2,a3,a4);
      sprintf(label, "%s \n",a5);
    


TFile *f = new TFile("Histo.root","RECREATE");
gStyle->SetOptStat(0000);
  double JetpTMin=200;
  double JetpTMax=4000;
   int nPVi=20;
   int nPVf=40;
  int noOfGroomingtech=0;
  int NoOfhistos=nPVf-nPVi+1;
  static const int anPVi=20;
  static const int bnPVf=40; 
  static const int siz=bnPVf-anPVi+1;
  //cout<<"siz = "<<siz<<endl;
  char histname[100];
  double rangeMin=0;
  double rangeMax=1000;
  int NoOfbins=100;
  TH1F *AvJetmass_h[bnPVf];

  ////////////////////////////////trial start1
  TH1F *Jetmass_h[41][7];
  
    

for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"Aveargejetmass%i_%i",jj,g0);
        Jetmass_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
       
      }
      }

/////////////////////////////////trial end1







  for (int j=nPVi;j<nPVf+1;j++) {
        sprintf(histname,"Aveargejetmass%i",j);
        AvJetmass_h[j] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        AvJetmass_h[j]->SetLineColor(3);
        AvJetmass_h[j]->SetLineWidth(3);
        AvJetmass_h[j]->GetXaxis()->SetTitle("Number of primary Vertices");
        AvJetmass_h[j]->GetXaxis()->SetTitleSize(0.055);
        AvJetmass_h[j]->GetXaxis()->SetLabelSize(0.055);
        AvJetmass_h[j]->GetYaxis()->SetTitle("Events");
        AvJetmass_h[j]->GetYaxis()->SetTitleSize(0.055);
        AvJetmass_h[j]->GetYaxis()->SetLabelSize(0.055);
      }






          if (fChain == 0) return;

          Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//event loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      // if (Cut(ientry) < 0) continue;
///////////////////////////////////////noOfevents to run
    if(jentry==noEvents){
       break;
     } 
//////////////////////////////////////noOfevents to run

 int nPV=n_PVertex;
    
   for (int jj=0;jj<4;jj++) {
      if(jj==0){
      if(nPV >= nPVi && nPV <= nPVf ){
      AvJetmass_h[nPV]->Fill(GroomedJet_AK0_8_mass_sd1_corr->at(jj));
         }
        }
           
    }

/////////////////////////////////////////////triAlstart2







  int sizj=GroomedJet_AK1_pt_sd1_corr->size();
for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
 for (int jj=0;jj<sizj;jj++) {
     
 if(jj==0 && g0==0){//sd1
      if(nPV >= nPVi && nPV <= nPVf ){
      Jetmass_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd1_corr->at(jj));
         }
        }//sd1


       






           
    }
}//no of grooming loop

//////////////////////////////////////trial end2











   }//event loop


TCanvas *c1=new TCanvas("mycanvas1","My Canvas1");
gStyle->SetOptFit(111);
TMultiGraph *mg = new TMultiGraph();
mg->SetTitle(" Pile up dependence of AK8 jet mass (1100 Gev < Jet pT < 1400 GeV); Number of Primary Vertices ;  Average Jet mass");
//int s=nPVf-nPVi+1;
double x_uncorr_mc[siz];
double y_uncorr_mc[siz];
double y_uncorr_mc_err[siz];


for(int i0=0;i0<nPVf-nPVi+1;i0++){
x_uncorr_mc[i0]=i0+nPVi;
y_uncorr_mc[i0]=AvJetmass_h[i0+nPVi]->GetMean();
y_uncorr_mc_err[i0]=AvJetmass_h[i0+nPVi]->GetRMS()/sqrt(AvJetmass_h[i0+nPVi]->GetEntries());
}



TGraphErrors *graph1  = new TGraphErrors(nPVf-nPVi+1,x_uncorr_mc,y_uncorr_mc,0,y_uncorr_mc_err);
mg->Add(graph1,"lp");//lp
c1->cd();
mg->Draw("AP");
leg = new TLegend(0.6,0.7,0.89,0.89);

leg->AddEntry(graph1,label,"lp");

leg->Draw();

f->Write();















}
