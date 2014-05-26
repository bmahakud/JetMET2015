#define AnalysisAvJetMass2_cxx
#include "AnalysisAvJetMass2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisAvJetMass2::Loop()
{
int noEvents=1000;
double JetpTMin=200;
double JetpTMax=4000;
TH1F *JetmassAK8_h[41][8];//firstone is no of diff PV,second one is no of grooming 
TH1F *JetmassAK10_h[41][8];
TH1F *JetmassAK12_h[41][8];

   int nPVi=20;
   int nPVf=40;
  int NoGrooming=8;

  int NoOfhistos=nPVf-nPVi+1;
  static const int anPVi=20;
  static const int bnPVf=40; 
  static const int siz=bnPVf-anPVi+1;
  //cout<<"siz = "<<siz<<endl;
  char histname[100];
  double rangeMin=0;
  double rangeMax=1000;
  int NoOfbins=100;


  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK8%i_%i",jj,g0);
        JetmassAK8_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK8_h[jj][g0]->SetLineColor(3);
        JetmassAK8_h[jj][g0]->SetLineWidth(3);
        JetmassAK8_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK8_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK8_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK8_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK8_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK8_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }

  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK10%i_%i",jj,g0);
        JetmassAK10_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK10_h[jj][g0]->SetLineColor(3);
        JetmassAK10_h[jj][g0]->SetLineWidth(3);
        JetmassAK10_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK10_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK10_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK10_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK10_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK10_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }

  for(int g0=0;g0<NoGrooming;g0++){
for (int jj=nPVi;jj<nPVf+1;jj++) {
        sprintf(histname,"AveargejetmassAK12%i_%i",jj,g0);
        JetmassAK12_h[jj][g0] = new TH1F(histname,histname,NoOfbins,rangeMin,rangeMax);
        JetmassAK12_h[jj][g0]->SetLineColor(3);
        JetmassAK12_h[jj][g0]->SetLineWidth(3);
        JetmassAK12_h[jj][g0]->GetXaxis()->SetTitle("Number of primary Vertices");
        JetmassAK12_h[jj][g0]->GetXaxis()->SetTitleSize(0.055);
        JetmassAK12_h[jj][g0]->GetXaxis()->SetLabelSize(0.055);
        JetmassAK12_h[jj][g0]->GetYaxis()->SetTitle("Events");
        JetmassAK12_h[jj][g0]->GetYaxis()->SetTitleSize(0.055);
        JetmassAK12_h[jj][g0]->GetYaxis()->SetLabelSize(0.055);
      }
      }











   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
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
/////////////////////////////////////////////////////////////////////////////AK8
  int sizak8=GroomedJet_AK0_8_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak8;jj++) {//ak8loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd1_corr->at(jj));
                         }//sd1
         

       if(jj==0 && g0==1){//sd2
       JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd2_corr->at(jj));
                         }//sd2
  
     if(jj==0 && g0==2){//sd3
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd3_corr->at(jj));
                         }//sd3

     if(jj==0 && g0==3){//sd4
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_sd4_corr->at(jj));
                         }//sd4

     if(jj==0 && g0==4){//tr
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_ft_corr->at(jj));
                         }//tr

      if(jj==0 && g0==5){//tr
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_tr_corr->at(jj));
                         }//tr

      if(jj==0 && g0==6){//pr
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_pr_corr->at(jj));
                         }//pr

     if(jj==0 && g0==7){//ungroomed
      JetmassAK8_h[nPV][g0]->Fill(GroomedJet_AK0_8_mass_corr->at(jj));
                         }//ungroomed


          
       }//nPV criteria
}//ak8 loop
}//no of grooming loop
///////////////////////////////////////////////////////////////////////////AK8


///////////////////////////////////////////////////////////////////////////AK10

int sizak10=GroomedJet_AK1_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak10;jj++) {//ak10loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd1_corr->at(jj));
                         }//sd1
         

       if(jj==0 && g0==1){//sd2
       JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd2_corr->at(jj));
                         }//sd2
  
     if(jj==0 && g0==2){//sd3
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd3_corr->at(jj));
                         }//sd3

     if(jj==0 && g0==3){//sd4
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_sd4_corr->at(jj));
                         }//sd4

     if(jj==0 && g0==4){//tr
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_ft_corr->at(jj));
                         }//tr

      if(jj==0 && g0==5){//tr
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_tr_corr->at(jj));
                         }//tr

      if(jj==0 && g0==6){//pr
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_pr_corr->at(jj));
                         }//pr

     if(jj==0 && g0==7){//ungroomed
      JetmassAK10_h[nPV][g0]->Fill(GroomedJet_AK1_mass_corr->at(jj));
                         }//ungroomed



       
       }//nPV criteria
}//ak10 loop
}//no of grooming loop

///////////////////////////////////////////////////////////////////////////AK10


///////////////////////////////////////////////////////////////////////////////AK12

int sizak12=GroomedJet_AK1_2_pt_sd1_corr->size();
      for(int g0=0;g0<NoGrooming;g0++){//no of grooming loop
      for (int jj=0;jj<sizak12;jj++) {//ak12loop
      
      if(nPV >= nPVi && nPV <= nPVf ){//nPV criteria
      if(jj==0 && g0==0){//sd1
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd1_corr->at(jj));
                         }//sd1
         

       if(jj==0 && g0==1){//sd2
       JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd2_corr->at(jj));
                         }//sd2
  
     if(jj==0 && g0==2){//sd3
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd3_corr->at(jj));
                         }//sd3

     if(jj==0 && g0==3){//sd4
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_sd4_corr->at(jj));
                         }//sd4

     if(jj==0 && g0==4){//tr
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_ft_corr->at(jj));
                         }//tr

      if(jj==0 && g0==5){//tr
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_tr_corr->at(jj));
                         }//tr

      if(jj==0 && g0==6){//pr
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_pr_corr->at(jj));
                         }//pr

     if(jj==0 && g0==7){//ungroomed
      JetmassAK12_h[nPV][g0]->Fill(GroomedJet_AK1_2_mass_corr->at(jj));
                         }//ungroomed








           
       }//nPV criteria
}//ak12loop
}//no of grooming loop


///////////////////////////////////////////////////////////////////////////////ak12

   }//event loop

TCanvas *c1=new TCanvas("mycanvas1","My Canvas1");
gStyle->SetOptFit(111);
TMultiGraph *mg = new TMultiGraph();
mg->SetTitle(" Pile up dependence of AK8 jet mass (1100 Gev < Jet pT < 1400 GeV); Number of Primary Vertices ;  Average Jet mass in GeV");

double x_corr_mc_nPV[siz];
double y_corr_mc_ak8_sd1[siz];
double y_corr_mc_ak8_sd1_err[siz];

double y_corr_mc_ak8_sd2[siz];
double y_corr_mc_ak8_sd2_err[siz];

double y_corr_mc_ak8_sd3[siz];
double y_corr_mc_ak8_sd3_err[siz];

double y_corr_mc_ak8_sd4[siz];
double y_corr_mc_ak8_sd4_err[siz];

double y_corr_mc_ak8_ft[siz];
double y_corr_mc_ak8_ft_err[siz];

double y_corr_mc_ak8_tr[siz];
double y_corr_mc_ak8_tr_err[siz];

double y_corr_mc_ak8_pr[siz];
double y_corr_mc_ak8_pr_err[siz];

double y_corr_mc_ak8_ung[siz];
double y_corr_mc_ak8_ung_err[siz];
//////////////////////////////////ak10start
double y_corr_mc_ak10_sd1[siz];
double y_corr_mc_ak10_sd1_err[siz];

double y_corr_mc_ak10_sd2[siz];
double y_corr_mc_ak10_sd2_err[siz];

double y_corr_mc_ak10_sd3[siz];
double y_corr_mc_ak10_sd3_err[siz];

double y_corr_mc_ak10_sd4[siz];
double y_corr_mc_ak10_sd4_err[siz];

double y_corr_mc_ak10_ft[siz];
double y_corr_mc_ak10_ft_err[siz];

double y_corr_mc_ak10_tr[siz];
double y_corr_mc_ak10_tr_err[siz];

double y_corr_mc_ak10_pr[siz];
double y_corr_mc_ak10_pr_err[siz];

double y_corr_mc_ak10_ung[siz];
double y_corr_mc_ak10_ung_err[siz];
////////////////////////////////ak10end

////////////////////////////////ak12
double y_corr_mc_ak12_sd1[siz];
double y_corr_mc_ak12_sd1_err[siz];

double y_corr_mc_ak12_sd2[siz];
double y_corr_mc_ak12_sd2_err[siz];

double y_corr_mc_ak12_sd3[siz];
double y_corr_mc_ak12_sd3_err[siz];

double y_corr_mc_ak12_sd4[siz];
double y_corr_mc_ak12_sd4_err[siz];

double y_corr_mc_ak12_ft[siz];
double y_corr_mc_ak12_ft_err[siz];

double y_corr_mc_ak12_tr[siz];
double y_corr_mc_ak12_tr_err[siz];

double y_corr_mc_ak12_pr[siz];
double y_corr_mc_ak12_pr_err[siz];

double y_corr_mc_ak12_ung[siz];
double y_corr_mc_ak12_ung_err[siz];

////////////////////////////////ak12

for(int i0=0;i0<nPVf-nPVi+1;i0++){
x_corr_mc_nPV[i0]=i0+nPVi;
y_corr_mc_ak8_sd1[i0]=JetmassAK8_h[i0+nPVi][0]->GetMean();
y_corr_mc_ak8_sd1_err[i0]=JetmassAK8_h[i0+nPVi][0]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][0]->GetEntries());

y_corr_mc_ak8_sd2[i0]=JetmassAK8_h[i0+nPVi][1]->GetMean();
y_corr_mc_ak8_sd2_err[i0]=JetmassAK8_h[i0+nPVi][1]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][1]->GetEntries());

y_corr_mc_ak8_sd3[i0]=JetmassAK8_h[i0+nPVi][2]->GetMean();
y_corr_mc_ak8_sd3_err[i0]=JetmassAK8_h[i0+nPVi][2]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][2]->GetEntries());

y_corr_mc_ak8_sd4[i0]=JetmassAK8_h[i0+nPVi][3]->GetMean();
y_corr_mc_ak8_sd4_err[i0]=JetmassAK8_h[i0+nPVi][3]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][3]->GetEntries());

y_corr_mc_ak8_ft[i0]=JetmassAK8_h[i0+nPVi][4]->GetMean();
y_corr_mc_ak8_ft_err[i0]=JetmassAK8_h[i0+nPVi][4]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][4]->GetEntries());


y_corr_mc_ak8_tr[i0]=JetmassAK8_h[i0+nPVi][5]->GetMean();
y_corr_mc_ak8_tr_err[i0]=JetmassAK8_h[i0+nPVi][5]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][5]->GetEntries());

y_corr_mc_ak8_pr[i0]=JetmassAK8_h[i0+nPVi][6]->GetMean();
y_corr_mc_ak8_pr_err[i0]=JetmassAK8_h[i0+nPVi][6]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][6]->GetEntries());

y_corr_mc_ak8_ung[i0]=JetmassAK8_h[i0+nPVi][7]->GetMean();
y_corr_mc_ak8_ung_err[i0]=JetmassAK8_h[i0+nPVi][7]->GetRMS()/sqrt(JetmassAK8_h[i0+nPVi][7]->GetEntries());

}


TGraphErrors *graph1  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd1,0,y_corr_mc_ak8_sd1_err);
        graph1->SetLineColor(3);
        graph1->SetLineWidth(3);
        graph1->GetXaxis()->SetTitle("Number of primary Vertices");
        graph1->GetXaxis()->SetTitleSize(0.055);
        graph1->GetXaxis()->SetLabelSize(0.055);
        graph1->GetYaxis()->SetTitle("Events");
        graph1->GetYaxis()->SetTitleSize(0.055);
        graph1->GetYaxis()->SetLabelSize(0.055);
TGraphErrors *graph2  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd2,0,y_corr_mc_ak8_sd2_err);
graph2->SetLineColor(4);

TGraphErrors *graph3  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd3,0,y_corr_mc_ak8_sd3_err);
graph3->SetLineColor(5);

TGraphErrors *graph4  = new TGraphErrors(nPVf-nPVi+1,x_corr_mc_nPV,y_corr_mc_ak8_sd4,0,y_corr_mc_ak8_sd4_err);
graph4->SetLineColor(6);
//graphs.back()->GetXaxis()->SetRangeUser(0, 360.0);






















mg->Add(graph1,"lp");//lp
mg->Add(graph2,"lp");//lp
mg->Add(graph3,"lp");//lp
mg->Add(graph4,"lp");//lp


c1->cd();
mg->Draw("AP");
leg = new TLegend(0.6,0.7,0.89,0.89);

leg->AddEntry(graph1,"AK8 softdrop #beta = 0","lp");
leg->AddEntry(graph2,"AK8 softdrop #beta = 2","lp");
leg->AddEntry(graph3,"AK8 softdrop #beta = -1","lp");
leg->AddEntry(graph4,"AK8 softdrop #beta = 1","lp");


leg->Draw();






}//main programme
