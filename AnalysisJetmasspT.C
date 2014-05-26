#define AnalysisJetmasspt_cxx
#include "AnalysisJetmasspt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalysisJetmasspt::Loop()
{

TFile *f = new TFile("JetmassplotsQCDpu40pt300to470.root","RECREATE");
gStyle->SetOptStat(0000);


double JetptcutFormass=100;

double JetptcutMin=200;
double JetptcutMax=4000;

TCanvas *c1=new TCanvas("My canvas1","Canvas1");
c1->Divide(2,2);

TCanvas *c2=new TCanvas("My canvas2","Canvas2");
c2->Divide(2,2);



////////////////////////nPV,Jet mass,Jet pt
TH1F *h1=new TH1F("h1","Event Primary Vertices",100,0,100);
h1->SetLineColor(2);
h1->SetLineWidth(3);
h1->GetXaxis()->SetTitle("Number of primary Vertices");
h1->GetXaxis()->SetTitleSize(0.055);
h1->GetXaxis()->SetLabelSize(0.055);
h1->GetYaxis()->SetTitle("Events");
h1->GetYaxis()->SetTitleSize(0.055);
h1->GetYaxis()->SetLabelSize(0.055);

TH1F *h2=new TH1F("h2","Leading Jet Pt distribution, AK8 jets",250,0,2000);
h2->SetLineWidth(3);
h2->SetLineColor(2);
h2->GetXaxis()->SetTitle("Jet pT in GeV");
h2->GetXaxis()->SetTitleSize(0.055);
h2->GetXaxis()->SetLabelSize(0.055);
h2->GetYaxis()->SetTitle("Events/8 GeV");
h2->GetYaxis()->SetTitleSize(0.055);
h2->GetYaxis()->SetLabelSize(0.055);





TH1F *h3=new TH1F("h3ft","Groomed jet mass in all n_PV for AK8 jets ",100,0,300);
h3->SetLineColor(5);
h3->SetLineWidth(3);
h3->GetXaxis()->SetTitle("Jet mass in GeV");
h3->GetXaxis()->SetTitleSize(0.055);
h3->GetXaxis()->SetLabelSize(0.055);
h3->GetYaxis()->SetTitle("Normalized to Unit");
h3->GetYaxis()->SetTitleSize(0.055);
h3->GetYaxis()->SetLabelSize(0.055);



TH1F *h4=new TH1F("h4tr","Trimmed jet mass",100,0,300);
h4->SetLineColor(3);
h4->SetLineWidth(3);
h4->GetXaxis()->SetTitle("Jet mass in GeV");
h4->GetXaxis()->SetTitleSize(0.055);
h4->GetXaxis()->SetLabelSize(0.055);
h4->GetYaxis()->SetTitle("Normalized to Unit");
h4->GetYaxis()->SetTitleSize(0.055);
h4->GetYaxis()->SetLabelSize(0.055);



TH1F *h5=new TH1F("h5pr","Prunned jet mass",100,0,300);
h5->SetLineColor(4);
h5->SetLineWidth(3);
h5->GetXaxis()->SetTitle("Jet mass in GeV");
h5->GetXaxis()->SetTitleSize(0.055);
h5->GetXaxis()->SetLabelSize(0.055);
h5->GetYaxis()->SetTitle("Normalized to Unit");
h5->GetYaxis()->SetTitleSize(0.055);
h5->GetYaxis()->SetLabelSize(0.055);



TH1F *h6=new TH1F("h6","Soft drop #beta=0",100,0,300);
h6->SetLineColor(2);
h6->SetLineWidth(3);
//h6->GetYaxis()->SetRangeUser(0,0.12);
h6->GetXaxis()->SetTitle("Jet mass in GeV");
h6->GetXaxis()->SetTitleSize(0.055);
h6->GetXaxis()->SetLabelSize(0.055);
h6->GetYaxis()->SetTitle("Normalized to Unit");
h6->GetYaxis()->SetTitleSize(0.055);
h6->GetYaxis()->SetLabelSize(0.055);



TH1F *h7=new TH1F("h7","Softdrop #beta=2",100,0,300);
h7->SetLineColor(6);
h7->SetLineWidth(3);
h7->GetXaxis()->SetTitle("Jet mass in GeV");
h7->GetXaxis()->SetTitleSize(0.055);
h7->GetXaxis()->SetLabelSize(0.055);
h7->GetYaxis()->SetTitle("Normalized to Unit");
h7->GetYaxis()->SetTitleSize(0.055);
h7->GetYaxis()->SetLabelSize(0.055);



TH1F *h8=new TH1F("h8","Softdrop #beta=-1",100,0,300);
h8->SetLineColor(8);
h8->SetLineWidth(3);
h8->GetXaxis()->SetTitle("Jet mass in GeV");
h8->GetXaxis()->SetTitleSize(0.055);
h8->GetXaxis()->SetLabelSize(0.055);
h8->GetYaxis()->SetTitle("Normalized to Unit");
h8->GetYaxis()->SetTitleSize(0.055);
h8->GetYaxis()->SetLabelSize(0.055);




TH1F *h9=new TH1F("h9","Softdrop #beta=1",100,0,300);
h9->SetLineColor(7);
h9->SetLineWidth(3);
h9->GetXaxis()->SetTitle("Jet mass in GeV");
h9->GetXaxis()->SetTitleSize(0.055);
h9->GetXaxis()->SetLabelSize(0.055);
h9->GetYaxis()->SetTitle("Normalized to Unit");
h9->GetYaxis()->SetTitleSize(0.055);
h9->GetYaxis()->SetLabelSize(0.055);




TH1F *h10=new TH1F("h10Ung","Ungroomed jet mass",100,0,300);
h10->SetLineColor(9);
h10->SetLineWidth(3);
h10->GetXaxis()->SetTitle("Jet mass in GeV");
h10->GetXaxis()->SetTitleSize(0.055);
h10->GetXaxis()->SetLabelSize(0.055);
h10->GetYaxis()->SetTitle("Events/ 3 GeV");
h10->GetYaxis()->SetTitleSize(0.055);
h10->GetYaxis()->SetLabelSize(0.055);

//////////////////////////////////average jet mass plot

double rangeAK10=450;
//////////////////////////////////////plots for AK10
TH1F *h2t=new TH1F("h2t","Leading Jet Pt distribution, AK10 jets",250,0,2000);
h2t->SetLineWidth(3);
h2t->SetLineColor(2);
h2t->GetXaxis()->SetTitle("Jet pT in GeV");
h2t->GetXaxis()->SetTitleSize(0.055);
h2t->GetXaxis()->SetLabelSize(0.055);
h2t->GetYaxis()->SetTitle("Events/8 GeV");
h2t->GetYaxis()->SetTitleSize(0.055);
h2t->GetYaxis()->SetLabelSize(0.055);





TH1F *h3t=new TH1F("h3ftt","Groomed jet mass in all n_PV for AK10 jets ",100,0,rangeAK10);
h3t->SetLineColor(5);
h3t->SetLineWidth(3);
h3t->GetXaxis()->SetTitle("Jet mass in GeV");
h3t->GetXaxis()->SetTitleSize(0.055);
h3t->GetXaxis()->SetLabelSize(0.055);
h3t->GetYaxis()->SetTitle("Normalized to Unit");
h3t->GetYaxis()->SetTitleSize(0.055);
h3t->GetYaxis()->SetLabelSize(0.055);



TH1F *h4t=new TH1F("h4trt","Trimmed jet mass",100,0,rangeAK10);
h4t->SetLineColor(3);
h4t->SetLineWidth(3);
h4t->GetXaxis()->SetTitle("Jet mass in GeV");
h4t->GetXaxis()->SetTitleSize(0.055);
h4t->GetXaxis()->SetLabelSize(0.055);
h4t->GetYaxis()->SetTitle("Normalized to Unit");
h4t->GetYaxis()->SetTitleSize(0.055);
h4t->GetYaxis()->SetLabelSize(0.055);



TH1F *h5t=new TH1F("h5prt","Prunned jet mass",100,0,rangeAK10);
h5t->SetLineColor(4);
h5t->SetLineWidth(3);
h5t->GetXaxis()->SetTitle("Jet mass in GeV");
h5t->GetXaxis()->SetTitleSize(0.055);
h5t->GetXaxis()->SetLabelSize(0.055);
h5t->GetYaxis()->SetTitle("Normalized to Unit");
h5t->GetYaxis()->SetTitleSize(0.055);
h5t->GetYaxis()->SetLabelSize(0.055);



TH1F *h6t=new TH1F("h6t","Soft drop #beta=0",100,0,rangeAK10);
h6t->SetLineColor(2);
h6t->SetLineWidth(3);
//h6->GetYaxis()->SetRangeUser(0,0.12);
h6t->GetXaxis()->SetTitle("Jet mass in GeV");
h6t->GetXaxis()->SetTitleSize(0.055);
h6t->GetXaxis()->SetLabelSize(0.055);
h6t->GetYaxis()->SetTitle("Normalized to Unit");
h6t->GetYaxis()->SetTitleSize(0.055);
h6t->GetYaxis()->SetLabelSize(0.055);



TH1F *h7t=new TH1F("h7t","Softdrop #beta=2",100,0,rangeAK10);
h7t->SetLineColor(6);
h7t->SetLineWidth(3);
h7t->GetXaxis()->SetTitle("Jet mass in GeV");
h7t->GetXaxis()->SetTitleSize(0.055);
h7t->GetXaxis()->SetLabelSize(0.055);
h7t->GetYaxis()->SetTitle("Normalized to Unit");
h7t->GetYaxis()->SetTitleSize(0.055);
h7t->GetYaxis()->SetLabelSize(0.055);



TH1F *h8t=new TH1F("h8t","Softdrop #beta=-1",100,0,rangeAK10);
h8t->SetLineColor(8);
h8t->SetLineWidth(3);
h8t->GetXaxis()->SetTitle("Jet mass in GeV");
h8t->GetXaxis()->SetTitleSize(0.055);
h8t->GetXaxis()->SetLabelSize(0.055);
h8t->GetYaxis()->SetTitle("Normalized to Unit");
h8t->GetYaxis()->SetTitleSize(0.055);
h8t->GetYaxis()->SetLabelSize(0.055);




TH1F *h9t=new TH1F("h9t","Softdrop #beta=1",100,0,rangeAK10);
h9t->SetLineColor(7);
h9t->SetLineWidth(3);
h9t->GetXaxis()->SetTitle("Jet mass in GeV");
h9t->GetXaxis()->SetTitleSize(0.055);
h9t->GetXaxis()->SetLabelSize(0.055);
h9t->GetYaxis()->SetTitle("Normalized to Unit");
h9t->GetYaxis()->SetTitleSize(0.055);
h9t->GetYaxis()->SetLabelSize(0.055);




TH1F *h10t=new TH1F("h10Ungt","Ungroomed jet mass",100,0,rangeAK10);
h10t->SetLineColor(9);
h10t->SetLineWidth(3);
h10t->GetXaxis()->SetTitle("Jet mass in GeV");
h10t->GetXaxis()->SetTitleSize(0.055);
h10t->GetXaxis()->SetLabelSize(0.055);
h10t->GetYaxis()->SetTitle("Events/ 3 GeV");
h10t->GetYaxis()->SetTitleSize(0.055);
h10t->GetYaxis()->SetLabelSize(0.055);


/////////////////////////////////////plots for AK10


//////////////////////////////////////////////plots for AK12

double rangeAK12=500;
//////////////////////////////////////plots for AK12
TH1F *h2te=new TH1F("h2te","Leading Jet Pt distribution, AK12 jets",250,0,2000);
h2te->SetLineWidth(3);
h2te->SetLineColor(2);
h2te->GetXaxis()->SetTitle("Jet pT in GeV");
h2te->GetXaxis()->SetTitleSize(0.055);
h2te->GetXaxis()->SetLabelSize(0.055);
h2te->GetYaxis()->SetTitle("Events/8 GeV");
h2te->GetYaxis()->SetTitleSize(0.055);
h2te->GetYaxis()->SetLabelSize(0.055);





TH1F *h3te=new TH1F("h3ftte","Groomed jet mass in all n_PV for AK12 jets ",100,0,rangeAK12);
h3te->SetLineColor(5);
h3te->SetLineWidth(3);
h3te->GetXaxis()->SetTitle("Jet mass in GeV");
h3te->GetXaxis()->SetTitleSize(0.055);
h3te->GetXaxis()->SetLabelSize(0.055);
h3te->GetYaxis()->SetTitle("Normalized to Unit");
h3te->GetYaxis()->SetTitleSize(0.055);
h3te->GetYaxis()->SetLabelSize(0.055);



TH1F *h4te=new TH1F("h4trte","Trimmed jet mass",100,0,rangeAK12);
h4te->SetLineColor(3);
h4te->SetLineWidth(3);
h4te->GetXaxis()->SetTitle("Jet mass in GeV");
h4te->GetXaxis()->SetTitleSize(0.055);
h4te->GetXaxis()->SetLabelSize(0.055);
h4te->GetYaxis()->SetTitle("Normalized to Unit");
h4te->GetYaxis()->SetTitleSize(0.055);
h4te->GetYaxis()->SetLabelSize(0.055);



TH1F *h5te=new TH1F("h5prte","Prunned jet mass",100,0,rangeAK12);
h5te->SetLineColor(4);
h5te->SetLineWidth(3);
h5te->GetXaxis()->SetTitle("Jet mass in GeV");
h5te->GetXaxis()->SetTitleSize(0.055);
h5te->GetXaxis()->SetLabelSize(0.055);
h5te->GetYaxis()->SetTitle("Normalized to Unit");
h5te->GetYaxis()->SetTitleSize(0.055);
h5te->GetYaxis()->SetLabelSize(0.055);



TH1F *h6te=new TH1F("h6te","Soft drop #beta=0",100,0,rangeAK12);
h6te->SetLineColor(2);
h6te->SetLineWidth(3);
//h6->GetYaxis()->SetRangeUser(0,0.12);
h6te->GetXaxis()->SetTitle("Jet mass in GeV");
h6te->GetXaxis()->SetTitleSize(0.055);
h6te->GetXaxis()->SetLabelSize(0.055);
h6te->GetYaxis()->SetTitle("Normalized to Unit");
h6te->GetYaxis()->SetTitleSize(0.055);
h6te->GetYaxis()->SetLabelSize(0.055);



TH1F *h7te=new TH1F("h7te","Softdrop #beta=2",100,0,rangeAK12);
h7te->SetLineColor(6);
h7te->SetLineWidth(3);
h7te->GetXaxis()->SetTitle("Jet mass in GeV");
h7te->GetXaxis()->SetTitleSize(0.055);
h7te->GetXaxis()->SetLabelSize(0.055);
h7te->GetYaxis()->SetTitle("Normalized to Unit");
h7te->GetYaxis()->SetTitleSize(0.055);
h7te->GetYaxis()->SetLabelSize(0.055);



TH1F *h8te=new TH1F("h8te","Softdrop #beta=-1",100,0,rangeAK12);
h8te->SetLineColor(8);
h8te->SetLineWidth(3);
h8te->GetXaxis()->SetTitle("Jet mass in GeV");
h8te->GetXaxis()->SetTitleSize(0.055);
h8te->GetXaxis()->SetLabelSize(0.055);
h8te->GetYaxis()->SetTitle("Normalized to Unit");
h8te->GetYaxis()->SetTitleSize(0.055);
h8te->GetYaxis()->SetLabelSize(0.055);




TH1F *h9te=new TH1F("h9te","Softdrop #beta=1",100,0,rangeAK12);
h9te->SetLineColor(7);
h9te->SetLineWidth(3);
h9te->GetXaxis()->SetTitle("Jet mass in GeV");
h9te->GetXaxis()->SetTitleSize(0.055);
h9te->GetXaxis()->SetLabelSize(0.055);
h9te->GetYaxis()->SetTitle("Normalized to Unit");
h9te->GetYaxis()->SetTitleSize(0.055);
h9te->GetYaxis()->SetLabelSize(0.055);




TH1F *h10te=new TH1F("h10Ungte","Ungroomed jet mass",100,0,rangeAK12);
h10te->SetLineColor(9);
h10te->SetLineWidth(3);
h10te->GetXaxis()->SetTitle("Jet mass in GeV");
h10te->GetXaxis()->SetTitleSize(0.055);
h10te->GetXaxis()->SetLabelSize(0.055);
h10te->GetYaxis()->SetTitle("Events/ 3 GeV");
h10te->GetYaxis()->SetTitleSize(0.055);
h10te->GetYaxis()->SetLabelSize(0.055);




//   In a ROOT session, you can do:
//      Root > .L AnalysisJetmasspt.C
//      Root > AnalysisJetmasspt t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

h1->Fill(n_PVertex);

//////////////////////////////////////AK8
h2->Fill(GroomedJet_AK0_8_pt_corr->at(0));

if(GroomedJet_AK0_8_pt_ft_corr->at(0) > JetptcutFormass){
h3->Fill(GroomedJet_AK0_8_mass_ft_corr->at(0));
}

if(GroomedJet_AK0_8_pt_tr_corr->at(0) > JetptcutFormass){
h4->Fill(GroomedJet_AK0_8_mass_tr_corr->at(0));
}

if(GroomedJet_AK0_8_pt_pr_corr->at(0) > JetptcutFormass){
h5->Fill(GroomedJet_AK0_8_mass_pr_corr->at(0));
}


if(GroomedJet_AK0_8_pt_sd1_corr->at(0) > JetptcutFormass){
h6->Fill(GroomedJet_AK0_8_mass_sd1_corr->at(0));
}

if(GroomedJet_AK0_8_pt_sd2_corr->at(0) > JetptcutFormass){
h7->Fill(GroomedJet_AK0_8_mass_sd2_corr->at(0));
}

if(GroomedJet_AK0_8_pt_sd3_corr->at(0) > JetptcutFormass){
for(int i=0;i<GroomedJet_AK0_8_pt_sd3_corr->size();i++){
if(i==0){
h8->Fill(GroomedJet_AK0_8_mass_sd3_corr->at(i));
}

}
}

if(GroomedJet_AK0_8_pt_sd4_corr->at(0) > JetptcutFormass){
h9->Fill(GroomedJet_AK0_8_mass_sd4_corr->at(0));
}

if(GroomedJet_AK0_8_pt_corr->at(0) > JetptcutFormass){
h10->Fill(GroomedJet_AK0_8_mass_corr->at(0));
}
/////////////////////////////////////////////AK8


///////////////////////////AK10

for(int ii=0;ii<GroomedJet_AK1_pt_corr->size();ii++){
if(ii==0){
h2t->Fill(GroomedJet_AK1_pt_corr->at(ii));
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_ft_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_ft_corr->at(ii) > JetptcutFormass){//pt cut
h3t->Fill(GroomedJet_AK1_mass_ft_corr->at(ii));
}//pt cut
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_tr_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_tr_corr->at(ii) > JetptcutFormass){//pt cut
h4t->Fill(GroomedJet_AK1_mass_tr_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_pr_corr->size();ii++){
if(ii==0){

if(GroomedJet_AK1_pt_pr_corr->at(ii) > JetptcutFormass){//pt cut
h5t->Fill(GroomedJet_AK1_mass_pr_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_sd1_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_sd1_corr->at(ii) > JetptcutFormass){//pt cut
h6t->Fill(GroomedJet_AK1_mass_sd1_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_sd2_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_sd2_corr->at(ii) > JetptcutFormass){//pt cut
h7t->Fill(GroomedJet_AK1_mass_sd2_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_pt_sd3_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_sd3_corr->at(ii) > JetptcutFormass){//pt cut
h8t->Fill(GroomedJet_AK1_mass_sd3_corr->at(ii));
}
}
}


for(int ii=0;ii<GroomedJet_AK1_pt_sd4_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_sd4_corr->at(ii) > JetptcutFormass){//pt cut
h9t->Fill(GroomedJet_AK1_mass_sd4_corr->at(ii));
}
}
}


for(int ii=0;ii<GroomedJet_AK1_pt_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_pt_corr->at(ii) > JetptcutFormass){//pt cut
h10t->Fill(GroomedJet_AK1_mass_corr->at(ii));
}
}
}






//////////////////////////AK10

/////////////////////////////AK12
for(int ii=0;ii<GroomedJet_AK1_2_pt_corr->size();ii++){
if(ii==0){
h2te->Fill(GroomedJet_AK1_2_pt_corr->at(ii));
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_ft_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_ft_corr->at(ii) > JetptcutFormass){//pt cut
h3te->Fill(GroomedJet_AK1_2_mass_ft_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_tr_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_tr_corr->at(ii) > JetptcutFormass){//pt cut
h4te->Fill(GroomedJet_AK1_2_mass_tr_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_pr_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_pr_corr->at(ii) > JetptcutFormass){//pt cut
h5te->Fill(GroomedJet_AK1_2_mass_pr_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_sd1_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_sd1_corr->at(ii) > JetptcutFormass){//pt cut
h6te->Fill(GroomedJet_AK1_2_mass_sd1_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_sd2_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_sd2_corr->at(ii) > JetptcutFormass){//pt cut
h7te->Fill(GroomedJet_AK1_2_mass_sd2_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_sd3_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_sd3_corr->at(ii) > JetptcutFormass){//pt cut
h8te->Fill(GroomedJet_AK1_2_mass_sd3_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_sd4_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_sd4_corr->at(ii) > JetptcutFormass){//pt cut
h9te->Fill(GroomedJet_AK1_2_mass_sd4_corr->at(ii));
}
}
}

for(int ii=0;ii<GroomedJet_AK1_2_pt_corr->size();ii++){
if(ii==0){
if(GroomedJet_AK1_2_pt_corr->at(ii) > JetptcutFormass){//pt cut
h10te->Fill(GroomedJet_AK1_2_mass_corr->at(ii));
}
}
}






/////////////////////////////AK12












   }//event loop


Double_t scale3 = 1/h3->Integral();
Double_t scale4= 1/h4->Integral();
Double_t scale5= 1/h5->Integral();
Double_t scale6= 1/h6->Integral();
Double_t scale7= 1/h7->Integral();
Double_t scale8= 1/h8->Integral();
Double_t scale9= 1/h9->Integral();
Double_t scale10= 1/h10->Integral();


h3->Scale(scale3);
h4->Scale(scale4);
h5->Scale(scale5);
h6->Scale(scale6);
h7->Scale(scale7);
h8->Scale(scale8);
h9->Scale(scale9);
h10->Scale(scale10);








c1->cd(1);
h1->Draw();

c1->cd(2);
h2->Draw();

c1->cd(3);

//h3->Draw();

c1->cd(3);
//h4->Draw("SAME");

c1->cd(3);
h5->Draw("SAME");

c1->cd(3);
h6->Draw("SAME");

c1->cd(3);
h7->Draw("SAME");

c1->cd(3);
h8->Draw("SAME");

c1->cd(3);
h9->Draw("SAME");

c1->cd(3);


h10->Draw("SAME");
leg = new TLegend(0.6,0.7,0.89,0.89);
//leg->AddEntry(h3,"Filtered","lp");
//leg->AddEntry(h4,"Trimmed","lp");
leg->AddEntry(h5,"Prunned","lp");
leg->AddEntry(h6,"Softdrop #beta=0(=MMDT) ","lp");
leg->AddEntry(h7,"Softdrop #beta=2","lp");
leg->AddEntry(h8,"Softdrop #beta=-1","lp");
leg->AddEntry(h9,"Softdrop #beta=1","lp");
leg->AddEntry(h10,"Ungroomed","lp");

leg->Draw();

///////////////////////////////////AK10

Double_t scale3t = 1/h3t->Integral();
Double_t scale4t= 1/h4t->Integral();
Double_t scale5t= 1/h5t->Integral();
Double_t scale6t= 1/h6t->Integral();
Double_t scale7t= 1/h7t->Integral();
Double_t scale8t= 1/h8t->Integral();
Double_t scale9t= 1/h9t->Integral();
Double_t scale10t= 1/h10t->Integral();


//h3t->Scale(scale3t);
//h4t->Scale(scale4t);
h5t->Scale(scale5t);
h6t->Scale(scale6t);
h7t->Scale(scale7t);
h8t->Scale(scale8t);
h9t->Scale(scale9t);
h10t->Scale(scale10t);










c2->cd(1);
h2t->Draw();

c2->cd(2);

//h3t->Draw();

c2->cd(2);
//h4t->Draw("SAME");

c2->cd(2);
h5t->Draw("SAME");

c2->cd(2);
h6t->Draw("SAME");

c2->cd(2);
h7t->Draw("SAME");

c2->cd(2);
h8t->Draw("SAME");

c2->cd(2);
h9t->Draw("SAME");

c2->cd(2);


h10t->Draw("SAME");
legt = new TLegend(0.6,0.7,0.89,0.89);
//legt->AddEntry(h3t,"Filtered","lp");
//legt->AddEntry(h4t,"Trimmed","lp");
legt->AddEntry(h5t,"Prunned","lp");
legt->AddEntry(h6t,"Softdrop #beta=0(=MMDT) ","lp");
legt->AddEntry(h7t,"Softdrop #beta=2","lp");
legt->AddEntry(h8t,"Softdrop #beta=-1","lp");
legt->AddEntry(h9t,"Softdrop #beta=1","lp");
legt->AddEntry(h10t,"Ungroomed","lp");

legt->Draw();








///////////////////////////////////////AK10



////////////////////////////AK12

Double_t scale3te = 1/h3te->Integral();
Double_t scale4te= 1/h4te->Integral();
Double_t scale5te= 1/h5te->Integral();
Double_t scale6te= 1/h6te->Integral();
Double_t scale7te= 1/h7te->Integral();
Double_t scale8te= 1/h8te->Integral();
Double_t scale9te= 1/h9te->Integral();
Double_t scale10te= 1/h10te->Integral();


h3te->Scale(scale3te);
h4te->Scale(scale4te);
h5te->Scale(scale5te);
h6te->Scale(scale6te);
h7te->Scale(scale7te);
h8te->Scale(scale8te);
h9te->Scale(scale9te);
h10te->Scale(scale10te);










c2->cd(3);
h2te->Draw();

c2->cd(4);

//h3te->Draw();

c2->cd(4);
//h4te->Draw("SAME");

c2->cd(4);
h5te->Draw("SAME");

c2->cd(4);
h6te->Draw("SAME");

c2->cd(4);
h7te->Draw("SAME");

c2->cd(4);
h8te->Draw("SAME");

c2->cd(4);
h9te->Draw("SAME");

c2->cd(4);


h10te->Draw("SAME");
legte = new TLegend(0.6,0.7,0.89,0.89);
//legte->AddEntry(h3te,"Filtered","lp");
//legte->AddEntry(h4te,"Trimmed","lp");
legte->AddEntry(h5te,"Prunned","lp");
legte->AddEntry(h6te,"Softdrop #beta=0(=MMDT) ","lp");
legte->AddEntry(h7te,"Softdrop #beta=2","lp");
legte->AddEntry(h8te,"Softdrop #beta=-1","lp");
legte->AddEntry(h9te,"Softdrop #beta=1","lp");
legte->AddEntry(h10te,"Ungroomed","lp");

legte->Draw();


///////////////////////////AK12








//c1->cd(4);
//h8->Draw();

f->Write();
















}//main
