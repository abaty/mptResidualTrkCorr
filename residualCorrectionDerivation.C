#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"
#include "TCut.h"
#include "TStyle.h"

void residualCorrection()
{
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  const int najBins = 4;
  const int nptBins = 5;

  const double ajBins[najBins+1] = {0,0.11,0.22,0.33,1};
  const double ptBins[nptBins+1] = {0.5,1,2,4,8,300};

  TFile * f= new TFile("/export/d00/scratch/abaty/trackingEff/closure_ntuples/Correction_Vs3PF_ntuple_dijet_etaLT16_NoFFJEC.root","read");
  TTree * nt_track = (TTree*)f->Get("nt_track");
  TTree * nt_particle = (TTree*)f->Get("nt_particle");
  
  TH2D * leadCorr = new TH2D("leadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr = new TH2D("subleadCorr",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  TH2D * leadCorr_gen = new TH2D("leadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);
  TH2D * subleadCorr_gen = new TH2D("subleadCorr_gen",":A_{j}:p_{t}^{track}",najBins,ajBins,nptBins,ptBins);

  TCut cutReco_lead = "trackselect && pt>0.5 && r_lead<0.2";
  TCut cutReco_sublead = "trackselect && pt>0.5 && r_sublead<0.2";
  nt_track->Draw("pt:asym>>leadCorr","((1-fake)/eff)*weight"*(cutReco_lead),"colz");
  nt_track->Draw("pt:asym>>subleadCorr","((1-fake)/eff)*weight"*(cutReco_sublead),"colz");

  TCut cutGen_lead = "pt>0.5 && r_lead<0.2";
  TCut cutGen_sublead = "pt>0.5 && r_sublead<0.2";
  nt_particle->Draw("pt:asym>>leadCorr_gen","weight"*(cutGen_lead),"colz");
  nt_particle->Draw("pt:asym>>subleadCorr_gen","weight"*(cutGen_sublead),"colz");

  leadCorr->Divide(leadCorr_gen);
  subleadCorr->Divide(subleadCorr_gen);

  TH2D * diff = (TH2D*) leadCorr->Clone("diff");
  diff->Add(subleadCorr,-1);

  TLegend * l1 = new TLegend(0.37,0.5,0.77,0.9);
  l1->AddEntry((TObject*)0,"r_{lead}<0.2 Closure","");
  l1->SetBorderSize(0);  
  l1->SetTextSize(0.04);

  TLegend * l2 = new TLegend(0.37,0.5,0.77,0.9);
  l2->AddEntry((TObject*)0,"r_{sublead}<0.2 Closure","");
  l2->SetBorderSize(0);
  l2->SetTextSize(0.04);

  TLegend * l3 = new TLegend(0.37,0.5,0.77,0.9);
  l3->AddEntry((TObject*)0,"Difference in Closures","");
  l3->AddEntry((TObject*)0,"No Frag. JEC Applied","");
  l3->AddEntry((TObject*)0,"akVs3PF jets","");
  l3->SetBorderSize(0);
  l3->SetTextSize(0.04);

  TCanvas * c2 = new TCanvas("c2","c2",500,500);
  c2->SetRightMargin(0.15);
  c2->SetLogy();
  leadCorr->GetYaxis()->SetTitle("p_{t}");
  leadCorr->GetXaxis()->SetTitle("A_{j}");
  leadCorr->Draw("colz");
  l1->Draw("same");
  c2->SaveAs("plots/leading_correction_AkVs3PF_NoFFJEC.pdf");
  c2->SaveAs("plots/leading_correction_AkVs3PF_NoFFJEC.png");

  TCanvas * c3 = new TCanvas("c3","c3",500,500);
  c3->SetRightMargin(0.15);
  c3->SetLogy();
  subleadCorr->GetYaxis()->SetTitle("p_{t}");
  subleadCorr->GetXaxis()->SetTitle("A_{j}");
  subleadCorr->Draw("colz");
  l2->Draw("same");
  c3->SaveAs("plots/subleading_correction_AkVs3PF_NoFFJEC.pdf");
  c3->SaveAs("plots/subleading_correction_AkVs3PF_NoFFJEC.png");

  TCanvas * c4 = new TCanvas("c4","c4",500,500);
  c4->SetRightMargin(0.15);
  c4->SetLogy();
  diff->GetYaxis()->SetTitle("p_{t}");
  diff->GetXaxis()->SetTitle("A_{j}");
  diff->Draw("colz");
  l3->Draw("same");
  c4->SaveAs("plots/difference_AkVs3PF_NoFFJEC.pdf");
  c4->SaveAs("plots/difference_AkVs3PF_NoFFJEC.png");

  TFile * outf = new TFile("rootFiles/ResidualTrkCorr_AkVs3PF_NoFFJEC.root","recreate");
  leadCorr->SetDirectory(0);
  subleadCorr->SetDirectory(0);
  leadCorr->Write();
  subleadCorr->Write();
  diff->Write();
}
