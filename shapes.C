  
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

int dolog=0;
void shapes() 
{ 
  char* atitle = "alpha2dsig ";
  char* hname1 = "ham2dfdk";
  char* hname2 = "ham2dfud";
  char* hname3 = "ham2dfb";
    //  char* lgd1 = "after preselection";
    //char* lgd2 = "after final selection";
  char* lgd1 = "dark quarks";
  char* lgd2 = "ud quarts";
  char* lgd3 = "b quarts";

  //TFile *f1 = new TFile("SumHistsQCD.root");
  //TFile *f1 = new TFile("SumHistsDebug.root");
  TFile *f2 = new TFile("SumHistsQQCD.root");
  //TFile *f1 = new TFile("SumHistsWMCtSkim.root");  
  //TFile *f1 = new TFile("SumHistsWSkim.root");  
  //TFile *f1 = new TFile("SumHists80.root");
  TFile *f1 = new TFile("SumHistsModelA.root");  
  //TFile *f1 = new TFile("SumHistsModelB.root");  
  //TFile *f1 = new TFile("SumHistsDATA.root");  
  //TFile *f1 = new TFile("SumHists80.root");  

 
  gStyle->SetOptStat(0);
 
  TString canvName = "Fig_";
  canvName += "hptdp_A_B";
  
  if( writeExtraText ) canvName += "-prelim";
  //if( iPos%10==0 ) canvName += "-out";
  //else if( iPos%10==1 ) canvName += "-left";
  //else if( iPos%10==2 )  canvName += "-center";
  //else if( iPos%10==3 )  canvName += "-right";
  int W = 800;
  int H = 600;
  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  // references for T, B, L, R
  float T = 0.08*H;
  float B = 0.12*H; 
  float L = 0.12*W;
  float R = 0.04*W;
  
  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  
  if (dolog) canv->SetLogy();

  //TH1* h_pt = new TH1F("h_pt"," ",100,0,500);
  //h_pt->GetXaxis()->SetNdivisions(6,5,0);
  //h_pt->GetXaxis()->SetTitle("Dark Pion p_{T} (GeV)");  
  //h_pt->GetXaxis()->SetTitleSize(0.05);  
  //h_pt->GetYaxis()->SetNdivisions(6,5,0);
  //h_pt->GetYaxis()->SetTitleOffset(1);
  //h_pt->GetYaxis()->SetTitle("Events / 5 GeV");  
  //h_pt->GetYaxis()->SetTitleSize(0.05);  
  
  //int max=  test->GetMaximum() + test->GetMaximum()*0.2; 
  //int max=  2000.;
  //h_pt->SetMaximum(max);
  //cout << max << endl;
 
  //h_pt->Draw();
  
  // int histLineColor = kOrange+7;
  //int histFillColor = kOrange-2;
  //float markerSize  = 1.0;

  TLatex latex;
  
  int n_ = 2;
  
  float x1_l = 1.0;
  //  float x1_l = 0.75;
  float y1_l = 0.90;
  
  float dx_l = 0.30;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  
 TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(32); lgd->SetFillColor(0);


  std::cout<<"getting first"<<std::endl;
  TH1F *A_pt = static_cast<TH1F*>(f1->Get(hname1)->Clone());
 A_pt->SetDirectory(0);
 //A_pt->Rebin(5);
  double aaA = A_pt->Integral();
std::cout<<" first entries is "<<aaA<<std::endl;
  A_pt->Scale(1./aaA);


  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt = static_cast<TH1F*>(f2->Get(hname2)->Clone());
  std::cout<<"ha"<<std::endl;
  B_pt->SetDirectory(0);
  //B_pt->Rebin(5);
  double aaB = B_pt->Integral();
std::cout<<" second entries is "<<aaB<<std::endl;
  B_pt->Scale(1/aaB);


  std::cout<<"getting third"<<std::endl;
  TH1F *C_pt = static_cast<TH1F*>(f2->Get(hname3)->Clone());
  std::cout<<"ha"<<std::endl;
  C_pt->SetDirectory(0);
  //C_pt->Rebin(5);
  double aaC = C_pt->Integral();
std::cout<<" third entries is "<<aaC<<std::endl;
  C_pt->Scale(1/aaC);


  double max = std::max(A_pt->GetMaximum(),B_pt->GetMaximum());
  max = std::max(max,C_pt->GetMaximum());
  A_pt->SetMaximum(max*1.3);

  A_pt->GetYaxis()->SetTitle(" percent  ");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  
  A_pt->GetXaxis()->SetTitle(atitle);  
  A_pt->GetXaxis()->SetTitleSize(0.05);  



  A_pt->SetLineColor(1);
  A_pt->SetLineWidth(3);
  A_pt->SetStats(0);
  A_pt->Draw("");

  

  B_pt->SetLineColor(2);
  B_pt->SetLineWidth(3);
  B_pt->SetStats(0);
  B_pt->Draw("same");

  C_pt->SetLineColor(3);
  C_pt->SetLineWidth(3);
  C_pt->SetStats(0);
  C_pt->Draw("same");


  //B_pt->Draw("esame");

  /*  
  std::cout<<"getting third"<<std::endl;
  TH1F *C_pt = static_cast<TH1F*>(f3->Get("haMgj")->Clone());
  C_pt->SetDirectory(0);
  double aaC = C_pt->Integral();
std::cout<<" third entries is "<<aaC<<std::endl;
C_pt->Scale(1/aaC);
  

  C_pt->SetLineColor(4);
  C_pt->SetLineWidth(3);
  C_pt->SetStats(0);
  
  C_pt->Draw("same");
  */


 
  //lgd->AddEntry(A_pt, "Monte Carlo QCD", "l");
  //lgd->AddEntry(B_pt, "Monte Carlo W to mu", "l");
 // lgd->AddEntry(C_pt, "data W to mu", "l");


  lgd->AddEntry(A_pt, lgd1, "l");
  lgd->AddEntry(B_pt, lgd2, "l");
  lgd->AddEntry(C_pt, lgd3, "l");
  //lgd->AddEntry(C_pt, "ModelBx500", "l");

 lgd->Draw();
    // Writing the lumi information and the CMS "logo"
   // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
   
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
   int iPos  = 11;
  CMS_lumi( canv, iPeriod, iPos );
  
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  lgd->Draw();

 
  if (dolog) {
    canv->Print(canvName+"_log.pdf",".pdf");
    canv->Print(canvName+"_log.png",".png");}
  else{ 
    canv->Print(canvName+".pdf",".pdf");
    canv->Print(canvName+".png",".png");}
  return;
}


