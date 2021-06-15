enum markers {fullCircle=20,fullSquare=21,fullDiamond=33,openCircle=24,openSquare=25,openDiamond=27};


void editHist(TH1F *h, double lineWidth, Int_t lineStyle, Color_t lineColor, double markerSize, Int_t markerStyle, Color_t markerColor, bool isbSat);

int createPlots()
{
  double tmin = -0.004; double tmax = -0.003;
  //  double tmin = -0.007; double tmax = -0.006;
  //  double tmin = -0.010; double tmax = -0.009;
  //  double tmin = -0.013; double tmax = -0.012;
  //-------------------------------------------
  // Parameters
  //-------------------------------------------
  TString dir = Form("/sphenix/user/gregtom3/ephenix-sbu/analysis/gregory_matousek/vlad_dvmp_but_fast/scaling/100M_T_%0.3f/",-tmin);
  
  TString outputdir = "./plots/";

  TString smearSystem = "Handbook";

  double dt   = tmax - tmin;

  int lumi = 100; // fb^-1

  //-------------------------------------------
  // TFiles
  // F1L, F1T, F3
  // xpomhigh, xpomlow
  //-------------------------------------------

  TFile *F1L_xpomhigh = new TFile(Form("%stScale_F1L_xpomhigh.root",dir.Data()),"READ");
  TFile *F1L_xpomlow = new TFile(Form("%stScale_F1L_xpomlow.root",dir.Data()),"READ");
  TFile *F1T_xpomhigh = new TFile(Form("%stScale_F1T_xpomhigh.root",dir.Data()),"READ");
  TFile *F1T_xpomlow = new TFile(Form("%stScale_F1T_xpomlow.root",dir.Data()),"READ");
  TFile *F3_xpomhigh = new TFile(Form("%stScale_F3_xpomhigh.root",dir.Data()),"R'EAD");
  TFile *F3_xpomlow = new TFile(Form("%stScale_F3_xpomlow.root",dir.Data()),"READ");

 
  //-------------------------------------------
  // Plot 1
  // F1 Phi->kaon (bSat model, xpomhigh)
  //-------------------------------------------
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TH1F *h1L_eAu_eCa = (TH1F*)F1L_xpomhigh->Get("h_divide_eAu_eCa_Phi_kaon_bSat_F1L_xpomhigh");
  TH1F *h1L_eAu_ep  = (TH1F*)F1L_xpomhigh->Get("h_divide_eAu_ep_Phi_kaon_bSat_F1L_xpomhigh");
  TH1F *h1T_eAu_eCa = (TH1F*)F1T_xpomhigh->Get("h_divide_eAu_eCa_Phi_kaon_bSat_F1T_xpomhigh");
  TH1F *h1T_eAu_ep  = (TH1F*)F1T_xpomhigh->Get("h_divide_eAu_ep_Phi_kaon_bSat_F1T_xpomhigh");
  editHist(h1L_eAu_eCa,1,1,kBlue,1,openCircle,kBlue,false);
  editHist(h1L_eAu_ep,1,1,kBlue,2,openDiamond,kBlue,false);
  editHist(h1T_eAu_eCa,1,1,kBlack,1,openCircle,kBlack,false);
  editHist(h1T_eAu_ep,1,1,kBlack,2,openDiamond,kBlack,false);
 
  TH1F *h1[4];
  for(int i = 0 ; i < 4 ; i++)
    h1[i]=new TH1F();
  editHist(h1[0],1,1,kBlue,2,openCircle,kBlue,false);
  editHist(h1[1],1,1,kBlue,3,openDiamond,kBlue,false);
  editHist(h1[2],1,1,kBlack,2,openCircle,kBlack,false);
  editHist(h1[3],1,1,kBlack,3,openDiamond,kBlack,false);
  
  TString lt1[4]; // Legend text
  lt1[0] = "#phi (T) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}";
  lt1[1] = "#phi (L) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}";
  lt1[2] = "#phi (T) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}";
  lt1[3] = "#phi (L) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}";
  
  TLegend *l1 = new TLegend(0.12,0.9,0.9,0.99);
  l1->SetNColumns(2);
  l1->AddEntry(h1[2],lt1[0],"p");
  l1->AddEntry(h1[0],lt1[1],"p");
  l1->AddEntry(h1[3],lt1[2],"p");
  l1->AddEntry(h1[1],lt1[3],"p");

  TPaveText *pt1 = new TPaveText(5,0,20,0.4);
  pt1->AddText("Simulations from Sartre 1.33");
  pt1->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt1->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt1->AddText("IPSat #gamma^{*} + A -> #phi + A");
  pt1->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt1->SetFillStyle(0);
  pt1->SetTextAlign(11);
  pt1->SetBorderSize(0);

  TString xlabel1 = "Q^{2}[GeV^{2}]";
  TString ylabel1 = "G(A_{1},t)^{-1}d#sigma^{A_{1}}/dt/G(A_{2},t)^{-1}d#sigma^{A_{2}}/dt";

  gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h1L_eAu_eCa->GetXaxis()->SetTitle(xlabel1);
  h1L_eAu_eCa->GetYaxis()->SetTitle(ylabel1);
  h1L_eAu_eCa->Draw("E1");
  h1L_eAu_ep->Draw("E1 same");
  h1T_eAu_eCa->Draw("E1 same");
  h1T_eAu_ep->Draw("E1 same");
  l1->Draw("same");
  pt1->Draw("same");
  gPad->RedrawAxis();
  c1->SaveAs(Form("%sF1_phi_kaon_bSat_xhigh_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  // c1->Close();

 //-------------------------------------------
  // Plot 2
  // F1 Phi->kaon (bSat model, xpomlow)
  //-------------------------------------------
  
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  TH1F *h2L_eAu_eCa = (TH1F*)F1L_xpomlow->Get("h_divide_eAu_eCa_Phi_kaon_bSat_F1L_xpomlow");
  TH1F *h2L_eAu_ep  = (TH1F*)F1L_xpomlow->Get("h_divide_eAu_ep_Phi_kaon_bSat_F1L_xpomlow");
  TH1F *h2T_eAu_eCa = (TH1F*)F1T_xpomlow->Get("h_divide_eAu_eCa_Phi_kaon_bSat_F1T_xpomlow");
  TH1F *h2T_eAu_ep  = (TH1F*)F1T_xpomlow->Get("h_divide_eAu_ep_Phi_kaon_bSat_F1T_xpomlow");
  editHist(h2L_eAu_eCa,1,1,kBlue,1,fullCircle,kBlue,false);
  editHist(h2L_eAu_ep,1,1,kBlue,2,fullDiamond,kBlue,false);
  editHist(h2T_eAu_eCa,1,1,kBlack,1,fullCircle,kBlack,false);
  editHist(h2T_eAu_ep,1,1,kBlack,2,fullDiamond,kBlack,false);
 
  TH1F *h2[4];
  for(int i = 0 ; i < 4 ; i++)
    h2[i]=new TH1F();
  editHist(h2[0],1,1,kBlue,2,fullCircle,kBlue,false);
  editHist(h2[1],1,1,kBlue,3,fullDiamond,kBlue,false);
  editHist(h2[2],1,1,kBlack,2,fullCircle,kBlack,false);
  editHist(h2[3],1,1,kBlack,3,fullDiamond,kBlack,false);
  
  TString lt2[4]; // Legend text
  lt2[0] = "#phi (T) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}";
  lt2[1] = "#phi (L) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}";
  lt2[2] = "#phi (T) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}";
  lt2[3] = "#phi (L) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}";
  
  TLegend *l2 = new TLegend(0.12,0.9,0.9,0.99);
  l2->SetNColumns(2);
  l2->AddEntry(h2[2],lt2[0],"p");
  l2->AddEntry(h2[0],lt2[1],"p");
  l2->AddEntry(h2[3],lt2[2],"p");
  l2->AddEntry(h2[1],lt2[3],"p");

  TPaveText *pt2 = new TPaveText(5,32,1.4,1.95);
  pt2->AddText("Simulations from Sartre 1.33");
  pt2->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt2->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt2->AddText("IPSat #gamma^{*} + A -> #phi + A");
  pt2->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt2->SetFillStyle(0);
  pt2->SetTextAlign(11);
  pt2->SetBorderSize(0);

  TString xlabel2 = "Q^{2}[GeV^{2}]";
  TString ylabel2 = "G(A_{1},t)^{-1}d#sigma^{A_{1}}/dt/G(A_{2},t)^{-1}d#sigma^{A_{2}}/dt";

  gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h2L_eAu_eCa->GetXaxis()->SetTitle(xlabel2);
  h2L_eAu_eCa->GetYaxis()->SetTitle(ylabel2);
  h2L_eAu_eCa->Draw("E1");
  h2L_eAu_ep->Draw("E1 same");
  h2T_eAu_eCa->Draw("E1 same");
  h2T_eAu_ep->Draw("E1 same");
  l2->Draw("same");
  pt2->Draw("same");
  gPad->RedrawAxis();
  c2->SaveAs(Form("%sF1_phi_kaon_bSat_xlow_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  //  c2->Close();

  //-------------------------------------------
  // Plot 3
  // F1 JPsi->ee (bSat model, xpomhigh)
  //-------------------------------------------
  
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  TH1F *h3L_eAu_eCa = (TH1F*)F1L_xpomhigh->Get("h_divide_eAu_eCa_JPsi_ee_bSat_F1L_xpomhigh");
  TH1F *h3L_eAu_ep  = (TH1F*)F1L_xpomhigh->Get("h_divide_eAu_ep_JPsi_ee_bSat_F1L_xpomhigh");
  TH1F *h3T_eAu_eCa = (TH1F*)F1T_xpomhigh->Get("h_divide_eAu_eCa_JPsi_ee_bSat_F1T_xpomhigh");
  TH1F *h3T_eAu_ep  = (TH1F*)F1T_xpomhigh->Get("h_divide_eAu_ep_JPsi_ee_bSat_F1T_xpomhigh");
  editHist(h3L_eAu_eCa,1,1,kBlue,1,openCircle,kBlue,false);
  editHist(h3L_eAu_ep,1,1,kBlue,2,openDiamond,kBlue,false);
  editHist(h3T_eAu_eCa,1,1,kBlack,1,openCircle,kBlack,false);
  editHist(h3T_eAu_ep,1,1,kBlack,2,openDiamond,kBlack,false);
 
  TH1F *h3[4];
  for(int i = 0 ; i < 4 ; i++)
    h3[i]=new TH1F();
  editHist(h3[0],1,1,kBlue,2,openCircle,kBlue,false);
  editHist(h3[1],1,1,kBlue,3,openDiamond,kBlue,false);
  editHist(h3[2],1,1,kBlack,2,openCircle,kBlack,false);
  editHist(h3[3],1,1,kBlack,3,openDiamond,kBlack,false);
  
  TString lt3[4]; // Legend text
  lt3[0] = "J/#Psi (T) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}";
  lt3[1] = "J/#Psi (L) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}";
  lt3[2] = "J/#Psi (T) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}";
  lt3[3] = "J/#Psi (L) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}";
  
  TLegend *l3 = new TLegend(0.12,0.9,0.9,0.99);
  l3->SetNColumns(2);
  l3->AddEntry(h3[2],lt3[0],"p");
  l3->AddEntry(h3[0],lt3[1],"p");
  l3->AddEntry(h3[3],lt3[2],"p");
  l3->AddEntry(h3[1],lt3[3],"p");

  TPaveText *pt3 = new TPaveText(5,32,1.4,1.95);
  pt3->AddText("Simulations from Sartre 1.33");
  pt3->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt3->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt3->AddText("IPSat #gamma^{*} + A -> J/#Psi + A");
  pt3->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt3->SetFillStyle(0);
  pt3->SetTextAlign(11);
  pt3->SetBorderSize(0);

  TString xlabel3 = "Q^{2}[GeV^{2}]";
  TString ylabel3 = "G(A_{1},t)^{-1}d#sigma^{A_{1}}/dt/G(A_{2},t)^{-1}d#sigma^{A_{2}}/dt";

  gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h3L_eAu_eCa->GetXaxis()->SetTitle(xlabel3);
  h3L_eAu_eCa->GetYaxis()->SetTitle(ylabel3);
  h3L_eAu_eCa->Draw("E1");
  h3L_eAu_ep->Draw("E1 same");
  h3T_eAu_eCa->Draw("E1 same");
  h3T_eAu_ep->Draw("E1 same");
  l3->Draw("same");
  pt3->Draw("same");
  gPad->RedrawAxis();
  c3->SaveAs(Form("%sF1_JPsi_ee_bSat_xhigh_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  // c3->Close();


  //-------------------------------------------
  // Plot 4
  // F1 JPsi->ee (bSat model, xpomlow)
  //-------------------------------------------
  
  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  TH1F *h4L_eAu_eCa = (TH1F*)F1L_xpomlow->Get("h_divide_eAu_eCa_JPsi_ee_bSat_F1L_xpomlow");
  TH1F *h4L_eAu_ep  = (TH1F*)F1L_xpomlow->Get("h_divide_eAu_ep_JPsi_ee_bSat_F1L_xpomlow");
  TH1F *h4T_eAu_eCa = (TH1F*)F1T_xpomlow->Get("h_divide_eAu_eCa_JPsi_ee_bSat_F1T_xpomlow");
  TH1F *h4T_eAu_ep  = (TH1F*)F1T_xpomlow->Get("h_divide_eAu_ep_JPsi_ee_bSat_F1T_xpomlow");
  editHist(h4L_eAu_eCa,1,1,kBlue,1,openCircle,kBlue,false);
  editHist(h4L_eAu_ep,1,1,kBlue,2,openDiamond,kBlue,false);
  editHist(h4T_eAu_eCa,1,1,kBlack,1,openCircle,kBlack,false);
  editHist(h4T_eAu_ep,1,1,kBlack,2,openDiamond,kBlack,false);
 
  TH1F *h4[4];
  for(int i = 0 ; i < 4 ; i++)
    h4[i]=new TH1F();
  editHist(h4[0],1,1,kBlue,2,openCircle,kBlue,false);
  editHist(h4[1],1,1,kBlue,3,openDiamond,kBlue,false);
  editHist(h4[2],1,1,kBlack,2,openCircle,kBlack,false);
  editHist(h4[3],1,1,kBlack,3,openDiamond,kBlack,false);
  
  TString lt4[4]; // Legend text
  lt4[0] = "J/#Psi (T) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}";
  lt4[1] = "J/#Psi (L) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}";
  lt4[2] = "J/#Psi (T) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}";
  lt4[3] = "J/#Psi (L) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}";
  
  TLegend *l4 = new TLegend(0.12,0.9,0.9,0.99);
  l4->SetNColumns(2);
  l4->AddEntry(h4[2],lt4[0],"p");
  l4->AddEntry(h4[0],lt4[1],"p");
  l4->AddEntry(h4[3],lt4[2],"p");
  l4->AddEntry(h4[1],lt4[3],"p");

  TPaveText *pt4 = new TPaveText(5,32,1.4,1.95);
  pt4->AddText("Simulations from Sartre 1.33");
  pt4->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt4->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt4->AddText("IPSat #gamma^{*} + A -> J/#Psi + A");
  pt4->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt4->SetFillStyle(0);
  pt4->SetTextAlign(11);
  pt4->SetBorderSize(0);

  TString xlabel4 = "Q^{2}[GeV^{2}]";
  TString ylabel4 = "G(A_{1},t)^{-1}d#sigma^{A_{1}}/dt/G(A_{2},t)^{-1}d#sigma^{A_{2}}/dt";

  gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h4L_eAu_eCa->GetXaxis()->SetTitle(xlabel4);
  h4L_eAu_eCa->GetYaxis()->SetTitle(ylabel4);
  h4L_eAu_eCa->Draw("E1");
  h4L_eAu_ep->Draw("E1 same");
  h4T_eAu_eCa->Draw("E1 same");
  h4T_eAu_ep->Draw("E1 same");
  l4->Draw("same");
  pt4->Draw("same");
  gPad->RedrawAxis();
  c4->SaveAs(Form("%sF1_JPsi_ee_bSat_xlow_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  //  c4->Close();



 //-------------------------------------------
  // Plot 5
  // F3 (bSat model, xpomhigh)
  //-------------------------------------------
  
  TCanvas *c5 = new TCanvas("c5","c5",800,600);
  TH1F *h5_eAu_JPsi = (TH1F*)F3_xpomhigh->Get("h_tScale_eAu_JPsi_ee_bSat_F3_xpomhigh");
  TH1F *h5_eCa_JPsi = (TH1F*)F3_xpomhigh->Get("h_tScale_eCa_JPsi_ee_bSat_F3_xpomhigh");
  TH1F *h5_ep_JPsi = (TH1F*)F3_xpomhigh->Get("h_tScale_ep_JPsi_ee_bSat_F3_xpomhigh");
  TH1F *h5_eAu_Phi = (TH1F*)F3_xpomhigh->Get("h_tScale_eAu_Phi_kaon_bSat_F3_xpomhigh");
  TH1F *h5_eCa_Phi = (TH1F*)F3_xpomhigh->Get("h_tScale_eCa_Phi_kaon_bSat_F3_xpomhigh");
  TH1F *h5_ep_Phi = (TH1F*)F3_xpomhigh->Get("h_tScale_ep_Phi_kaon_bSat_F3_xpomhigh");
  editHist(h5_eAu_JPsi,1,1,kBlack,1,fullSquare,kBlack,true);
  editHist(h5_eCa_JPsi,1,1,kBlack,1,fullCircle,kBlack,true);
  editHist(h5_ep_JPsi,1,1,kBlack,2,fullDiamond,kBlack,true);
  editHist(h5_eAu_Phi,1,1,kBlue,1,fullSquare,kBlue,true);
  editHist(h5_eCa_Phi,1,1,kBlue,1,fullCircle,kBlue,true);
  editHist(h5_ep_Phi,1,1,kBlue,2,fullDiamond,kBlue,true);
 
  TH1F *h5[6];
  for(int i = 0 ; i < 6 ; i++)
    h5[i]=new TH1F();
  editHist(h5[0],1,1,kBlack,2,fullSquare,kBlack,true);
  editHist(h5[1],1,1,kBlack,2,fullCircle,kBlack,true);
  editHist(h5[2],1,1,kBlack,3,fullDiamond,kBlack,true);
  editHist(h5[3],1,1,kBlue,2,fullSquare,kBlue,true);
  editHist(h5[4],1,1,kBlue,2,fullCircle,kBlue,true);
  editHist(h5[5],1,1,kBlue,3,fullDiamond,kBlue,true);
 
  TString lt5[6]; // Legend text
  lt5[0] = "#phi (L) , A=197";
  lt5[1] = "#phi (L) , A=40";
  lt5[2] = "#phi (L) , A=1";
  lt5[3] = "J/#Psi (L) , A=197";
  lt5[4] = "J/#Psi (L) , A=40";
  lt5[5] = "J/#Psi (L) , A=1";

  
  TLegend *l5 = new TLegend(0.12,0.9,0.9,0.99);
  l5->SetNColumns(3);
  l5->AddEntry(h5[3],lt5[0],"p");
  l5->AddEntry(h5[4],lt5[1],"p");
  l5->AddEntry(h5[5],lt5[2],"p");
  l5->AddEntry(h5[0],lt5[3],"p");
  l5->AddEntry(h5[1],lt5[4],"p");
  l5->AddEntry(h5[2],lt5[5],"p");

  TPaveText *pt5 = new TPaveText(5,32,1.4,1.95);
  pt5->AddText("Simulations from Sartre 1.33");
  pt5->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt5->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt5->AddText("IPSat #gamma^{*} + A -> V + A");
  pt5->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt5->SetFillStyle(0);
  pt5->SetTextAlign(11);
  pt5->SetBorderSize(0);

  TString xlabel5 = "Q^{2}[GeV^{2}]";
  TString ylabel5 = "G(A,t)^{-1}Q^{6}d#sigma^{A}/dt [nb GeV^{4}]";

  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h5_eAu_JPsi->GetXaxis()->SetTitle(xlabel5);
  h5_eAu_JPsi->GetYaxis()->SetTitle(ylabel5);
  h5_eAu_JPsi->GetYaxis()->SetRangeUser(0.001,2*pow(10,3));
  h5_eAu_JPsi->Draw("E1");
  h5_eCa_JPsi->Draw("E1 same");
  h5_ep_JPsi->Draw("E1 same");
  h5_eAu_Phi->Draw("E1 same");
  h5_eCa_Phi->Draw("E1 same");
  h5_ep_Phi->Draw("E1 same");
  l5->Draw("same");
  pt5->Draw("same");
  gPad->RedrawAxis();
  c5->SaveAs(Form("%sF3_bSat_xhigh_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  //  c5->Close();

//-------------------------------------------
  // Plot 6
  // F3 (bSat model, xpomlow)
  //-------------------------------------------
  
  TCanvas *c6 = new TCanvas("c6","c6",800,600);
  TH1F *h6_eAu_JPsi = (TH1F*)F3_xpomlow->Get("h_tScale_eAu_JPsi_ee_bSat_F3_xpomlow");
  TH1F *h6_eCa_JPsi = (TH1F*)F3_xpomlow->Get("h_tScale_eCa_JPsi_ee_bSat_F3_xpomlow");
  TH1F *h6_ep_JPsi = (TH1F*)F3_xpomlow->Get("h_tScale_ep_JPsi_ee_bSat_F3_xpomlow");
  TH1F *h6_eAu_Phi = (TH1F*)F3_xpomlow->Get("h_tScale_eAu_Phi_kaon_bSat_F3_xpomlow");
  TH1F *h6_eCa_Phi = (TH1F*)F3_xpomlow->Get("h_tScale_eCa_Phi_kaon_bSat_F3_xpomlow");
  TH1F *h6_ep_Phi = (TH1F*)F3_xpomlow->Get("h_tScale_ep_Phi_kaon_bSat_F3_xpomlow");
  editHist(h6_eAu_JPsi,1,1,kBlack,1,openSquare,kBlack,true);
  editHist(h6_eCa_JPsi,1,1,kBlack,1,openCircle,kBlack,true);
  editHist(h6_ep_JPsi,1,1,kBlack,2,openDiamond,kBlack,true);
  editHist(h6_eAu_Phi,1,1,kBlue,1,openSquare,kBlue,true);
  editHist(h6_eCa_Phi,1,1,kBlue,1,openCircle,kBlue,true);
  editHist(h6_ep_Phi,1,1,kBlue,2,openDiamond,kBlue,true);
 
  TH1F *h6[6];
  for(int i = 0 ; i < 6 ; i++)
    h6[i]=new TH1F();
  editHist(h6[0],1,1,kBlack,2,openSquare,kBlack,true);
  editHist(h6[1],1,1,kBlack,2,openCircle,kBlack,true);
  editHist(h6[2],1,1,kBlack,3,openDiamond,kBlack,true);
  editHist(h6[3],1,1,kBlue,2,openSquare,kBlue,true);
  editHist(h6[4],1,1,kBlue,2,openCircle,kBlue,true);
  editHist(h6[5],1,1,kBlue,3,openDiamond,kBlue,true);
 
  TString lt6[6]; // Legend text
  lt6[0] = "#phi (L) , A=197";
  lt6[1] = "#phi (L) , A=40";
  lt6[2] = "#phi (L) , A=1";
  lt6[3] = "J/#Psi (L) , A=197";
  lt6[4] = "J/#Psi (L) , A=40";
  lt6[5] = "J/#Psi (L) , A=1";

  
  TLegend *l6 = new TLegend(0.12,0.9,0.9,0.99);
  l6->SetNColumns(3);
  l6->AddEntry(h6[3],lt6[0],"p");
  l6->AddEntry(h6[4],lt6[1],"p");
  l6->AddEntry(h6[5],lt6[2],"p");
  l6->AddEntry(h6[0],lt6[3],"p");
  l6->AddEntry(h6[1],lt6[4],"p");
  l6->AddEntry(h6[2],lt6[5],"p");

  TPaveText *pt6 = new TPaveText(5,32,1.4,1.95);
  pt6->AddText("Simulations from Sartre 1.33");
  pt6->AddText(Form("Detector smearing based on %s",smearSystem.Data()));
  pt6->AddText(Form("Integrated Luminosity = %dfb^{-1}/A",lumi));
  pt6->AddText("IPSat #gamma^{*} + A -> V + A");
  pt6->AddText(Form("%0.3f < |t| < %0.3f",-tmax,-tmin));
  pt6->SetFillStyle(0);
  pt6->SetTextAlign(11);
  pt6->SetBorderSize(0);

  TString xlabel6 = "Q^{2}[GeV^{2}]";
  TString ylabel6 = "G(A,t)^{-1}Q^{6}d#sigma^{A}/dt [nb GeV^{4}]";

  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  h6_eAu_JPsi->GetXaxis()->SetTitle(xlabel6);
  h6_eAu_JPsi->GetYaxis()->SetTitle(ylabel6);
  h6_eAu_JPsi->GetYaxis()->SetRangeUser(0.001,pow(10,3));
  h6_eAu_JPsi->Draw("E1");
  h6_eCa_JPsi->Draw("E1 same");
  h6_ep_JPsi->Draw("E1 same");
  h6_eAu_Phi->Draw("E1 same");
  h6_eCa_Phi->Draw("E1 same");
  h6_ep_Phi->Draw("E1 same");
  l6->Draw("same");
  pt6->Draw("same");
  gPad->RedrawAxis();
  c6->SaveAs(Form("%sF3_bSat_xlow_T_%.3f_%.3f.pdf",outputdir.Data(),-tmax,-tmin));
  //  c6->Close();
  

  return 0;
}




void editHist(TH1F *h, double lineWidth, Int_t lineStyle, Color_t lineColor, double markerSize, Int_t markerStyle, Color_t markerColor, bool isbSat)
{
  h->SetTitle("");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("");
  h->SetLineWidth(lineWidth);
  h->SetLineStyle(lineStyle);
  h->SetLineColor(lineColor);
  h->SetMarkerSize(markerSize);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerColor(markerColor);

  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);

  h->GetXaxis()->SetTitleOffset(0.5);
  h->GetYaxis()->SetTitleOffset(0.85);
  
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  if(isbSat)
    {
    }
  else
    {
      h->GetYaxis()->SetRangeUser(0,1.5);
    }
}
