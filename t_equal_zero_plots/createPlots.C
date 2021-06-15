void editHist(TH1F *h, TString title, TString xaxis, TString yaxis, double lineWidth, Int_t lineStyle, Color_t lineColor, double markerSize, Int_t markerStyle, Color_t markerColor, bool isbSat);

int createPlots()
{
  // ********* USE THIS FOR FIGURES 1, 3, 4, 5, 6 ********** // 
  //TString inputDir = "/sphenix/user/gregtom3/ephenix-sbu/analysis/gregory_matousek/vlad_dvmp_but_fast/t_equal_zero_plots/100M_test";
  // ******************************************************* //

  // ********* USE THIS FOR FIGURE 2 ********* // 
  TString inputDir = "/sphenix/user/gregtom3/ephenix-sbu/analysis/gregory_matousek/vlad_dvmp_but_fast/t_equal_zero_plots/100M_test_figure2";
  // ******************************************************* //

  TString outputDir = "./plots";
 
  int numFigure = 6; // 1, 2, 3, 4, 5 or 6

  // If 1 or 2 is selected, then pick either {vm, decay} = {1,1} or {3,2}
  // If 3,4,5,6 is selected, the below options do not matter
  int vm        = 1; // 1 = J/Psi
                     // 2 = Rho
                     // 3 = Phi
  int decay     = 1; // 1 = ee
                     // 2 = kaon
  
  TString smearSystem = "Handbook";
  TString tempvmStrings[3]={"J/#Psi","#rho","#phi"};
  TString tempvmStrings2[3]={"JPsi","Rho","Phi"};
  TString vmString = tempvmStrings[vm-1];
  TString vmString2 = tempvmStrings2[vm-1];
  TString decayString_arr[2] = {"ee","kaon"};
  TString decayString = decayString_arr[decay-1];
  int lumi = 10;
  // -------------------------------------------------
  
  TString figtitles[8] = {"F1T","F1L","F2T","F2L","F3","F4","F5","F6"};
  TString xpomtitles[2] = {"_xpomlow","_xpomhigh"};
  TFile *inFile[8][2];
  for(int i = 0 ; i < 8 ; i++) {
    for(int j = 0 ; j < 2 ; j++){
      inFile[i][j] = new TFile(Form("%s/%s%s.root",inputDir.Data(),figtitles[i].Data(),xpomtitles[j].Data()),"READ");
    }
  }

  // -------------------------------------------------
 
  //Figure3_Plot(inputDir,10,"Perfect");
  //return 0;
  
  if(numFigure==1||numFigure==2)
    {
      int idx = (numFigure-1)*2;
      // Grab the 4 relevant files for figure 1 or 2 plots
      TFile *f[4];
      f[0]=inFile[idx][0]; // T xlow
      f[1]=inFile[idx][1]; // T xhigh
      f[2]=inFile[idx+1][0]; // L xlow
      f[3]=inFile[idx+1][1]; // L xhigh
      
      // Grab the 16 relevant histograms
      // h[8][2] , where second index marks either xlow or xhigh
      TH1F *h[8][2];
      if(numFigure==1)
	{
	  h[0][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F1T_xpomlow",decayString.Data()));
	  h[1][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F1T_xpomlow",decayString.Data()));
	  h[2][0]=(TH1F*)f[0]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F1T_xpomlow",decayString.Data()));
	  h[3][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F1L_xpomlow",decayString.Data()));
	  h[4][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F1L_xpomlow",decayString.Data()));
	  h[5][0]=(TH1F*)f[2]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F1L_xpomlow",decayString.Data()));
	  h[6][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F1T_xpomlow",decayString.Data()));
	  h[7][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F1L_xpomlow",decayString.Data()));
	  h[0][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F1T_xpomhigh",decayString.Data()));
	  h[1][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F1T_xpomhigh",decayString.Data()));
	  h[2][1]=(TH1F*)f[1]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F1T_xpomhigh",decayString.Data()));
	  h[3][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F1L_xpomhigh",decayString.Data()));
	  h[4][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F1L_xpomhigh",decayString.Data()));
	  h[5][1]=(TH1F*)f[3]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F1L_xpomhigh",decayString.Data()));
	  h[6][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F1T_xpomhigh",decayString.Data()));
	  h[7][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F1L_xpomhigh",decayString.Data()));	
	}
      else
	{
	  h[0][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F2T_xpomlow",decayString.Data()));
	  h[1][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F2T_xpomlow",decayString.Data()));
	  h[2][0]=(TH1F*)f[0]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F2T_xpomlow",decayString.Data()));
	  h[3][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F2L_xpomlow",decayString.Data()));
	  h[4][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F2L_xpomlow",decayString.Data()));
	  h[5][0]=(TH1F*)f[2]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F2L_xpomlow",decayString.Data()));
	  h[6][0]=(TH1F*)f[0]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F2T_xpomlow",decayString.Data()));
	  h[7][0]=(TH1F*)f[2]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F2L_xpomlow",decayString.Data()));
	  h[0][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F2T_xpomhigh",decayString.Data()));
	  h[1][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F2T_xpomhigh",decayString.Data()));
	  h[2][1]=(TH1F*)f[1]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F2T_xpomhigh",decayString.Data()));
	  h[3][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bSat_F2L_xpomhigh",decayString.Data()));
	  h[4][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_eCa_"+vmString2+"_%s_bSat_F2L_xpomhigh",decayString.Data()));
	  h[5][1]=(TH1F*)f[3]->Get(Form("h_divide_eCa_ep_"+vmString2+"_%s_bSat_F2L_xpomhigh",decayString.Data()));
	  h[6][1]=(TH1F*)f[1]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F2T_xpomhigh",decayString.Data()));
	  h[7][1]=(TH1F*)f[3]->Get(Form("h_divide_eAu_ep_"+vmString2+"_%s_bNonSat_F2L_xpomhigh",decayString.Data()));	
	}  
 
      editHist(h[0][0],"","","",1,1,kBlack,2,33,kBlack,false);
      editHist(h[3][0],"","","",1,1,kBlue,2,33,kBlue,false);
      editHist(h[0][1],"","","",1,1,kBlack,2,27,kBlack,false);
      editHist(h[3][1],"","","",1,1,kBlue,2,27,kBlue,false);
      
      editHist(h[1][0],"","","",1,1,kBlack,1,20,kBlack,false);
      editHist(h[4][0],"","","",1,1,kBlue,1,20,kBlue,false);
      editHist(h[1][1],"","","",1,1,kBlack,1,24,kBlack,false);
      editHist(h[4][1],"","","",1,1,kBlue,1,24,kBlue,false);
      
      editHist(h[2][0],"","","",1,1,kBlack,1,21,kBlack,false);
      editHist(h[5][0],"","","",1,1,kBlue,1,21,kBlue,false);
      editHist(h[2][1],"","","",1,1,kBlack,1,25,kBlack,false);
      editHist(h[5][1],"","","",1,1,kBlue,1,25,kBlue,false);

      editHist(h[6][0],"","","",1,1,kBlack,2,33,kBlack,false);
      editHist(h[7][0],"","","",1,1,kBlue,2,33,kBlue,false);
      editHist(h[6][1],"","","",1,1,kBlack,2,27,kBlack,false);
      editHist(h[7][1],"","","",1,1,kBlue,2,27,kBlue,false);
      // Make dummy histograms to store large markers
      TH1F *hdummy[8][2];
      for(int i = 0 ; i < 8 ; i++){
	for(int j = 0 ; j < 2 ; j++){
	  hdummy[i][j]=new TH1F();}}

      editHist(hdummy[0][0],"","","",1,1,kBlack,3,33,kBlack,true);
      editHist(hdummy[3][0],"","","",1,1,kBlue,3,33,kBlue,true);
      editHist(hdummy[0][1],"","","",1,1,kBlack,3,27,kBlack,true);
      editHist(hdummy[3][1],"","","",1,1,kBlue,3,27,kBlue,true);
      
      editHist(hdummy[1][0],"","","",1,1,kBlack,2,20,kBlack,true);
      editHist(hdummy[4][0],"","","",1,1,kBlue,2,20,kBlue,true);
      editHist(hdummy[1][1],"","","",1,1,kBlack,2,24,kBlack,true);
      editHist(hdummy[4][1],"","","",1,1,kBlue,2,24,kBlue,true);
      
      editHist(hdummy[2][0],"","","",1,1,kBlack,2,21,kBlack,true);
      editHist(hdummy[5][0],"","","",1,1,kBlue,2,21,kBlue,true);
      editHist(hdummy[2][1],"","","",1,1,kBlack,2,25,kBlack,true);
      editHist(hdummy[5][1],"","","",1,1,kBlue,2,25,kBlue,true);
      
      editHist(hdummy[6][0],"","","",1,1,kBlack,3,33,kBlack,true);
      editHist(hdummy[7][0],"","","",1,1,kBlue,3,33,kBlue,true);
      editHist(hdummy[6][1],"","","",1,1,kBlack,3,27,kBlack,true);
      editHist(hdummy[7][1],"","","",1,1,kBlue,3,27,kBlue,true);
      // Histogram text positioning
      // {bSat_x1, bSat_x2, bSat_y1, bSat_y2,
      // bNonSat_x1, bNonSat_x2, bNonSat_y1, bNonSat_y2}
      double textSettings[8];
      if(vm==1)
	{
	  double temptextSettings[8]={5,32,1.4,1.95,5,32,1.4,1.95};
	  for(int i = 0 ; i < 8 ; i++){textSettings[i]=temptextSettings[i];}
	}
      else if(vm==2)
	{
	  double temptextSettings[8]={5,32,1.4,1.95,5,32,1.4,1.95};
	  for(int i = 0 ; i < 8 ; i++){textSettings[i]=temptextSettings[i];}
	}
      else if(vm==3)
	{
	  double temptextSettings[8]={5,32,1.4,1.95,5,32,1.4,1.95};
	  for(int i = 0 ; i < 8 ; i++){textSettings[i]=temptextSettings[i];}
	}
      
      
      // Setting up the Plots
      TCanvas *c[6]; double cwidth=800; double cheight=600;
      TString ct[6]={"c_bSat_eAu_ep","c_bSat_eAu_eCa","c_bSat_eCa_ep",
		     "c_bNonSat_eAu_ep","c_bSat_xlow","c_bSat_xhigh"};
      TLegend *l[6]; double lData[4]={0.12,0.9,0.9,0.99};
      TString lt[6][4];
      lt[0][0]=Form("%s (T) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[0][1]=Form("%s (L) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[0][2]=Form("%s (T) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[0][3]=Form("%s (L) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());

      lt[1][0]=Form("%s (T) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[1][1]=Form("%s (L) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[1][2]=Form("%s (T) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[1][3]=Form("%s (L) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}",vmString.Data());

      lt[2][0]=Form("%s (T) , A=40/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[2][1]=Form("%s (L) , A=40/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[2][2]=Form("%s (T) , A=40/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[2][3]=Form("%s (L) , A=40/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());

      lt[3][0]=Form("%s (T) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[3][1]=Form("%s (L) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[3][2]=Form("%s (T) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[3][3]=Form("%s (L) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      
      // 12/19/2020 I have swapped, for the below 8 lines, (T) <-> (L) 
      // I will comment again if I change it back
      // 1/18/2021 It has been changed back
      lt[4][0]=Form("%s (T) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[4][1]=Form("%s (L) , A=197/A=40, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[4][2]=Form("%s (T) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[4][3]=Form("%s (L) , A=197/A=1, low x_{#lower[-0.25]{#Rho}}",vmString.Data());

      lt[5][0]=Form("%s (T) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[5][1]=Form("%s (L) , A=197/A=40, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[5][2]=Form("%s (T) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());
      lt[5][3]=Form("%s (L) , A=197/A=1, high x_{#lower[-0.25]{#Rho}}",vmString.Data());

      TPaveText *pt[6]; 
      TString t[6];
      t[0]="Simulations from Sartre 1.33";
      t[1]=Form("Detector smearing based on %s",smearSystem.Data());
      t[2]=Form("Integrated Luminosity = %dfb^{-1}/A",lumi);
      t[3]=Form("IPSat #gamma^{*} + A -> %s + A",vmString.Data());
      t[4]=Form("IPNonSat #gamma^{*} + A -> %s + A",vmString.Data());
      if(numFigure==1)
	t[5]="0 < |t| < 0.001";
      else if(numFigure==2)
	t[5]="0 < |t| < 0.5";
      
      TString xlabel = "Q^{2}[GeV^{2}]";
      TString ylabel = "";
      if(numFigure==1)
	ylabel = "A_{1}^{-2}d#sigma^{A_{1}}/dt/A_{2}^{-2}d#sigma^{A_{2}}/dt";
      else
	ylabel = "A_{1}^{-2}#sigma^{A_{1}}/A_{2}^{-2}#sigma^{A_{2}}";
  
      for(int i = 0 ; i < 4; i++)
	{
	  c[i]=new TCanvas(ct[i],ct[i],cwidth,cheight);
	  gPad->SetLogx();
	  gPad->SetLeftMargin(0.12);
	  gPad->SetBottomMargin(0.12);
	  gStyle->SetOptStat(0);
	  gStyle->SetErrorX(0);
	  l[i] = new TLegend(lData[0],lData[1],lData[2],lData[3]);
	  l[i]->SetNColumns(2);
	  if(i==3)
	    {
	      l[i]->AddEntry(hdummy[6][0],lt[i][0],"p");
	      l[i]->AddEntry(hdummy[7][0],lt[i][1],"p");
	      l[i]->AddEntry(hdummy[6][1],lt[i][2],"p");
	      l[i]->AddEntry(hdummy[7][1],lt[i][3],"p");
	    }
	  else
	    {
	      l[i]->AddEntry(hdummy[i][0],lt[i][0],"p");
	      l[i]->AddEntry(hdummy[i+3][0],lt[i][1],"p");
	      l[i]->AddEntry(hdummy[i][1],lt[i][2],"p");
	      l[i]->AddEntry(hdummy[i+3][1],lt[i][3],"p");
	    }
	 
	  if(i==3)
	    {
	      pt[i]=new TPaveText(textSettings[4],textSettings[6],textSettings[5],textSettings[7]);
	      pt[i]->AddText(t[0]);
	      pt[i]->AddText(t[1]);
	      pt[i]->AddText(t[2]);
	      pt[i]->AddText(t[4]);
	      pt[i]->AddText(t[5]);
	    }
	  else
	    {
	      pt[i]=new TPaveText(textSettings[0],textSettings[2],textSettings[1],textSettings[3]);
	      pt[i]->AddText(t[0]);
	      pt[i]->AddText(t[1]);
	      pt[i]->AddText(t[2]);
	      pt[i]->AddText(t[3]);
	      pt[i]->AddText(t[5]);
	    }
	  pt[i]->SetFillStyle(0);
	  pt[i]->SetTextAlign(11);
	  pt[i]->SetBorderSize(0);
	  if(i==3)
	    {
	      h[6][0]->GetXaxis()->SetTitle(xlabel);
	      h[6][0]->GetYaxis()->SetTitle(ylabel);
	      h[6][0]->GetYaxis()->SetTitleOffset(0.85);
	      h[6][0]->Draw("E1");
	      h[7][0]->Draw("E1 same");
	      h[6][1]->Draw("E1 same");
	      h[7][1]->Draw("E1 same");
	    }
	  else
	    {
	      h[i][0]->GetXaxis()->SetTitle(xlabel);
	      h[i][0]->GetYaxis()->SetTitle(ylabel);
	      h[i][0]->GetYaxis()->SetTitleOffset(0.85);
	      h[i][0]->Draw("E1");
	      h[i+3][0]->Draw("E1 same");
	      h[i][1]->Draw("E1 same");
	      h[i+3][1]->Draw("E1 same");
	    }
	  l[i]->Draw("same");
	  pt[i]->Draw("same");
	  gPad->RedrawAxis();
	  c[i]->SaveAs(Form("%s/fig%d_%s_%s.pdf",outputDir.Data(),numFigure,ct[i].Data(),vmString2.Data()));  
	}

      // Figures for Vlad's Paper
      for(int i = 4 ; i < 6; i++)
	{
	  c[i]=new TCanvas(ct[i],ct[i],cwidth,cheight);
	  gPad->SetLogx();
	  gPad->SetLeftMargin(0.12);
	  gPad->SetBottomMargin(0.12);
	  gStyle->SetOptStat(0);
	  gStyle->SetErrorX(0);
	  l[i] = new TLegend(lData[0],lData[1],lData[2],lData[3]);
	  l[i]->SetNColumns(2);
	  if(i==4)
	    {
	      l[i]->AddEntry(hdummy[1][0],lt[i][0],"p");
	      l[i]->AddEntry(hdummy[4][0],lt[i][1],"p");
	      l[i]->AddEntry(hdummy[0][0],lt[i][2],"p");
	      l[i]->AddEntry(hdummy[3][0],lt[i][3],"p");
	    }
	  else
	    {
	      l[i]->AddEntry(hdummy[1][1],lt[i][0],"p");
	      l[i]->AddEntry(hdummy[4][1],lt[i][1],"p");
	      l[i]->AddEntry(hdummy[0][1],lt[i][2],"p");
	      l[i]->AddEntry(hdummy[3][1],lt[i][3],"p");
	    }
	 
	  pt[i]=new TPaveText(textSettings[0],textSettings[2],textSettings[1],textSettings[3]);
	  pt[i]->AddText(t[0]);
	  pt[i]->AddText(t[1]);
	  pt[i]->AddText(t[2]);
	  pt[i]->AddText(t[3]);
	  pt[i]->AddText(t[5]);
	  
	  pt[i]->SetFillStyle(0);
	  pt[i]->SetTextAlign(11);
	  pt[i]->SetBorderSize(0);
	  if(i==4)
	    {
	      h[1][0]->GetXaxis()->SetTitle(xlabel);
	      h[1][0]->GetYaxis()->SetTitle(ylabel);
	      h[1][0]->GetYaxis()->SetTitleOffset(0.85);
	      h[1][0]->Draw("E1");
	      h[4][0]->Draw("E1 same");
	      h[0][0]->Draw("E1 same");
	      h[3][0]->Draw("E1 same");
	    }
	  else
	    {
	      h[1][1]->GetXaxis()->SetTitle(xlabel);
	      h[1][1]->GetYaxis()->SetTitle(ylabel);
	      h[1][1]->GetYaxis()->SetTitleOffset(0.85);
	      h[1][1]->Draw("E1");
	      h[4][1]->Draw("E1 same");
	      h[0][1]->Draw("E1 same");
	      h[3][1]->Draw("E1 same");
	    }
	  l[i]->Draw("same");
	  pt[i]->Draw("same");
	  gPad->RedrawAxis();
	  c[i]->SaveAs(Form("%s/fig%d_%s_%s.pdf",outputDir.Data(),numFigure,ct[i].Data(),vmString2.Data()));  
	}






    }
  else // Figures 3 -6
    {
      int idx = numFigure+1;
      // Grab the 2 relevant files
      TFile *f[2];
      f[0]=inFile[idx][0]; // xlow
      f[1]=inFile[idx][1]; // xhigh
      
      
      // Grab the 30 relevant histograms
      // h[15][2] , where second index marks either xlow or xhigh
      // As for the first index...
      //   0 - 2 = eAu bSat
      //    0 = JPsi
      //    1 = Rho
      //    2 = Phi
      //   3 - 5 = eCa bSat
      //   6 - 8 = ep bSat
      //   9 - 11 = eAu bNonSat
      //   12 - 14 = ep bNonSat
      TString coltypes[15] = {"eAu","eAu","eAu",
			      "eCa","eCa","eCa",
			      "ep","ep","ep",
			      "eAu","eAu","eAu",
			      "ep","ep","ep"};
      TString vmtypes[15] = {"JPsi","Rho","Phi",
			     "JPsi","Rho","Phi",
			     "JPsi","Rho","Phi",
			     "JPsi","Rho","Phi",
			     "JPsi","Rho","Phi"}; 
      /*  TString vmtypes[15] = {"JPsi","Phi","Rho",
			     "JPsi","Phi","Rho",
			     "JPsi","Phi","Rho",
			     "JPsi","Phi","Rho",
			     "JPsi","Phi","Rho"};*/

      /*      TString decaytypes[15] = {"ee","kaon","ee",
				"ee","kaon","ee",
				"ee","kaon","ee",
				"ee","kaon","ee",
				"ee","kaon","ee"};*/

      TString decaytypes[15] = {"ee","ee","kaon",
				"ee","ee","kaon",
				"ee","ee","kaon",
				"ee","ee","kaon",
				"ee","ee","kaon"};

      TString modeltypes[15] = {"bSat","bSat","bSat",
				"bSat","bSat","bSat",
				"bSat","bSat","bSat",
				"bNonSat","bNonSat","bNonSat",
				"bNonSat","bNonSat","bNonSat"};
      TString figureString = Form("F%d",numFigure);
      double markers_xlow[15] = {25,25,25,
				 24,24,24,
				 27,27,27,
				 25,25,25,
				 27,27,27};

      double markers_xhigh[15] = {21,21,21,
				 20,20,20,
				 33,33,33,
				 21,21,21,
				 33,33,33};
      double marker_size[15] = {1,1,1,
				1,1,1,
				2,2,2,
				1,1,1,
				2,2,2};
      int colors[15] = {1,4,4,
			   1,4,4,
			   1,4,4,
			   1,4,4,
			   1,4,4};
      TH1F *h[15][2];
      for(int i = 0 ; i < 15 ; i++){
	for(int j = 0 ; j < 2 ; j++){
	  TString hname = Form("h_%s_%s_%s_%s_%s%s",coltypes[i].Data(),vmtypes[i].Data(),decaytypes[i].Data(),modeltypes[i].Data(),figureString.Data(),xpomtitles[j].Data());
	  h[i][j]=(TH1F*)f[j]->Get(hname);
	  if(j==0)
	    editHist(h[i][j],"","","",1,1,colors[i],marker_size[i],markers_xlow[i],colors[i],true);
	  else
	    editHist(h[i][j],"","","",1,1,colors[i],marker_size[i],markers_xhigh[i],colors[i],true);
	}
      }
      //      h[0][0]->Draw("hist");
      //      h[0][1]->Draw("same");
      //      return 0;
      // Make dummy histograms to store large markers
      // First Index
      //  0 - 2 = eAu
      //    0 = JPsi
      //    1 = Rho
      //    2 = Phi
      //  3 - 5 = eCa
      //  6 - 8 = ep
      // Then for bNonSat
      // Second index is either{xlow, xhigh}
      
      TH1F *hdummy[15][2];
      for(int i = 0 ; i < 15 ; i++)
	{
	  for(int j = 0 ; j < 2 ; j++)
	    {
	      hdummy[i][j]=new TH1F();
	      if(j==0)
		editHist(hdummy[i][j],"","","",1,1,colors[i],marker_size[i]+1,markers_xlow[i],colors[i],true);
	      else
		editHist(hdummy[i][j],"","","",1,1,colors[i],marker_size[i]+1,markers_xhigh[i],colors[i],true);
	    }
	}
      
      

      // Each of the 4 possible "Figure" will have 4 plots (xlow xhigh bSat bNonSat)
      // Each plot will have a location for the histogram text
      // Each plot will have a location for the legend

      double textSatSettings[8];
      double textNonSatSettings[8];
      if(numFigure==3){
	double TEMPtextSatSettings[8] = {3,30,2,100,3,30,2,100};
        double TEMPtextNonSatSettings[8] = {3,30,2,100,3,30,2,100};
	for(int i = 0 ; i < 8 ; i++){
	  textSatSettings[i]=TEMPtextSatSettings[i];
	  textNonSatSettings[i]=TEMPtextNonSatSettings[i];
	}
      }else if(numFigure==4){
	double TEMPtextSatSettings[8] = {1.2,15,pow(10,5),0.93*pow(10,7),3,30,20,1100};
        double TEMPtextNonSatSettings[8] = {1.2,15,pow(10,5),0.93*pow(10,7),3,30,20,1100};
	for(int i = 0 ; i < 8 ; i++){
	  textSatSettings[i]=TEMPtextSatSettings[i];
	  textNonSatSettings[i]=TEMPtextNonSatSettings[i];
	}
      }else if(numFigure==5){
	double TEMPtextSatSettings[8] = {5,35,1.5*pow(10,3),0.9*pow(10,5),5,35,1.5*pow(10,3),0.9*pow(10,5)};
        double TEMPtextNonSatSettings[8] = {5,35,1.5*pow(10,3),0.9*pow(10,5),5,35,1.5*pow(10,3),0.9*pow(10,5)};
	for(int i = 0 ; i < 8 ; i++){
	  textSatSettings[i]=TEMPtextSatSettings[i];
	  textNonSatSettings[i]=TEMPtextNonSatSettings[i];
	}
      }else{
	double TEMPtextSatSettings[8] = {4,35,1500,.8*pow(10,5),1.1,10,0.02,2};
        double TEMPtextNonSatSettings[8] = {6,35,1500,.8*pow(10,5),1.1,10,0.02,2};
	for(int i = 0 ; i < 8 ; i++){
	  textSatSettings[i]=TEMPtextSatSettings[i];
	  textNonSatSettings[i]=TEMPtextNonSatSettings[i];
	}
      }
      
      // Setting up the Canvases
      // Each "Figure" has 4 plots to be made
      TCanvas *c[4]; double cwidth=800; double cheight=600;
      TString ct[4]={"c_bSat_xlow","c_bSat_xhigh","c_bNonSat_xlow",
		     "c_bNonSat_xhigh"};
      

      // Initializing the ymin and ymax of plots
      double ymin[4][4] = { {1,1,1,1} ,
			    {10,10,10,10} ,
			    {pow(10,-2),pow(10,-2),pow(10,-2),pow(10,-2)},
			    {pow(10,-2),pow(10,-2),pow(10,-2),pow(10,-2)}};
      
      double ymax[4][4] ={ { pow(10,6), pow(10,6), pow(10,6), pow(10,6) },
			   { pow(10,7), pow(10,7), pow(10,7), pow(10,7) },
			   { pow(10,5), pow(10,5), pow(10,5), pow(10,5)},
			   { pow(10,5), pow(10,5), pow(10,5), pow(10,5) }};
      // TLegend and its input
      TLegend *l[4]; double lData[4]={0.12,0.9,0.9,0.99};
      for(int i = 0 ; i < 4 ; i++)
	{
	  l[i]=new TLegend(0.12,0.9,0.9,0.99);
	  if(i<2)
	    {
	      l[i]->SetNColumns(3);
	    }
	  else
	    {
	      l[i]->SetNColumns(2);
	    }
	}
      if(numFigure==3||numFigure==5)
	{
	  l[0]->AddEntry(hdummy[2][0],"#phi (L) , A=197","p");
	  l[0]->AddEntry(hdummy[5][0],"#phi (L) , A=40","p");
	  l[0]->AddEntry(hdummy[8][0],"#phi (L) , A=1","p");
	  l[0]->AddEntry(hdummy[0][0],"J/#Psi (L) , A=197","p");
	  l[0]->AddEntry(hdummy[3][0],"J/#Psi (L) , A=40","p");
	  l[0]->AddEntry(hdummy[6][0],"J/#Psi (L) , A=1","p");

	  l[1]->AddEntry(hdummy[2][1],"#phi (L) , A=197","p");
	  l[1]->AddEntry(hdummy[5][1],"#phi (L) , A=40","p");
	  l[1]->AddEntry(hdummy[8][1],"#phi (L) , A=1","p");
	  l[1]->AddEntry(hdummy[0][1],"J/#Psi (L) , A=197","p");
	  l[1]->AddEntry(hdummy[3][1],"J/#Psi (L) , A=40","p");
	  l[1]->AddEntry(hdummy[6][1],"J/#Psi (L) , A=1","p");
      
	  l[2]->AddEntry(hdummy[11][0],"#phi (L) , A=197","p");
	  l[2]->AddEntry(hdummy[14][0],"#phi (L) , A=1","p");
	  l[2]->AddEntry(hdummy[9][0],"J/#Psi (L) , A=197","p");
	  l[2]->AddEntry(hdummy[12][0],"J/#Psi (L) , A=1","p");

	  l[3]->AddEntry(hdummy[11][1],"#phi (L) , A=197","p");
	  l[3]->AddEntry(hdummy[14][1],"#phi (L) , A=1","p");
	  l[3]->AddEntry(hdummy[9][1],"J/#Psi (L) , A=197","p");
	  l[3]->AddEntry(hdummy[12][1],"J/#Psi (L) , A=1","p");
	}
      else
	{
	  l[0]->AddEntry(hdummy[2][0],"#phi (T) , A=197","p");
	  l[0]->AddEntry(hdummy[5][0],"#phi (T) , A=40","p");
	  l[0]->AddEntry(hdummy[8][0],"#phi (T) , A=1","p");
	  l[0]->AddEntry(hdummy[0][0],"J/#Psi (T) , A=197","p");
	  l[0]->AddEntry(hdummy[3][0],"J/#Psi (T) , A=40","p");
	  l[0]->AddEntry(hdummy[6][0],"J/#Psi (T) , A=1","p");

	  l[1]->AddEntry(hdummy[2][1],"#phi (T) , A=197","p");
	  l[1]->AddEntry(hdummy[5][1],"#phi (T) , A=40","p");
	  l[1]->AddEntry(hdummy[8][1],"#phi (T) , A=1","p");
	  l[1]->AddEntry(hdummy[0][1],"J/#Psi (T) , A=197","p");
	  l[1]->AddEntry(hdummy[3][1],"J/#Psi (T) , A=40","p");
	  l[1]->AddEntry(hdummy[6][1],"J/#Psi (T) , A=1","p");
      
	  l[2]->AddEntry(hdummy[11][0],"#phi (T) , A=197","p");
	  l[2]->AddEntry(hdummy[14][0],"#phi (T) , A=1","p");
	  l[2]->AddEntry(hdummy[9][0],"J/#Psi (T) , A=197","p");
	  l[2]->AddEntry(hdummy[12][0],"J/#Psi (T) , A=1","p");

	  l[3]->AddEntry(hdummy[11][1],"#phi (T) , A=197","p");
	  l[3]->AddEntry(hdummy[14][1],"#phi (T) , A=1","p");
	  l[3]->AddEntry(hdummy[9][1],"J/#Psi (T) , A=197","p");
	  l[3]->AddEntry(hdummy[12][1],"J/#Psi (T) , A=1","p");
	}
    
      TPaveText *pt[4]; 
      pt[0]=new TPaveText(textSatSettings[0],textSatSettings[2],textSatSettings[1],textSatSettings[3]);
      pt[1]=new TPaveText(textSatSettings[4],textSatSettings[6],textSatSettings[5],textSatSettings[7]);
      pt[2]=new TPaveText(textNonSatSettings[0],textNonSatSettings[2],textNonSatSettings[1],textNonSatSettings[3]);
      pt[3]=new TPaveText(textNonSatSettings[4],textNonSatSettings[6],textNonSatSettings[5],textNonSatSettings[7]);

      TString t[6];
      t[0]="Simulations from Sartre 1.33";
      t[1]=Form("Detector smearing based on %s",smearSystem.Data());
      t[2]=Form("Integrated Luminosity = %dfb^{-1}/A",lumi);
      t[3]="IPSat #gamma^{*} + A -> VM + A";
      t[4]="IPNonSat #gamma^{*} + A -> VM + A";
      t[5]="0 < |t| < 0.001";
      
      TString xlabel = "Q^{2}[GeV^{2}]";
      TString ylabel[4] = {"Q^{6}A^{-2}d#sigma^{#gamma A->V A}/dt [nb GeV^{4}]",
			   "Q^{8}A^{-2}d#sigma^{#gamma A->V A}/dt [nb GeV^{6}]",
			   "Q^{-2}A^{-4/3}d#sigma^{#gamma A->V A}/dt [nb GeV^{-4}]",
			   "Q^{0}A^{-4/3}d#sigma^{#gamma A->V A}/dt [nb GeV^{-2}]"};

      for(int i = 0 ; i < 4; i++)
	{
	  c[i]=new TCanvas(ct[i],ct[i],cwidth,cheight);
	  gPad->SetLeftMargin(0.12);
	  gPad->SetBottomMargin(0.12);
	  gStyle->SetOptStat(0);
	  gStyle->SetErrorX(0);
	  pt[i]->AddText(t[0]);
	  pt[i]->AddText(t[1]);
	  pt[i]->AddText(t[2]);
	  if(i<2)//bSat
	    {
	      pt[i]->AddText(t[3]);
	    }
	  else	    
	    {
	      pt[i]->AddText(t[4]);
	    }
	  pt[i]->AddText(t[5]);
	  pt[i]->SetFillStyle(0);
	  pt[i]->SetTextAlign(11);
	  pt[i]->SetBorderSize(0);

	  if(i==0){
	    h[0][0]->GetXaxis()->SetTitle(xlabel);
	    h[0][0]->GetYaxis()->SetTitle(ylabel[numFigure-3]);
	    h[0][0]->GetXaxis()->SetRangeUser(1,20);
	    h[0][0]->GetYaxis()->SetRangeUser(ymin[numFigure-3][i],ymax[numFigure-3][i]);
	    h[0][0]->GetYaxis()->SetTitleOffset(0.85); 
	    h[0][0]->Draw("E1");
	    h[2][0]->Draw("E1 same");
	    h[3][0]->Draw("E1 same");
	    h[5][0]->Draw("E1 same");
	    h[6][0]->Draw("E1 same");
	    h[8][0]->Draw("E1 same");
	    cout << h[0][0]->GetName() << endl;
	    cout << h[2][0]->GetName() << endl;
	    cout << h[3][0]->GetName() << endl;
	    cout << h[5][0]->GetName() << endl;
	    cout << h[6][0]->GetName() << endl;
	    cout << h[8][0]->GetName() << endl;
	  }else if(i==1){
	    h[0][1]->GetXaxis()->SetTitle(xlabel);
	    h[0][1]->GetYaxis()->SetTitle(ylabel[numFigure-3]);
	    h[0][1]->GetXaxis()->SetRangeUser(1,20);
	    h[0][1]->GetYaxis()->SetRangeUser(ymin[numFigure-3][i],ymax[numFigure-3][i]);
	    h[0][1]->GetYaxis()->SetTitleOffset(0.85);
	    h[0][1]->Draw("E1");
	    h[2][1]->Draw("E1 same");
	    h[3][1]->Draw("E1 same");
	    h[5][1]->Draw("E1 same");
	    h[6][1]->Draw("E1 same");
	    h[8][1]->Draw("E1 same");
	  }else if(i==2){
	    h[9][0]->GetXaxis()->SetTitle(xlabel);
	    h[9][0]->GetYaxis()->SetTitle(ylabel[numFigure-3]);
	    h[9][0]->GetXaxis()->SetRangeUser(1,20);
	    h[9][0]->GetYaxis()->SetRangeUser(ymin[numFigure-3][i],ymax[numFigure-3][i]);
	    h[9][0]->GetYaxis()->SetTitleOffset(0.85);
	    h[9][0]->Draw("E1");
	    h[11][0]->Draw("E1 same");
	    h[12][0]->Draw("E1 same");
	    h[14][0]->Draw("E1 same");
	  }else if(i==3){
	    h[9][1]->GetXaxis()->SetTitle(xlabel);
	    h[9][1]->GetYaxis()->SetTitle(ylabel[numFigure-3]);
	    h[9][1]->GetXaxis()->SetRangeUser(1,20);
	    h[9][1]->GetYaxis()->SetRangeUser(ymin[numFigure-3][i],ymax[numFigure-3][i]);
	    h[9][1]->GetYaxis()->SetTitleOffset(0.85);
	    h[9][1]->Draw("E1");
	    h[11][1]->Draw("E1 same");
	    h[12][1]->Draw("E1 same");
	    h[14][1]->Draw("E1 same");
	  }
	  l[i]->Draw("same");
	  pt[i]->Draw("same");
	  gPad->SetLogx();
	  gPad->SetLogy();
	  gPad->RedrawAxis();
	  c[i]->SaveAs(Form("%s/fig%d_%s.pdf",outputDir.Data(),numFigure,ct[i].Data()));  
	}
    }
  return 0;
}

void editHist(TH1F *h, TString title, TString xaxis, TString yaxis, double lineWidth, Int_t lineStyle, Color_t lineColor, double markerSize, Int_t markerStyle, Color_t markerColor, bool isbSat)
{
  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  h->SetLineWidth(lineWidth);
  h->SetLineStyle(lineStyle);
  h->SetLineColor(lineColor);
  h->SetMarkerSize(markerSize);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerColor(markerColor);

  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);

  h->GetXaxis()->SetTitleOffset(0.5);
  h->GetYaxis()->SetTitleOffset(0.65);
  
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  if(isbSat)
    {
      //  h->GetYaxis()->SetRangeUser(0,2);
    }
  else
    {
      h->GetYaxis()->SetRangeUser(0,1.5);
    }
}
