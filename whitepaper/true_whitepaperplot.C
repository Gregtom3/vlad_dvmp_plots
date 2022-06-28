int true_whitepaperplot()
{
  int N = 100000000;
  double lumi = 10; // fb^-1
  double lumi_save = lumi;
  int tsave = 4; // %-age  (0 -> 0.01 | 1 -> 0.05 | 2 -> 0.1 | 3 -> 0.2 | 4 -> Method L



  string tsmear[5] = {"#delta t = 0.01t", "#delta t = 0.05t", "#delta t = 0.1t" ,"#delta t = 0.2t", "#delta t ~ Method L"};
  string Mytsmear = tsmear[tsave];

  TString path_to_sims = "./data/";
  TString path_to_truth_sims = "./has_truth_data/";
  
  TFile *f[8];
  TFile *ftruth[8];
  TTree *t[8];
  TTree *ttruth[8];
  TH1F  *h[8][2];
  TH1F  *htruth[8][2]; 

  TGraphErrors *tge[8][2];
  TGraphErrors *tgetruth[8][2];
  
  TString fnames[8] = {"eAu_JPsi_ee_bSat",
		       "eAu_JPsi_ee_bSat",
		       "eAu_Phi_kaon_bSat",
		       "eAu_Phi_kaon_bSat",
		       "eAu_JPsi_ee_bNonSat",
		       "eAu_JPsi_ee_bNonSat",
		       "eAu_Phi_kaon_bNonSat",
		       "eAu_Phi_kaon_bNonSat"};

  double branching_ratio[8] = {0.0594,1.,0.0002954,0.489,0.0594,1.,0.0002954,0.489};
  TString hnames[8][2];
  TString gnames[8][2];
  int marker_colors[8][2] = {{9,9},{9,9},{9,9},{9,9},{1,1},{1,1},{1,1},{1,1}};
  int marker_styles[8][2] = {{21,20},{21,20},{21,20},{21,20},{25,24},{25,24},{25,24},{25,24}};

  double tmin = 0.0;
  double tmax = 0.18;
  const int bins = 50;
  double xsec[8];
  double xsec_temp=1.0;
  for(int i = 0 ; i < 8 ; i++)
    {
      hnames[i][0]=Form("h_whitepaper_%s_coherent",fnames[i].Data());
      hnames[i][1]=Form("h_whitepaper_%s_incoherent",fnames[i].Data());
      gnames[i][0]=Form("g_%s_coherent",fnames[i].Data());
      gnames[i][1]=Form("g_%s_incoherent",fnames[i].Data());
      f[i]=new TFile(Form("%s%s.root",path_to_sims.Data(),fnames[i].Data()),"READ");
      ftruth[i]=new TFile(Form("%s%s.root",path_to_truth_sims.Data(),fnames[i].Data()),"READ");

      h[i][0]=(TH1F*)f[i]->Get(Form("%s_%d",hnames[i][0].Data(),tsave));
      h[i][1]=(TH1F*)f[i]->Get(Form("%s_%d",hnames[i][1].Data(),tsave));

      htruth[i][0]=(TH1F*)ftruth[i]->Get(Form("%s_0",hnames[i][0].Data()));
      htruth[i][1]=(TH1F*)ftruth[i]->Get(Form("%s_0",hnames[i][1].Data()));
    
      tge[i][0]=new TGraphErrors(bins);
      tge[i][1]=new TGraphErrors(bins);

      tge[i][0]->GetYaxis()->SetTitle("d#sigma^{(e+Au->e'+Au'+VM)}/dt (nb/GeV^{2})");
      tge[i][1]->GetYaxis()->SetTitle("d#sigma^{(e+Au->e'+Au'+VM)}/dt (nb/GeV^{2})");
      tge[i][0]->GetXaxis()->SetTitle("|t| (GeV^{2})");
      tge[i][1]->GetXaxis()->SetTitle("|t| (GeV^{2})");
      tge[i][0]->SetMarkerStyle(marker_styles[i][0]);
      tge[i][1]->SetMarkerStyle(marker_styles[i][1]);
      tge[i][0]->SetMarkerColor(marker_colors[i][0]);
      tge[i][1]->SetMarkerColor(marker_colors[i][1]);

      htruth[i][0]->SetLineColor(marker_colors[i][0]);
      htruth[i][1]->SetLineColor(marker_colors[i][1]);
    }
  
  int nevents[8][2];
  double scale[8][2];

  TGraphErrors *tge_clone[8][2];
  for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  for(int k = 1; k <= bins ; k++)
	    {
	      tge[i][j]->SetPoint(k-1,h[i][j]->GetBinCenter(k),h[i][j]->GetBinContent(k));
	      tge[i][j]->SetPointError(k-1,0,h[i][j]->GetBinError(k));
	    }
	  tge[i][j]->GetXaxis()->SetRangeUser(tmin,tmax);
	  tge[i][j]->GetYaxis()->SetRangeUser(0.0001,100000.0);  //0.001, 300000.0
	  tge[i][j]->GetXaxis()->SetLabelSize(0.045);
	  tge[i][j]->GetXaxis()->SetTitleOffset(0.8);
	  tge[i][j]->GetXaxis()->SetTitleSize(0.06);
	  tge[i][j]->GetYaxis()->SetLabelSize(0.045);
	  tge[i][j]->GetYaxis()->SetTitleOffset(0.85);
	  tge[i][j]->GetYaxis()->SetTitleSize(0.06);
	  tge_clone[i][j]=(TGraphErrors*)tge[i][j]->Clone();
	  tge_clone[i][j]->SetMarkerSize(2.5);
	}
    }


  /* Draw only specific TGraphErrors from the Whitepaper */

  TCanvas *c[2];
  TPaveText *pt[2][2]; 
  TLegend *legend[2];
  for(int i = 0 ; i < 2 ; i++)
    {
      c[i] = new TCanvas(Form("c%d",i),Form("c%d",i),800,600);
      if(i==0)
	{
	  pt[i][0] = new TPaveText(0.015, 1.4*1200, 0.150, 1.4*30000.0);  //0.05,800,0.17,20000.0
	  pt[i][1] = new TPaveText(0.110, 1.5*150, 0.165, 1.5*30000.0);  //0.17,100,0.24,20000.0
	  
	}
      else
	{
	  pt[i][0] = new TPaveText(0.015, (1.4/1.5)*0.0012, 0.150, (1.4/1.5)*0.03);  //0.05,8000,0.17,200000.0
	  pt[i][1] = new TPaveText(0.110, 0.00015, 0.165, 0.03);  //0.17,1000,0.24,200000.0
	}
      legend[i] = new TLegend(.12,.9,.9,.99);
      pt[i][0]->SetFillStyle(0);
      pt[i][0]->SetTextAlign(11);
      pt[i][0]->SetBorderSize(0);
      pt[i][1]->SetFillStyle(0);
      pt[i][1]->SetTextAlign(11);
      pt[i][1]->SetBorderSize(0);
      pt[i][0]->AddText("Simulations from Sartre 1.33");
      pt[i][0]->AddText(Form("#scale[0.6]{#int} Ldt = %.0f fb^{-1}/A",lumi_save));
      pt[i][0]->AddText("1 < Q^{2} < 10 GeV^{2}");
      pt[i][1]->AddText("x < 0.01,  y > 0.01");
      if(i==0)
	{
	  pt[0][1]->AddText("|#eta(e_{#lower[-0.4]{decay}})| < 3.5");
	  pt[0][1]->AddText("p_{T}(e_{#lower[-0.4]{decay}}) > 0.3 GeV/c");
	}
      else
	{
	  pt[1][1]->AddText("|#eta(K_{#lower[-0.4]{decay}})| < 3.5");
	  pt[1][1]->AddText("p_{T}(K_{#lower[-0.4]{decay}}) > 0.3 GeV/c");
	}
      
  
      pt[i][1]->AddText(Form("%s", Mytsmear.c_str()));

      legend[i]->SetNColumns(2);
      if(i==0)
	{
	  legend[i]->AddEntry(tge_clone[0][0],"coherent - IPSat","P");
	  legend[i]->AddEntry(tge_clone[0][1],"incoherent - IPSat","P");
	  legend[i]->AddEntry(tge_clone[4][0],"coherent - IPNonSat","P");
	  legend[i]->AddEntry(tge_clone[4][1],"incoherent - IPNonSat","P");
	}
      else
	{
	  legend[i]->AddEntry(tge_clone[3][0],"coherent - IPSat","P");
	  legend[i]->AddEntry(tge_clone[3][1],"incoherent - IPSat","P");
	  legend[i]->AddEntry(tge_clone[7][0],"coherent - IPNonSat","P");
	  legend[i]->AddEntry(tge_clone[7][1],"incoherent - IPNonSat","P");
	}
   
      gPad->SetLogy();
      gPad->SetLeftMargin(0.12);
      gPad->SetBottomMargin(0.12);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0);
      if(i==0)
	{
	  tge[0][0]->GetYaxis()->SetTitle("d#sigma^{(e+Au#rightarrowe'+Au'+J/#Psi)}/dt  [nb/GeV^{ 2}]"); //new
	  tge[0][0]->GetYaxis()->SetRangeUser(0.0001,100000.0);  //0.001,30000.0
          tge[0][0]->GetXaxis()->SetTitle("|t| [GeV^{ 2}]");  //new
	  tge[0][0]->SetTitle("");
	  tge[0][0]->GetXaxis()->SetNdivisions(9,5,0);
	  tge[0][0]->Draw("AP");
	  tge[0][1]->Draw("P same");
	  htruth[0][0]->Draw("hist c same");
	  htruth[0][1]->Draw("hist c same");
	  tge[4][0]->Draw("P same");
	  tge[4][1]->Draw("P same");
	  htruth[4][0]->Draw("hist c same");
	  htruth[4][1]->Draw("hist c same");
	}
      else
	{
	  tge[3][0]->GetYaxis()->SetTitle("d#sigma^{(e+Au#rightarrowe'+Au'+#phi)}/dt  [nb/GeV^{ 2}]"); //new
          tge[3][0]->GetXaxis()->SetTitle("|t| [GeV^{ 2}]");  //new
	  tge[3][0]->SetTitle("");
	  tge[3][0]->GetXaxis()->SetNdivisions(9,5,0);
	  tge[3][0]->Draw("AP");
	  tge[3][1]->Draw("P same");
	  htruth[3][0]->Draw("hist c same");
	  htruth[3][1]->Draw("hist c same");
	  tge[7][0]->Draw("P same");
	  tge[7][1]->Draw("P same");
	  htruth[7][0]->Draw("hist c same");
	  htruth[7][1]->Draw("hist c same");	  
	}

      legend[i]->Draw("same");
      pt[i][1]->Draw("same");
      pt[i][0]->Draw("same");

      gPad->RedrawAxis();
      if(i==0)
	c[i]->SaveAs(Form("./plots/tsmear_%d_jpsi_lumi_%.0f.pdf",tsave,lumi_save));
      else
	c[i]->SaveAs(Form("./plots/tsmear_%d_phi_lumi_%.0f.pdf",tsave,lumi_save));
      //      c[i]->SaveAs(Form("./new_plots/t_jpsi_%.0f_%d.eps",lumi_original,i));
    }
  
  
  

  
  return 0;
}
