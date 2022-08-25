

//==============================================================================
//  sartreEICTREE.cpp
//==============================================================================
#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TFile.h"
#include "Sartre.h"
#include "DipoleModelParameters.h"
#include "TGenPhaseSpace.h"
#include "Settings.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include <fstream>
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>
#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;

struct rootSartreEvent {
    double t;
    double Q2;
    double x;
    double s;
    double y;
    double W;
    double xpom;
    int    iEvent;
    int    pol;      // 0=transverse or 1=longitudinal
    int    dmode;    // 0=coherent, 1=Incoherent
};

double smear_E(double E, double eta, bool isElectron);
double photon_flux_T(double t, double xpom, double Q2, double s, double vm);
double photon_flux_L(double t, double xpom, double Q2, double s, double vm);
double differential(double t, double xpom, double Q2, double s, double vm);
double tScale(double A, double b0j, double t);

void dQ2(TH1F *h[3][2]);    // divide bin by dQ2
void dxpom(TH1F *h[3][2]);  // divide bin by dxpom
void dt(TH1F *h[3][2], double delta_t); // divide bin by dt
void dQ2_2(TH1F *h[8][2]);    // divide bin by dQ2
void dxpom_2(TH1F *h[8][2]);  // divide bin by dxpom
void dt_2(TH1F *h[8][2], double delta_t); // divide bin by dt
void scale_error(TH1F *h[3][2], double entries, double xsec, double scale, double BR); // scale histogram to lumi
void scale_error_2(TH1F *h[8][2], double entries, double xsec, double scale, double BR); // scale histogram to lumi

void log_bin(TH1F *h); // Log bin the x-axis

rootSartreEvent myRootSartreEvent;
void randomlyReverseBeams(Event* );  // for UPC mode
void addParticle(int index, int status, int id, double px, double py, double pz, double energy, double mass, std::ofstream * myfile);
int main(int argc, char *argv[])
{

  double lumi = 10.;

    //
    //  Check command line arguments
    //
    string runcard = argv[1];
    TString rootfile = TString(argv[2]);
    //    int nbatches = std::stoi(argv[4]);

    //
    //  Create the generator and initialize it.
    //  Once initialized you cannot (should not) change
    //  the settings w/o re-initialing sartre.
    //
    Sartre sartre;
    bool ok = sartre.init(runcard);
    if (!ok) {
        cerr << "Initialization of sartre failed." << endl;
        return 1;
    }

    EventGeneratorSettings* settings = sartre.runSettings();

    //
    //  ROOT file
    //  Use the one from command line unless not given
    //  in which case the one in the runcard is used.
    //
    string textfile;
    if (argc == 3) {
        textfile = argv[2];
        settings->setRootfile(argv[2]);
    }
    else
        textfile = settings->rootfile();
    
  
    //
    //  Print out all Sartre settings
    //

    settings->list();
    //
    // Prepare output file 
    //
    int A = settings->A();
    TString str_A = "";
    double b0j = 1.;
    if(A==1)
      {
	str_A="ep";
	b0j = 0.8783 * 5.0684;
	lumi = 100;
      }
    else if(A==40)
      {
	str_A="eCa";
	b0j = 1.0053 * 5.0684;
      }
    else if(A==197)
      {
	str_A="eAu";
	b0j = 1.203 * 5.0684;
      }
    lumi = lumi*1000000.;
    int VM = settings->vectorMesonId();
    double vmmass = 0.;
    double BR = 0.;
    TString str_vm = "";
    TString str_decay = "";
    bool isElectron=false;
    if(VM == 443)
      {
	str_vm = "JPsi";
	str_decay = "ee";
	vmmass = 3.096916;
	BR = 0.0594;
	isElectron=true;
      }
    else if(VM == 333)
      {
	str_vm = "Phi";
	str_decay = "kaon";
	vmmass = 1.019461;
	BR = 0.492;
      }
    TString str_model = TString(settings->dipoleModelName());
    TString htitle = Form("h_whitepaper_%s_%s_%s_%s",str_A.Data(),str_vm.Data(),str_decay.Data(),str_model.Data());

    
    TH1F *h[5][2];
    gROOT->cd();

    const int nbins = 50;
    const double tmax = 0.18;
    h[0][0] = new TH1F(Form("%s_coherent_0",htitle.Data()),Form("%s_coherent_0",htitle.Data()),nbins,0,tmax);
    h[1][0] = new TH1F(Form("%s_coherent_1",htitle.Data()),Form("%s_coherent_1",htitle.Data()),nbins,0,tmax);
    h[2][0] = new TH1F(Form("%s_coherent_2",htitle.Data()),Form("%s_coherent_2",htitle.Data()),nbins,0,tmax);
    h[3][0] = new TH1F(Form("%s_coherent_3",htitle.Data()),Form("%s_coherent_3",htitle.Data()),nbins,0,tmax);
    h[4][0] = new TH1F(Form("%s_coherent_4",htitle.Data()),Form("%s_coherent_4",htitle.Data()),nbins,0,tmax);

    h[0][1] = new TH1F(Form("%s_incoherent_0",htitle.Data()),Form("%s_incoherent_0",htitle.Data()),nbins,0,tmax);
    h[1][1] = new TH1F(Form("%s_incoherent_1",htitle.Data()),Form("%s_incoherent_1",htitle.Data()),nbins,0,tmax);
    h[2][1] = new TH1F(Form("%s_incoherent_2",htitle.Data()),Form("%s_incoherent_2",htitle.Data()),nbins,0,tmax);
    h[3][1] = new TH1F(Form("%s_incoherent_3",htitle.Data()),Form("%s_incoherent_3",htitle.Data()),nbins,0,tmax);
    h[4][1] = new TH1F(Form("%s_incoherent_4",htitle.Data()),Form("%s_incoherent_4",htitle.Data()),nbins,0,tmax);
  
    double tbin = tmax/nbins;
    //
    //  Setup ROOT tree
    //
    TLorentzVector *eIn = new TLorentzVector;
    TLorentzVector *pIn = new TLorentzVector;
    TLorentzVector *vm = new TLorentzVector;
    TLorentzVector *eOut = new TLorentzVector;
    TLorentzVector *pOut = new TLorentzVector;
    TLorentzVector *PomOut = new TLorentzVector;
    TLorentzVector *gamma = new TLorentzVector;
    TLorentzVector *vmDecay1 = new TLorentzVector;
    TLorentzVector *vmDecay2 = new TLorentzVector;
    
    TClonesArray protons("TLorentzVector");
    TClonesArray neutrons("TLorentzVector");
    TClonesArray remnants("TLorentzVector");
   
    


    //
    //  Prepare event generation
    //
    
    TGenPhaseSpace *decay = new TGenPhaseSpace();  // for VM decays
    int daughterID = settings->userInt();
    double daughterMasses[2] = {0, 0};
    bool doPerformDecay = false;
    if (daughterID && settings->vectorMesonId() != 22) {
        doPerformDecay = true;
	daughterMasses[0] = settings->lookupPDG(daughterID)->Mass();
        daughterMasses[1] = settings->lookupPDG(-daughterID)->Mass();
        cout << "Will decay vector meson: ";
        cout << settings->lookupPDG(settings->vectorMesonId())->GetName();
        cout << " -> ";
        cout << settings->lookupPDG(daughterID)->GetName();
        cout << " ";
        cout << settings->lookupPDG(-daughterID)->GetName();
        cout << endl;
    }
        

    
    
    //-------------------------------------------------------
    //  ** RECENTLY ADDED ** 
    // Calculate the total number
    // of events to run based on
    // the cross section of the event
    // and integrated luminosity
    // (Now not done anymore 3/8/2020)
    //-------------------------------------------------------
    unsigned long maxEvents = settings->numberOfEvents();
    double totalCS=sartre.totalCrossSection();

    double scale = lumi*totalCS/(A*1.)*BR/(1.*maxEvents);
    cout << "Generating " << maxEvents << " events." << endl << endl;
        int nPrint;
    if (settings->timesToShow())
        nPrint = maxEvents/settings->timesToShow();
    else
        nPrint = 0;

    //
    //  Event generation
    //
    for (unsigned long iEvent = 0; iEvent < maxEvents; iEvent++) {
        //
        //  Generate one event
        //
        Event *event = sartre.generateEvent();
        if (nPrint && (iEvent+1)%nPrint == 0 && iEvent != 0) {
            cout << "processed " << iEvent+1 << " events" << endl;
        }
        
        //
        //  If Sartre is run in UPC mode, half of the events needs to be
        //  rotated around and axis perpendicular to z:
        //
        if(settings->UPC() and settings->A()==settings->UPCA()){
            randomlyReverseBeams(event);
        }
   
        if (iEvent < 10) event->list();
  
	// 0 = transverse, 1 = longitudinal
       	int pol = event->polarization==transverse ? 0 : 1;
	// 0 = coherent, 1 = incoherent
	int diff = event->diffractiveMode==coherent ? 0 : 1;

        eIn     = &event->particles[0].p;
        pIn     = &event->particles[1].p;	
        eOut    = &event->particles[2].p;
        pOut    = &event->particles[6].p;
        vm      = &event->particles[4].p;
        gamma   = &event->particles[3].p;
	PomOut = &event->particles[5].p;

	TLorentzVector *vmDaughter1 = 0;
	TLorentzVector *vmDaughter2 = 0;
	if (doPerformDecay)
	  {
	    if (decay->SetDecay(*vm, 2, daughterMasses))
	      {
		double weight = decay->Generate();  // weight is always 1 here
		if ((weight - 1) > FLT_EPSILON)
		  {
		    cout << "PHSartre: Warning decay weight != 1, weight = " << weight << endl;
		  }
     
		vmDaughter1 = decay->GetDecay(0);
	        vmDaughter2 = decay->GetDecay(1);
		
		event->particles[4].status = 2;  // set VM status
		
		Particle vmDC1;
		vmDC1.index = event->particles.size();
		vmDC1.pdgId = daughterID;
		vmDC1.status = 1;  // final state
		vmDC1.p = *vmDaughter1;
		vmDC1.parents.push_back(4);
		event->particles.push_back(vmDC1);
		vmDecay1 = &event->particles[event->particles.size() - 1].p;
		
		Particle vmDC2;
		vmDC2.index = event->particles.size();
		vmDC2.pdgId = -daughterID;
		vmDC2.status = 1;  // final state
		vmDC2.p = *vmDaughter2;
		vmDC2.parents.push_back(4);
		event->particles.push_back(vmDC2);
		vmDecay2 = &event->particles[event->particles.size() - 1].p;
	      }
	    else
	      {
		cout << "Sartre: WARNING: Kinematics of Vector Meson does not allow decay!" << endl;
	      }
	  }


	double t = -event->t;
	double t_0 = gRandom->Gaus(t, 0.0 * t);
	double t_1 = gRandom->Gaus(t, 0.05 * t);
	double t_2 = gRandom->Gaus(t, 0.1 * t);
	double t_3 = gRandom->Gaus(t, 0.2 * t);
	double t_4 = 0.0;

	double delta_t = 0.0;
	if(t>=0&&t<0.01)
	  {
	    delta_t = 0.048;
	  }
	else if(t>=0.01&&t<0.04)
	  {
	    delta_t = 0.016;
	  }
	else if(t>=0.04&&t<0.07)
	  {
	    delta_t = 0.009;
	  }
	else if(t>=0.07&&t<0.10)
	  {
	    delta_t = 0.007;
	  }
	else if(t>=0.10&&t<0.13)
	  {
	    delta_t = 0.006;
	  }
	else if(t>=0.13)
	  {
	    delta_t = 0.005;
	  }
	bool fail = false;
	t_4 = gRandom->Gaus(t,delta_t*t);

	double s = event->s;


	double eOut_E = eOut->E();
	double eOut_P = eOut->P();
	double eOut_Theta = eOut->Theta();
	if(eOut_Theta<0.)
	  eOut_Theta+=TMath::Pi();

	double eOut_Eta = -log((tan(eOut_Theta/2.)));

	double vmDaughter1_E = vmDaughter1->E();
	double vmDaughter1_P = vmDaughter1->P();
	double vmDaughter1_Theta = vmDaughter1->Theta();
	if(vmDaughter1_Theta<0.)
	  vmDaughter1_Theta+=TMath::Pi();

	double vmDaughter1_Eta = -log((tan(vmDaughter1_Theta/2.)));

	double vmDaughter2_E = vmDaughter2->E();
	double vmDaughter2_P = vmDaughter2->P();
	double vmDaughter2_Theta = vmDaughter2->Theta();
	if(vmDaughter2_Theta<0.)
	  vmDaughter2_Theta+=TMath::Pi();

	double vmDaughter2_Eta = -log((tan(vmDaughter2_Theta/2.)));
  

	double eOut_E_reco = 0.;
	double eOut_Theta_reco = 0.;
	double eOut_Eta_reco = 0.;

	double vmDaughter1_E_reco = 0.;
	double vmDaughter1_Theta_reco = 0.;
	double vmDaughter1_Eta_reco = 0.;

	double vmDaughter2_E_reco = 0.;
	double vmDaughter2_Theta_reco = 0.;
	double vmDaughter2_Eta_reco = 0.;

	/*


	// ----------
	// Smear Theta
	// ----------

	if(abs(eOut_Eta)<=3.5)
	  eOut_Theta_reco = gRandom->Gaus(eOut_Theta,0.001);
	else
	  fail = true;
  
	if(abs(vmDaughter1_Eta)<=3.5)
	  vmDaughter1_Theta_reco = gRandom->Gaus(vmDaughter1_Theta,0.001);
	else
	  fail = true;

	if(abs(vmDaughter2_Eta)<=3.5)
	  vmDaughter2_Theta_reco = gRandom->Gaus(vmDaughter2_Theta,0.001);
	else
	  fail = true;

	// ----------
	// Smear Energy
	// ----------

	eOut_E_reco = smear_E(eOut_E,eOut_Eta,true);
	vmDaughter1_E_reco = smear_E(vmDaughter1_E,vmDaughter1_Eta,isElectron);
	vmDaughter2_E_reco = smear_E(vmDaughter2_E,vmDaughter2_Eta,isElectron);
  
	if(eOut_E_reco==-1||vmDaughter1_E_reco==-1||vmDaughter2_E_reco==-1)
	  fail = true;
	
	if(eOut->Pt()<0.3||vmDaughter1->Pt()<0.3||vmDaughter2->Pt()<0.3)
	  fail = true;
	

	double y = 1.0 - (eOut_E_reco/(2.0*eIn->E()))*(1.0-cos(eOut_Theta_reco));
	double Q2 = 2.0*eOut_E_reco*eIn->E()*(1.+cos(eOut_Theta_reco));
	double x = eOut_E_reco*(1.+cos(eOut_Theta_reco))/(2.*y*pIn->E());
	double W2 = pow(0.938272,2)+Q2*(1./x-1.);
	// t > 0 so we do "Q2 + t" as opposed to "Q2 - t"
	double xpom = (vmmass*vmmass+Q2+t)/(W2+Q2-pow(0.938272,2));
	*/

	double y = 1.0 - (eOut_E/(2.0*eIn->E()))*(1.0-cos(eOut_Theta));
	double x = eOut_E*(1.+cos(eOut_Theta))/(2.*y*pIn->E());
	// CAN REMOVE
	// Added 2/1/2021
	if(event->y<0.01)
	  fail = true;
	if(event->x>0.01)
	  fail=true;
	// Skip event if cuts aren't made
	if(fail==true)
	  continue;

	if(diff==0)
	  {
	    h[0][0]->Fill(t_0,scale);
	    h[1][0]->Fill(t_1,scale);
	    h[2][0]->Fill(t_2,scale);
	    h[3][0]->Fill(t_3,scale);
	    h[4][0]->Fill(t_4,scale);
	  }
	else if(diff==1)
	  {
	    h[0][1]->Fill(t_0,scale);
	    h[1][1]->Fill(t_1,scale);
	    h[2][1]->Fill(t_2,scale);
	    h[3][1]->Fill(t_3,scale);
	    h[4][1]->Fill(t_4,scale);
	  }
	
    }
    
    cout << "All events processed\n" << endl;
    
    for(int i = 0 ; i < 5; i++)
      {
	for(int j = 0 ; j < 2 ; j++)
	  {
	    for(int k = 1 ; k <= nbins ; k++)
	      {
		h[i][j]->SetBinError(k,sqrt(h[i][j]->GetBinContent(k)));
	      }
	    h[i][j]->Scale(A/lumi/BR/tbin);
	  }
      }
    
    cout << "Total cross-section: " << totalCS << " nb" << endl;
 
    TFile *fout = new TFile(Form("%s/%s_%s_%s_%s.root",rootfile.Data(),str_A.Data(),str_vm.Data(),str_decay.Data(),str_model.Data()),"RECREATE");
    for(int i = 0 ; i < 5 ; i++)
      {
	h[i][0]->Write();
	h[i][1]->Write();
      }

    fout->Close();
    return 0;   
}   

// UPC only
void randomlyReverseBeams(Event* myEvent) {
    
    TRandom3 *random = EventGeneratorSettings::randomGenerator();
    
    if(random->Uniform(1) > 0.5){
        for(unsigned int i=0; i<myEvent->particles.size(); i++)
            myEvent->particles.at(i).p.RotateX(M_PI);
    }
}

void addParticle(int index, int status, int id, double px, double py, double pz, double energy, double mass, std::ofstream * myfile)
{
  *myfile << index;
  *myfile << "\t";
  
  *myfile << status;
  *myfile << "\t";
  
  *myfile << id;
  *myfile << "\t";
  
  if(index==4) // Pythia ASCII must denote scattered lepton with '3'
    {
      *myfile << "3\t 0\t 0\t";
    }
  else
    {
      *myfile << "0\t 0\t 0\t";
    }
  
  *myfile << px;
  *myfile << "\t";

  *myfile << py;
  *myfile << "\t";

  *myfile << pz;
  *myfile << "\t";

  *myfile << energy;
  *myfile << "\t";

  *myfile << mass;
  *myfile << "\t";

  *myfile << "0\t 0\t 0\n";

}



double differential(double t, double xpom, double Q2, double s, double vm)
{
  return(vm*vm+Q2-t)/(s*xpom*xpom);
}
double photon_flux_T(double t, double xpom, double Q2, double s, double vm)
{
  double alpha = 0.0072973525693;
  double diff = differential(t,xpom,Q2,s,vm);
  double y = (vm*vm+Q2-t)/(s*xpom);
  double gamma = alpha*(1.0+(1.0-y)*(1.0-y))/(2.0*TMath::Pi()*y*Q2);
  return (1.0)/(diff*gamma);
}
double photon_flux_L(double t, double xpom, double Q2, double s, double vm)
{
  double alpha = 0.0072973525693;
  double diff = differential(t,xpom,Q2,s,vm);
  double y = (vm*vm+Q2-t)/(s*xpom);
  double gamma = alpha*(2.0*(1.0-y))/(2.0*TMath::Pi()*y*Q2);
  return (1.0)/(diff*gamma);
}

void log_bin(TH1F *h)
{
  TAxis *xaxis = h->GetXaxis();
  int xbins = xaxis->GetNbins();
  Axis_t xmin = xaxis->GetXmin();
  Axis_t xmax = xaxis->GetXmax();
  Axis_t xwidth = (xmax - xmin) / xbins;
  Axis_t *new_xbins = new Axis_t[xbins+1];

  for(int i = 0 ; i <= xbins; i++)
    {
      new_xbins[i] = TMath::Power( 10 , xmin+i*xwidth );
    }
  xaxis->Set(xbins,new_xbins);
}

double tScale(double A, double b0j, double t)
{
  double delta = TMath::Power(t,0.5);
  double numerator = 2.0 - 2.0*TMath::Cos(TMath::Power(A,1.0/3.0)*delta*b0j) - 2.0*TMath::Power(A,1.0/3.0)*TMath::Sin(TMath::Power(A,1.0/3.0)*delta*b0j)*delta*b0j + TMath::Power(TMath::Power(A,1.0/3.0)*delta*b0j,2.0);
  double denominator = delta*delta*delta*delta;
  // std::cout << numerator << " / " << denominator << std::endl;
  return TMath::Power(A,2.0/3.0)*numerator/denominator;
  //return numerator/denominator;
}




void scale_error(TH1F *h[3][2], double entries, double xsec, double scale, double BR)
{
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < 2 ; j++){
      TAxis *xaxis = h[i][j]->GetXaxis();
      int xbins = xaxis->GetNbins();
      for(int k = 1 ; k<=xbins; k++)
	{
	  double nevents_bin = h[i][j]->GetBinError(k);
	  if(nevents_bin==0)continue;
	  double err_nevents_bin = sqrt(nevents_bin);
	  double average_fill_bin = h[i][j]->GetBinContent(k) / nevents_bin;
	  double y = nevents_bin * average_fill_bin * xsec/entries/scale;
	  double y_err = err_nevents_bin * average_fill_bin * xsec/entries/scale;
	  h[i][j]->SetBinContent(k,y);      
	  h[i][j]->SetBinError(k,y_err);
	}  
    }
  }
}



void scale_error_2(TH1F *h[8][2], double entries, double xsec, double scale, double BR)
{
  for(int i = 0 ; i < 8 ; i++){
    for(int j = 0 ; j < 2 ; j++){
      TAxis *xaxis = h[i][j]->GetXaxis();
      int xbins = xaxis->GetNbins();
      for(int k = 1 ; k<=xbins; k++)
	{
	  double nevents_bin = h[i][j]->GetBinError(k);
	  if(nevents_bin==0)continue;
	  double err_nevents_bin = sqrt(nevents_bin);
	  double average_fill_bin = h[i][j]->GetBinContent(k) / nevents_bin;
	  double y = nevents_bin * average_fill_bin * xsec/entries/scale;
	  double y_err = err_nevents_bin * average_fill_bin * xsec/entries/scale;
	  h[i][j]->SetBinContent(k,y);      
	  h[i][j]->SetBinError(k,y_err);
	}  
    }
  }
}
void dQ2(TH1F *h[3][2])
{
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < 2 ; j++){
      TAxis *xaxis = h[i][j]->GetXaxis();
      int xbins = xaxis->GetNbins();
      for(int k = 1 ; k<=xbins; k++)
	{
	  double val = h[i][j]->GetBinContent(k);
	  double err = h[i][j]->GetBinError(k);
	  h[i][j]->SetBinContent(k,val/h[i][j]->GetBinWidth(k));
	  h[i][j]->SetBinError(k,err/h[i][j]->GetBinWidth(k));
	}
    }
  }
}
void dxpom(TH1F *h[3][2])
{
  for(int i = 0 ; i < 3 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  TAxis *xaxis = h[i][j]->GetXaxis();
	  int xbins = xaxis->GetNbins();
	  for(int k = 1 ; k<= xbins ; k++)
	    {
	      double dxpom = 0.0;
	      if(j==0) dxpom = 0.0007;
	      else if(j==1) dxpom = 0.007;
	      double val = h[i][j]->GetBinContent(k);
	      double err = h[i][j]->GetBinError(k);
	      h[i][j]->SetBinContent(k,val/dxpom);
	      h[i][j]->SetBinError(k,err/dxpom); 
	    }
	  
	}
    }
}
void dt(TH1F *h[3][2], double delta_t)
{
  for(int i = 0 ; i < 3 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  TAxis *xaxis = h[i][j]->GetXaxis();
	  int xbins = xaxis->GetNbins();
	  for(int k = 1 ; k<= xbins ; k++)
	    {
	      double val = h[i][j]->GetBinContent(k);
	      double err = h[i][j]->GetBinError(k);
	      h[i][j]->SetBinContent(k,val/delta_t);
	      h[i][j]->SetBinError(k,err/delta_t); 
	    }
	  
	}
    }
}





void dQ2_2(TH1F *h[8][2])
{
  for(int i = 0 ; i < 8 ; i++){
    for(int j = 0 ; j < 2 ; j++){
      TAxis *xaxis = h[i][j]->GetXaxis();
      int xbins = xaxis->GetNbins();
      for(int k = 1 ; k<=xbins; k++)
	{
	  double val = h[i][j]->GetBinContent(k);
	  double err = h[i][j]->GetBinError(k);
	  h[i][j]->SetBinContent(k,val/h[i][j]->GetBinWidth(k));
	  h[i][j]->SetBinError(k,err/h[i][j]->GetBinWidth(k));
	}
    }
  }
}
void dxpom_2(TH1F *h[8][2])
{
 for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  TAxis *xaxis = h[i][j]->GetXaxis();
	  int xbins = xaxis->GetNbins();
	  for(int k = 1 ; k<= xbins ; k++)
	    {
	      double val = h[i][j]->GetBinContent(k);
	      double err = h[i][j]->GetBinError(k);
	      h[i][j]->SetBinContent(k,val/0.004);
	      h[i][j]->SetBinError(k,err/0.004); 
	    }
	  
	}
    }
}
void dt_2(TH1F *h[8][2], double delta_t)
{
for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 2 ; j++)
	{
	  TAxis *xaxis = h[i][j]->GetXaxis();
	  int xbins = xaxis->GetNbins();
	  for(int k = 1 ; k<= xbins ; k++)
	    {
	      double val = h[i][j]->GetBinContent(k);
	      double err = h[i][j]->GetBinError(k);
	      h[i][j]->SetBinContent(k,val/delta_t);
	      h[i][j]->SetBinError(k,err/delta_t); 
	    }
	  
	}
    }
}



double smear_E(double E, double eta, bool isElectron)
{
  if(isElectron)
    {
      if(eta>=-4.5&&eta<=-2.)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.01*E,2.)+pow(0.01,2.)*E)); 
	}
      else if(eta>=-2.&&eta<=-1.)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.02*E,2.)+pow(0.08,2.)*E));
	}
      else if(eta>=-1&&eta<=4.5)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.02*E,2.)+pow(0.12,2.)*E));
	}
      else
	return -1;
    }
  else
    {
      if(eta>=-3.5&&eta<=-1.)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.06*E,2.)+pow(0.45,2.)*E));
	}
      else if(eta>=-1.&&eta<=1.)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.07*E,2.)+pow(0.85,2.)*E));
	}
      else if(eta>=1.&&eta<=3.5)
	{
	  return gRandom->Gaus(E,sqrt(pow(0.06*E,2.)+pow(0.45,2.)*E));
	}
      else
	return -1;
    }
}
