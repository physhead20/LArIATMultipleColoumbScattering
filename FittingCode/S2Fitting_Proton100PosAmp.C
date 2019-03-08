//-------------------------------------------------------------------------
//--------------------------- Fitting Function ----------------------------
//-------------------------------------------------------------------------
Double_t fitf(Double_t *x,Double_t *par) 
{

float momInGeV = x[0] / 1000.;

float S2 = momInGeV*0.006084 + 12.09;

float X0 = 0.14; 
float l = 0.14;
float mass = 0.9383; //(For proton)

Double_t arg = 0;

float velocity_numerator   = (momInGeV/mass)*(momInGeV/mass);
float velocity_denomenator = 1 + ( ((momInGeV/mass)*(momInGeV/mass)) );
float velocity = velocity_numerator/velocity_denomenator;
   
float beta = velocity;

float Term1 = (par[0]) / (momInGeV * beta);
float Term2 = sqrt(l/X0);
float Term3 = Term2*(1 + 1 *log(l/X0));
Term3 = 1;

if (par[1]!=0)  
//arg = Term1*Term2*Term3 ;
arg = Term1*Term3;

Double_t fitval = arg;

return fitval;
}

// ###########################################
// ### Main Function for Fitting and Plots ###
// ###########################################
void S2Fitting_Proton100PosAmp() 
{
// *** Mass of the particle ****
float mass = 0.9383; //(For proton in GeV)

// *** Radiation Length ****
float X0 = 0.14; 
// *** Length of the track segment ****
float l = 0.14;

// **** Variables used in the TGraph ****
float Sigma[5000], P[5000];
int nHighlandPoints = 0;
float epsilon = 0.038;
// ######################################
// ### Nominal Highland Formula plot ####
// ######################################

for( int momentum = 100; momentum < 4000; momentum++)
   {
   
   // *** To make calculations easier, convert momentum to GeV ****
   float momInGeV = momentum / 1000.;
   
   // *** Calculate the velocity term ****
   float velocity_numerator   = (momInGeV/mass)*(momInGeV/mass);
   float velocity_denomenator = 1 + ( ((momInGeV/mass)*(momInGeV/mass)) );
   float velocity = velocity_numerator/velocity_denomenator;
   
   float beta = velocity;
   float term = momInGeV*momInGeV;
   
   float S2 = 13.6;
   
   float Term1 = (S2) / (momentum * beta);
   float Term2 = sqrt(l/X0);
   float Term3 = Term2*(1 + (epsilon *log(l/X0)));
   
   float TestSigma = Term1*Term3;
   
   Sigma[nHighlandPoints] = TestSigma * 1000;
   P[nHighlandPoints] = momentum;
   nHighlandPoints++;
   
   }//<---End momentum plot
TGraph *gHighLandNominal = new TGraph(nHighlandPoints, P, Sigma);
gHighLandNominal->SetLineColor(kBlue);
gHighLandNominal->SetLineWidth(2);
gHighLandNominal->SetFillColor(kWhite);
gHighLandNominal->SetTitle("#pi-MCS; Momentum (MeV); #sigma_{MCS}[milliradians]");

// ############################################
// ### Open a file and get a histogram ########
// ############################################
TFile *f = new TFile("../ROOTFILES/Highland_nTuple/100PosAmp_ProtonData_MomVsTheta3d2_Highland.root");
TH1D *hpx = (TH1D*)f->Get("HighlandFormula");
hpx->Scale(1000);
hpx->SetMarkerColor(kBlack);
hpx->SetLineColor(kBlack);
hpx->SetLineWidth(2);
hpx->SetMarkerStyle(20);


// Create a TF1 object using the function defined above.
// The last three parameters specify the number of parameters
// for the function.
TF1 *func = new TF1("fit",fitf,-1,1,1);

// give the parameters meaningful names
func->SetParNames ("S2");
//func->SetParNames ("Epsilon");
func->SetLineColor(kBlack);
func->SetLineWidth(2);
func->SetLineStyle(2);

hpx->Fit("fit","","",400, 1100);

// ##########################################
// ####		Making the plot 	  ###
// ##########################################
TCanvas *c01 = new TCanvas("c01","S2 Fit");
gHighLandNominal->Draw("");
hpx->Draw("same");
func->Draw("same");

gStyle->SetOptStat(0);
gStyle->SetOptFit(0211);

TLegend *leg = new TLegend();
leg = new TLegend(0.58,0.65,0.88,0.88);
leg->SetTextSize(0.04);
leg->SetTextAlign(12);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);
//leg->SetHeader("LArIAT Run-II");
leg->AddEntry(hpx,"Proton Data");
leg->AddEntry(func,"Highland Fit: S2"); 
leg->AddEntry(gHighLandNominal,"Nominal Highland Formula: S2 = 13.6 MeV");
leg->Draw();


}//<----End Main Function



