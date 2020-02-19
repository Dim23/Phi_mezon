#define MTREE_cxx
#include "MTREE.h"


int pn=23;TH1F *hmassW[23];TH1F *hmassE[23];TH1F *FitW[23];TH1F *FitE[23]; TH1F *invW;TH1F *invE;TH1F *invWE;TH1F *invALL;TH1F *normKW;TH1F *normKE;
static const double bin_w[24]={0.3,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8};

double *binP=new double[pn];double *binPer=new double[pn];
double *sigW=new double[pn];double *sigE=new double[pn];double *sigWer=new double[pn];double *sigEer=new double[pn];
double *sig2W=new double[pn];double *sig2E=new double[pn];double *sig2Wer=new double[pn];double *sig2Eer=new double[pn];
double *sig3W=new double[pn];double *sig3E=new double[pn];double *sig3Wer=new double[pn];double *sig3Eer=new double[pn];
float pt,eta,m,beta,dtime,charge_mom;
const float CC=299792458,Mpi=0.019479955,Mk=0.24371698,Mpr=0.880354499,pi=3.141592654;
char strE[20],strW[20];
int NN,NNTREE;
    TH1F *dtimW;
    TH1F *dtimE;
    TH1F *hm2E ;
    TH1F *hm2W  ;
    TH2F *hpidW;
    TH2F *hpidE;
    TH1F *htofW;
    TH1F *htofE;
	TH1F *hphiW;TH1F *htetaW; TH1F *hphiE;TH1F *htetaE;
    TFile *d_outfile;TFile *d_outfileinv;

void MTREE::Book_Hist(){
int NN=0,NNTREE=0;

for(int n=0;n<5;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",300,-0.2,1.2);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",300,-0.2,1.2);}

for(int n=5;n<10;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",150,-0.2,1.2);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",150,-0.2,1.2);}

for(int n=10;n<pn;n++){
sprintf(strE,"%s%d","mass2E",n);
sprintf(strW,"%s%d","mass2W",n);
const char *WestH=(char*)strW;
const char *EastH=(char*)strE;
hmassW[n]=new TH1F(WestH,"mass^{2} for Tof.West",150,-0.2,1.2);
hmassE[n]=new TH1F(EastH,"mass^{2} for Tof.East",150,-0.2,1.2);}

dtimW=new TH1F("dtimeW","tof - t(exp for #pi),ns  TOF West Arm",200,-2,10);
dtimE=new TH1F("dtimeE","tof - t(exp for #pi),ns  TOF East Arm ",200,-2,10);

hm2E  = new TH1F("hm2E","mass^{2} for Tof.East",300,-0.2,1.8);
hm2W  = new TH1F("hm2W","mass^{2} for Tof.West",300,-0.2,1.8);
invW  = new TH1F("invW","inv West",300,0.98,1.28);
invE  = new TH1F("invE","inv East",300,0.98,1.28);

hphiW = new TH1F("phiW","phi West",400,-4,4); htetaW = new TH1F("tetaW","teta for west",400,-4,4);
hphiE = new TH1F("phiE","phi east",400,-4,4); htetaE = new TH1F("tetaE","teta for east",400,-4,4);

invWE  = new TH1F("invWE","one Kaon in West and East",500,0,7);
invALL  = new TH1F("invALL","West plus East",300,0.98,1.28);

normKE  = new TH1F("normKE","NORM K East",50,-4,4);
normKW  = new TH1F("normKW","NORM K West",50,-4,4);

hpidW= new TH2F("hpidW","charge/momentum vs time of flight,  West Arm ",200,10,60,400,-2.2,2.2);
hpidE= new TH2F("hpidE","charge/momentum vs time of flight,  East Arm ",200,10,60,400,-2.2,2.2);

htofW= new TH1F("htofW","tof - t(exp for #pi),ns  TOF West Arm ",200,-2,3);
htofE= new TH1F("htofE","tof - t(exp for #pi),ns  TOF East Arm ",200,-3,3);
}

void MTREE::Loop()
{
//   In a ROOT session, you can do:
//      root> .L hTANA.C
//      root> hTANA t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
for(int n=0;n<pn;n++){ binP[n]=0.5*(bin_w[n]+bin_w[n+1]);}

Loop_imp();

/*if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
        
        Long64_t ientry = LoadTree(jentry);
        if(ientry%100000==0) cout << ientry<<endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

if(fabs(bbcz)<30 && cent>0 && cent<=80)
{
    for (int i=0;i<mh;i++)
    {
        pt=p[i]*sin(the0[i]);charge_mom=charge[i]/p[i];
    
        if (pt>1 && pt<1.5 && pltof[i]>0 && sigtof[i]<3) 
        {   
            m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1), 
            //m=p[i]*sqrt((pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1)),
            beta=p[i]/sqrt(pow(p[i],2)+pow(Mpi,2)),dtime=ttof[i]-pltof[i]*(1e+7)/(beta*CC);

            if (dcarm[i]==0 && etof[i]>0.002)
            {dtimE->Fill(dtime);htofE->Fill(dtime);hm2E->Fill(m);hpidE->Fill(ttof[i],charge_mom);}; 

            if (dcarm[i]==1 && etof[i]>60 && etof[i]<600)
            {dtimW->Fill(dtime);htofW->Fill(dtime);hm2W->Fill(m);hpidW->Fill(ttof[i],charge_mom);};
        }
    }
}}

*/
}

void MTREE::add_file(const char *file) {
  TFile *treefile = TFile::Open(file);
  TTree *tree = (TTree*)treefile->Get("mtree");
  if(tree == 0) {
    cout <<"htree is not found in "<<file<<endl;
    treefile->Close();
    return;
  }
  cout << file <<" is opened"<<endl;
  Init(tree);
  cout <<"one file processed"<<endl;
}

void MTREE::Loop_imp()
{
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast(),nbytes = 0, nb = 0;

for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {NN+=1;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if(ientry==0){NNTREE+=1;}
        if(ientry%100000==0) {cout << ientry<< " all entri is "<< NN <<" the namber of TRee is "<< NNTREE<<endl;}
        nb = fChain->GetEntry(jentry);  nbytes += nb;


if(fabs(bbcz)<30 && cent>20 && cent<60 )
{
    for (int i=0;i<mh;i++)
    {  
        pt=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);NN+=1;

for(int n=0;n<pn;n++){
 
        if ( pt>bin_w[n] && pt<bin_w[n+1] && sigtof[i]<3 && pltof[i]>0) 
        {   
            if (dcarm[i]==0 && etof[i]>0.002)

            {hmassE[n]->Fill(m);}; 

            if (dcarm[i]==1 && etof[i]>60 && etof[i]<600)

            {hmassW[n]->Fill(m);};
        }}
    }
}}
cout << " Histograms for fitting was lopped"<<endl;
}

void MTREE::INV()
{
float m12,mk,pti,ptk,pxi,pxk,pyi,pyk,pzi,pzk,Pt;
float PtMin=0.3,PtMax=3.8,Pz,Py,Px,Ei,Ek;
if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast(),nbytes = 0, nb = 0;

for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {NN+=1;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if(ientry==0){NNTREE+=1;cout << " all entri is "<< NN <<" the naber of TRee is "<< NNTREE<<endl;}
        nb = fChain->GetEntry(jentry);  nbytes += nb;


    if(fabs(bbcz)<30 && cent>20 && cent<60 )
        {
        for (int i=0;i<mh;i++)
            {  
            pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);
pxi=pti*cos(phi0[i]);pyi=pti*sin(phi0[i]);pzi=p[i]*cos(the0[i]);
hphiE->Fill(phi0[i]);htetaE->Fill(the0[i]);

            if (dcarm[i]==0 && etof[i]>0.002  && fabs(IsKaonE(m,pti))<3 && fabs(IsPionE(m,pti))>3 && pti>=PtMin &&pti<PtMax && m>0.12 && m<0.5)

            {normKE->Fill(IsKaonE(m,pti));hm2E->Fill(m);
            for (int k=i;k<mh;k++)
            {
        mk=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);ptk=p[k]*sin(the0[k]);
pxk=ptk*cos(phi0[k]);pyk=ptk*sin(phi0[k]);pzk=p[k]*cos(the0[k]);
Pt=sqrt(pow((pxi+pxk),2)+pow((pyi+pyk),2));Px=pxi+pxk;Py=pyi+pyk;Pz=pzi+pzk;
            if (charge[i]!=charge[k] && dcarm[i]==0 && etof[k]>0.002 && fabs(IsKaonE(mk,ptk))<3 && fabs(IsPionE(m,ptk))>3 && ptk>PtMin &&ptk<PtMax  && mk>0.12 && mk<0.5){

            Ei=sqrt(p[i]*p[i]+Mk);Ek=sqrt(p[k]*p[k]+Mk);
            m12=sqrt(pow((Ei+Ek),2)-(Px*Px+Py*Py+Pz*Pz));
if (Pt>1.2){invE->Fill(m12);invALL->Fill(m12);}}}}; 
        }   }


        if(fabs(bbcz)<30 && cent>20 && cent<60 )
        {
        for (int i=0;i<mh;i++)
            {  
            pti=p[i]*sin(the0[i]);m=pow(p[i],2)*(pow(ttof[i]*(1e-7)*CC/pltof[i],2)-1);pxi=pti*cos(phi0[i]);pyi=pti*sin(phi0[i]);
		hphiW->Fill(phi0[i]);htetaW->Fill(the0[i]);pzi=p[i]*cos(the0[i]);

            if (dcarm[i]==1 && etof[i]>60 && etof[i]<600 && fabs(IsKaonW(m,pti))<2 && fabs(IsPionW(m,pti))>3 && pti>=PtMin && pti<PtMax && m>0.18 && m<0.4)

            {normKW->Fill(IsKaonW(m,pti));hm2W->Fill(m);
            for (int k=i;k<mh;k++)
            {
        mk=pow(p[k],2)*(pow(ttof[k]*(1e-7)*CC/pltof[k],2)-1);ptk=p[k]*sin(the0[k]);pxk=ptk*cos(phi0[k]);pyk=ptk*sin(phi0[k]);
	Pt=sqrt(pow((pxi+pxk),2)+pow((pyi+pyk),2));pzk=p[k]*cos(the0[k]);Px=pxi+pxk;Py=pyi+pyk;Pz=pzi+pzk;

            if (charge[i]!=charge[k] && dcarm[k]==1 && etof[k]>60 && etof[k]<600 && fabs(IsKaonW(mk,ptk))<2 && fabs(IsPionW(m,ptk))>3 && ptk>PtMin &&ptk<PtMax && mk>0.15 && mk<0.4){

            Ei=sqrt(p[i]*p[i]+Mk);Ek=sqrt(p[k]*p[k]+Mk);
            m12=sqrt(pow((Ei+Ek),2)-(Px*Px+Py*Py+Pz*Pz));
if (Pt>1.2){invW->Fill(m12);invALL->Fill(m12);}}}}; 
        }   }

    }
cout << " Histograms for invariant mass was lopped"<<endl;
/*
TCanvas *MyW1 = new TCanvas ("CanSigma1","Test canvasqq1",1);MyW1->Divide(4,1);
MyW1 -> cd(1);invE->GetXaxis()->SetTitle("inv mass,Gev");invE->Draw();
c=invE->GetMaximum();TLine *line1=new TLine();
line1->DrawLine(1.02,0,1.02,c);cout <<gPad->GetUymax()<<endl;
MyW1 -> cd(2);invW->GetXaxis()->SetTitle("inv mass,Gev");invW->Draw();
c=invW->GetMaximum();TLine *line2=new TLine();
line2->DrawLine(1.02,0,1.02,c);cout <<MyW1->cd(2)->GetUymax()<<endl;
MyW1 -> cd(3);invWE->GetXaxis()->SetTitle("inv mass,Gev");invWE->Draw();
MyW1 -> cd(4);invALL->Draw();c=invALL->GetMaximum();TLine *line3=new TLine();
line3->DrawLine(1.02,0,1.02,c);*/
}

void MTREE::INVsave(const char *outfile)  {
TCanvas *MyW1 = new TCanvas ("INVmass","Test canvasqq1",1);MyW1->Divide(2,1);
MyW1 -> cd(1);invE->GetXaxis()->SetTitle("inv mass,Gev");invE->Draw();
TLine *line1=new TLine();
line1->DrawLine(1.02,0,1.02,invE->GetMaximum());
MyW1 -> cd(2);invW->GetXaxis()->SetTitle("inv mass,Gev");invW->Draw();
TLine *line2=new TLine();
d_outfileinv = new TFile(outfile,"recreate");
line2->DrawLine(1.02,0,1.02,invW->GetMaximum());
d_outfileinv->cd();
invE->Write();invW->Write();MyW1->Write();invALL->Write();normKW->Write();normKE->Write();
hm2E->Write();hm2W->Write();hphiW->Write();htetaW->Write();hphiE->Write();htetaE->Write();

d_outfileinv->Close();}             

               

//Функция для фита;
double FitGaus(double *x, double *par) 
{
double g1,g2,g3,pdf;
g1=exp(-pow(x[0]-par[0],2)/(2*par[1]*par[1]))/(par[1]*pow(2*pi,0.5));
g2=exp(-pow(x[0]-par[2],2)/(2*par[3]*par[3]))/(par[3]*pow(2*pi,0.5));
g3=exp(-pow(x[0]-par[4],2)/(2*par[5]*par[5]))/(par[5]*pow(2*pi,0.5));
pdf=par[6]*(par[7]*g1+par[8]*g2+(1-par[7]-par[8])*g3);
    return pdf;
};

//Функция для фита;
double FitGaus2(double *x, double *par) 
{
double g1,g2,g3,pdf;
g1=exp(-pow(x[0]-par[0],2)/(2*par[1]*par[1]))/(par[1]*pow(2*pi,0.5));
pdf=par[2]*g1;
    return pdf;
};

float Func(float pt,float a,float b,float c,float d){
return a+b*pt+c/pt+d*log(pt);
}
//Функция для опраксимации sigma и mean;
double Func2(double x,double a0,double a1,double a2,double a3,double a4) 
{
    return (a0+a1/x+a2/(x*x))*exp(a3*x+a4*x*x);
};

void MTREE::Fit_imp(){
TF1 *fitPi=new TF1("Pi",FitGaus2,-0.1,0.1,3);
TF1 *fitK=new TF1("K",FitGaus2,0.15,0.4,3);
TF1 *fitPr=new TF1("Pr",FitGaus2,0.7,1.1,3);

TF1 *EfitPi=new TF1("EPi",FitGaus2,-0.1,0.1,3);
TF1 *EfitK=new TF1("EK",FitGaus2,0.15,0.4,3);
TF1 *EfitPr=new TF1("EPr",FitGaus2,0.7,1.1,3);

TF1 *f1=new TF1("fc1",FitGaus,-0.2,1.2,9);
TF1 *f2=new TF1("fc2",FitGaus,-0.2,1.2,9);
double conW,meanm2W,meanm3W,sigmW,sigm2W,sigm3W;
double conE,meanm2E,meanm3E,sigmE,sigm2E,sigm3E;
for(int n=0;n<pn;n++){
conW=hmassW[n]->Integral(1, hmassW[n]->GetNbinsX(), "width");conE=hmassE[n]->Integral(1, hmassE[n]->GetNbinsX(), "width");

sigmW=Func(binP[n],-0.03422,0.0589257,-0.0121301,-0.0477698);sigmE=Func(binP[n],-0.0505325,0.0801884,-0.010992,-0.0545074);
sigm2W=Func(binP[n],-0.0534937,0.089851,-0.0133339,-0.0710717);sigm2E=Func(binP[n],-0.0729391,0.0970988,0.000729498,-0.0523558);
sigm3W=Func(binP[n],-0.0641154,0.0771587,0.0475797,-0.0230427);sigm3E=Func(binP[n],-0.100868,0.107671,0.052214,-0.0540008);

//f1->SetParameters(Mpi,sigmW,Mk,sigm2W,Mpr,sigm3W,conW,0.5,0.3);
fitPi->SetParameter(0,Mpi);fitPi->SetParameter(1,sigmW);

fitK->SetParameter(0,Mk);fitK->SetParameter(1,sigm2W);
fitPr->SetParameter(0,Mpr);fitPr->SetParameter(1,sigm3W);
EfitPi->SetParameter(0,Mpi);EfitPi->SetParameter(1,sigmE);

EfitK->SetParameter(0,Mk);EfitK->SetParameter(1,sigm2E);
EfitPr->SetParameter(0,Mpr);EfitPr->SetParameter(1,sigm3E);
//f2->SetParameters(Mpi,sigmE,Mk,sigm2E,Mpr,sigm3E,conE,0.5,0.3);
//f1->FixParameter(0,Mpi);f1->FixParameter(2,Mk);f1->FixParameter(4,Mpr);
//f2->FixParameter(0,Mpi);f2->FixParameter(2,Mk);f2->FixParameter(4,Mpr);

hmassW[n]->Fit(fitPi,"RM");
sigW[n]=fitPi->GetParameter(1);sigWer[n]=fitPi->GetParError(1);
hmassW[n]->Fit(fitPr,"RM");
sig3W[n]=fitPr->GetParameter(1);sig3Wer[n]=fitPr->GetParError(1);
hmassW[n]->Fit(fitK,"RM");
sig2W[n]=fitK->GetParameter(1);sig2Wer[n]=fitK->GetParError(1);

hmassE[n]->Fit(EfitPi,"RM");
sigE[n]=EfitPi->GetParameter(1);sigEer[n]=EfitPi->GetParError(1);
hmassE[n]->Fit(EfitPr,"RM");
sig3E[n]=EfitPr->GetParameter(1);sig3Eer[n]=EfitPr->GetParError(1);
hmassE[n]->Fit(EfitK,"RM");
sig2E[n]=EfitK->GetParameter(1);sig2Eer[n]=EfitK->GetParError(1);

//hmassE[n]->Fit(f2,"R");
//sigE[n]=f2->GetParameter(1);sigEer[n]=f2->GetParError(1);
//sig2E[n]=f2->GetParameter(3);sig2Eer[n]=f2->GetParError(3);
//sig3E[n]=f2->GetParameter(5);sig3Eer[n]=f2->GetParError(5);
binPer[n]=0.5*(bin_w[n+1]-bin_w[n]);
    }
}


void MTREE::ana_end(const char *outfile) {

Fit_imp();
TGraphErrors *grsig1W = new TGraphErrors (pn,binP, sigW, binPer, sigWer);grsig1W->SetName("gr_sigW"); grsig1W->GetYaxis()->SetTitle("<Sigma>");grsig1W->GetXaxis()->SetTitle("Pt,Gev/c");
TGraphErrors *grsig1E = new TGraphErrors (pn,binP, sigE,binPer, sigEer);grsig1E->SetName("gr_sigE"); grsig1E->GetYaxis()->SetTitle("<Sigma>");grsig1E->GetXaxis()->SetTitle("Pt,Gev/c");
grsig1W->SetLineColor(1);grsig1E->SetLineColor(3);grsig1W->SetTitle("sigma for Pion, West");grsig1E->SetTitle("sigma for Pion, East");

TGraphErrors *grsig2W = new TGraphErrors (pn,binP, sig2W,binPer, sig2Wer);grsig2W->SetName("gr_sig2W");
TGraphErrors *grsig3W = new TGraphErrors (pn,binP, sig3W,binPer, sig3Wer);grsig3W->SetName("gr_sig3W"); 
grsig2W->SetTitle("sigma for Kaon, West");grsig3W->SetTitle("sigma for Proton, West");
TGraphErrors *grsig2E = new TGraphErrors (pn,binP, sig2E,binPer, sig2Eer);grsig2E->SetName("gr_sig2E"); 
TGraphErrors *grsig3E = new TGraphErrors (pn,binP, sig3E,binPer, sig3Eer);grsig3E->SetName("gr_sig3E");
grsig2E->SetTitle("sigma for Kaon, East");grsig3E->SetTitle("sigma for Proton, East");
TCanvas *MyW = new TCanvas ("CanSigma","Test canvasqq",1);
MyW -> Divide(3,2);
  MyW -> cd(1);grsig1W ->Draw();
  MyW -> cd(4);grsig1E -> Draw();
  MyW -> cd(2);grsig2W ->Draw();
  MyW -> cd(5);grsig2E -> Draw();
  MyW -> cd(3);grsig3W ->Draw();
  MyW -> cd(6);grsig3E -> Draw();

/*TCanvas *MyC = new TCanvas ("CanData","Test canvas",1);
  MyC -> Divide(3,2);
  MyC -> cd(1);hm2E ->Draw();
  MyC -> cd(4);hm2W -> Draw();
  MyC -> cd(2);hpidE->Draw("colz");
  MyC -> cd(5);hpidW->Draw("colz");
  MyC -> cd(3);htofE -> Draw();
  MyC -> cd(6);htofW -> Draw();*/

TCanvas *MyC2W = new TCanvas ("CanFitW","Test canvasW",1);
  MyC2W -> Divide(3,2);
  MyC2W -> cd(1);hmassW[3] ->Draw();
  MyC2W -> cd(2);hmassW[6] -> Draw();
  MyC2W -> cd(3);hmassW[9]->Draw();
  MyC2W -> cd(4);hmassW[11]->Draw();
  MyC2W -> cd(5);hmassW[12] -> Draw();
  MyC2W -> cd(6);hmassW[13] -> Draw();

TCanvas *MyC2E = new TCanvas ("CanFitE","Test canvasE",1);
  MyC2E -> Divide(3,2);
  MyC2E -> cd(1);hmassE[3] ->Draw();
  MyC2E -> cd(2);hmassE[6] -> Draw();
  MyC2E -> cd(3);hmassE[9]->Draw();
  MyC2E -> cd(4);hmassE[11]->Draw();
  MyC2E -> cd(5);hmassE[12] -> Draw();
  MyC2E -> cd(6);hmassE[13] -> Draw();
d_outfile = new TFile(outfile,"recreate");
d_outfile->cd();
MyW->Write();
//MyC->Write();
MyC2W->Write();MyC2E->Write();
grsig1W->Write();grsig1E->Write();
grsig2W->Write(); grsig3W->Write();
//grmean2W->Write(); grmean3W->Write();
grsig2E->Write(); grsig3E->Write();
//grmean2E->Write(); grmean3E->Write();
//hpidW->Write();hpidE->Write();
//htofW->Write();htofE->Write();
//dtimW->Write();dtimE->Write();
//hm2W->Write();hm2E->Write();
//hbbcz->Write();
for(int n=0;n<pn;n++){
hmassW[n]->Write();
hmassE[n]->Write();}
d_outfile->Close();}

