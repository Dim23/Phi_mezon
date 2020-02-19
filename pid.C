
float FitFunc(float pt,float a,float b,float c,float d){
return a+b*pt+c/pt+d*log(pt);
}

static float IsPionW(float m2tof, float pt){
  float ispion=-9999.0;
  float sigma=1.0;
  float m2exp=0.01948816;
  sigma=FitFunc(pt,-0.0300609,0.0435236,-0.000149448,-0.0190025);
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}

static float IsPionE(float m2tof, float pt){
  float ispion=-9999.0;
  float sigma=1.0;
  float m2exp=0.01948816;
  sigma=FitFunc(pt,-0.0453771,0.0798237,-0.0162309,-0.0615215);
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}

static float IsKaonW(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.24371698;
  sigma=FitFunc(pt,-0.0372014,0.0823999,-0.0225657,-0.0770799);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

static float IsKaonE(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.24371698;
  sigma=FitFunc(pt,-0.0733706,0.148784,-0.0514034,-0.164131);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

static float IsProtonW(float m2tof, float pt){
  float isproton=-9999.0;
  float sigma=1.0;
  float m2exp=0.880;
  sigma=FitFunc(pt,-0.0640018,-0.00796045,0.130155,0.148396);
  isproton=(m2tof-m2exp)/sigma;
  return  isproton;
}

static float IsProtonE(float m2tof, float pt){
  float isproton=-9999.0;
  float sigma=1.0;
  float m2exp=0.880;
  sigma=FitFunc(pt,-0.0665157,0.119628,0.00906812,-0.114111);
  isproton=(m2tof-m2exp)/sigma;
  return  isproton;
}


TF1 *fitf=new TF1("funcFit","[0]+[1]*x+[2]/x+[3]*log(x)",0.55,2.1);

void NewPid(const char *file_adress){
float aPiW,bPiW,cPiW,dPiW;
float aPiE,bPiE,cPiE,dPiE;
float aKW,bKW,cKW,dKW;
float aKE,bKE,cKE,dKE;
float aPrW,bPrW,cPrW,dPrW;
float aPrE,bPrE,cPrE,dPrE;
TFile *file = new TFile(file_adress);
TGraphErrors *grsig1W = (TGraphErrors*)file->Get("gr_sigW");fitf->SetParameters(-0.03422,0.0589257,-0.0121301,-0.0477698);
grsig1W->Fit(fitf,"RM");
aPiW=fitf->GetParameter(0);bPiW=fitf->GetParameter(1);cPiW=fitf->GetParameter(2);dPiW=fitf->GetParameter(3);
std::cout <<"grsig1W parametrs are "<<aPiW<<","<<bPiW<<","<<cPiW<<","<<dPiW<<endl;

TGraphErrors *grsig1E = (TGraphErrors*)file->Get("gr_sigE");
grsig1E->Fit(fitf,"RM");
aPiE=fitf->GetParameter(0);bPiE=fitf->GetParameter(1);cPiE=fitf->GetParameter(2);dPiE=fitf->GetParameter(3);
std::cout <<"grsig1E parametrs are "<<aPiE<<","<<bPiE<<","<<cPiE<<","<<dPiE<<endl;

TGraphErrors *grsig2W = (TGraphErrors*)file->Get("gr_sig2W");
grsig2W->Fit(fitf,"RM");
aKW=fitf->GetParameter(0);bKW=fitf->GetParameter(1);cKW=fitf->GetParameter(2);dKW=fitf->GetParameter(3);
std::cout <<"grsig2W parametrs are "<<aKW<<","<<bKW<<","<<cKW<<","<<dKW<<endl;

TGraphErrors *grsig2E = (TGraphErrors*)file->Get("gr_sig2E");
grsig2E->Fit(fitf,"RM");
aKE=fitf->GetParameter(0);bKE=fitf->GetParameter(1);cKE=fitf->GetParameter(2);dKE=fitf->GetParameter(3);
std::cout <<"grsig2E parametrs are "<<aKE<<","<<bKE<<","<<cKE<<","<<dKE<<endl;

TGraphErrors *grsig3W = (TGraphErrors*)file->Get("gr_sig3W");fitf->SetParameters(-0.0641154,0.0771587,0.0475797,-0.0230427);
grsig3W->Fit(fitf,"RM");
aPrW=fitf->GetParameter(0);bPrW=fitf->GetParameter(1);cPrW=fitf->GetParameter(2);dPrW=fitf->GetParameter(3);
std::cout <<"grsig3W parametrs are "<<aPrW<<","<<bPrW<<","<<cPrW<<","<<dPrW<<endl;

TGraphErrors *grsig3E= (TGraphErrors*)file->Get("gr_sig3E");
grsig3E->Fit(fitf,"RM");
aPrE=fitf->GetParameter(0);bPrE=fitf->GetParameter(1);cPrE=fitf->GetParameter(2);dPrE=fitf->GetParameter(3);
std::cout <<"grsig3E parametrs are "<<aPrE<<","<<bPrE<<","<<cPrE<<","<<dPrE<<endl;

TCanvas *MyWFit = new TCanvas ("CanSigma","Test canvasqq",1);
  MyWFit -> Divide(3,2);
  MyWFit -> cd(1);grsig1W ->Draw();
  MyWFit -> cd(4);grsig1E -> Draw();
  MyWFit -> cd(2);grsig2W ->Draw();
  MyWFit -> cd(5);grsig2E -> Draw();
  MyWFit -> cd(3);grsig3W ->Draw();
  MyWFit -> cd(6);grsig3E -> Draw();

TFile *d_outfile = new TFile(file_adress,"update");
d_outfile->cd();
MyWFit->Write();
}

