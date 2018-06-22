#include "Helpers.hh"

#include <TApplication.h>
#include <TVirtualFitter.h>

double medianLightYield(TH1D* h)
{
  return TMath::Median(h->GetNbinsX(), &(h->GetArray()[1])); //full range
  //return TMath::Median(56, &(h->GetArray()[261])); //shortRange
}

double meanLightYield(TH1D* h)
{
  double mean=0;
  int n = h->GetNbinsX();
  double* channelLightYield = &(h->GetArray()[1]);

  // full range
  for (int i=0; i<n; ++i)
    if ((((i+1)%64)==0 || (i%64)==0))
      ; else mean+=channelLightYield[i];
  mean /= double(n-16);

  /*  // short range
  int count = 0;
  for (int i=0; i<n; ++i)
    if ((i>260) && (i<315)) {
      mean+=channelLightYield[i];
      count++;
    }
  mean /= double(count-1); */

  return mean;
}

double meanLightYieldError(TH1D* h)
{
  double mean = meanLightYield(h);
  int n = h->GetNbinsX();
  double* channelLightYield = &(h->GetArray()[1]);
  double variance = 0;

  // full range
  for (int i=0; i<n; ++i)
    if (((i+1)%64)==0 || (i%64)==0)
      ; else variance+=TMath::Power(mean-channelLightYield[i], 2);
  variance /= double(n-17);
  double rms= TMath::Sqrt(variance)/TMath::Sqrt(n-16);

  /* // short range
  int count = 0;
  for (int i=0; i<n; ++i)
    if ((i>260) && (i<315)) {
      variance+=TMath::Power(mean-channelLightYield[i], 2);
      count++;
    }
  variance /= double(count-1);
  double rms= TMath::Sqrt(variance)/TMath::Sqrt(double(count-1)); */

  return rms;
}

double medianLightYieldError(TH1D* h)
{
  return meanLightYieldError(h)*1.2533;
}

double MAD(TH1D* h)
{
  using namespace TMath;
  Double_t median = medianLightYield(h);
  Double_t* tempArray = new Double_t[h->GetNbinsX()];
  for (int i=0; i<h->GetNbinsX(); ++i) {
    tempArray[i] = Abs( h->GetArray()[i+1] -median);
  }
  return Median(h->GetNbinsX(), tempArray);
}

TF1* doubleExpFit(TGraphErrors* g)
{
  TF1* expLongComp = new TF1("expLongComp", "[0]*exp(-x/[1])", 100., 230.);
  expLongComp->SetParameters( 25., 700.);
  g->Fit("expLongComp", "R");

  TF1* expFit = new TF1("expFit", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))", 0., 230.);
  expFit->SetParameters( 25., 50., 0.7, 700.);
  expFit->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}");
  expFit->SetParLimits(1, 1., 1000.);
  expFit->SetParLimits(3, 1., 1000.);
  expFit->SetParLimits(2, 0, .99);
  g->Fit("expFit", "ME");
  expFit->Draw("SAME");
  return expFit;
}

TF1* doubleExpFit_Mirror(TGraphErrors* g)
{
  TF1* expFit = new TF1("expFit", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(242.4*2-x)/[1])+[2]*exp(-(242.4*2-x)/[3]))", 0., 200.);
  expFit->SetParameters( 20., 50., 0.7, 400., 0.78);
  expFit->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}", "r");
  expFit->SetParLimits(1, 1., 800.);
  expFit->SetParLimits(3, 1., 800.);
  expFit->SetParLimits(2, .6, .9);
  expFit->SetParLimits(4, .6, .98);
  g->Fit("expFit", "ME");
  expFit->Draw("SAME");
  return expFit;
}

TF1* doubleExpFit_Fresnel(TGraphErrors* g)
{
  TF1* expFit = new TF1("expFit", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(242.35*2-x)/[1])+[2]*exp(-(242.35*2-x)/[3]))", 0., 200.);
  expFit->SetParameters( 20., 50., 0.7, 400., 0.089);
  expFit->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}", "r");
  expFit->SetParLimits(1, 1., 800.);
  expFit->SetParLimits(3, 1., 800.);
  expFit->SetParLimits(2, .6, .9);
  expFit->FixParameter(4, 0.089);
  g->Fit("expFit", "MEQ");
  expFit->SetLineWidth(0.2);
  //expFit->Draw("SAME");

  TGraphErrors* grint = new TGraphErrors(250.);
  for (int i=0; i<250.; ++i)
    grint->SetPoint(i,i,0);
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(grint,0.68);
  grint->Draw("SAME Z");

  {
    double y1= grint->Eval(30.);
    double y2= grint->Eval(0.);
    double e1= grint->GetErrorY(30);
    double e2= grint->GetErrorY(0);
    //qDebug() << y1 << y2 << e1 << e2;
    qDebug() << " Expected loss on 30cm = " << y1/y2 << "pm  " << TMath::Sqrt((1/y1/y1)*e1*e1 + (y1/y2/y2)*(y1/y2/y2)*e2*e2);
  }

  gStyle->SetOptFit(0);
  TLegend* l = new TLegend(.55,.38,.89,.89);
  l->AddEntry(expFit, TString::Format("#chi^{2}/ndf  (%.2f/%i)", expFit->GetChisquare(), expFit->GetNDF()), "l");
  l->AddEntry((TObject*)0, TString::Format("I_{0}  (%.2f #pm %.2f)pixel", expFit->GetParameter(0), expFit->GetParError(0)), "");
  l->AddEntry((TObject*)0, TString::Format("f_{long}  (%.2f #pm %.2f)", expFit->GetParameter(2), expFit->GetParError(2)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{long}  (%.1f #pm %.1f)cm", expFit->GetParameter(3), expFit->GetParError(3)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{short}  (%.1f #pm %.1f)cm", expFit->GetParameter(1), expFit->GetParError(1)), "");
  l->AddEntry((TObject*)0, TString::Format("r_{Fresnel} [fix]  %.3f", expFit->GetParameter(4)), "");
  l->Draw("SAME");

  return expFit;
}

TF1* doubleExpFit_Fresnel_fixLong(TGraphErrors* g)
{
  TF1* preFit = new TF1("preFit", "[0]*exp(-x/[1])+[0]*[2]*exp(-(242.35-x)/[1])", 100., 200.);
  preFit->SetParameters(20., 400., 0.089);
  preFit->FixParameter(2, 0.089);
  preFit->SetLineColor(kGreen);
  g->Fit("preFit", "R");

  TF1* expFit = new TF1("expFit", "[0]*((1-[2])*exp(-x/[1])+[2]*exp(-x/[3]))+ [0]*[4]*((1-[2])*exp(-(242.35*2-x)/[1])+[2]*exp(-(242.35*2-x)/[3]))", 0., 200.);
  expFit->SetParameters( 30., 50., 0.75, 300., 0.089);
  expFit->SetParNames( "I_{0}", "#Lambda_{s}", "f_{l}", "#Lambda_{l}", "r");
  expFit->SetParLimits(1, 1., 800.);
  expFit->SetParLimits(2, .6, .9);
  expFit->FixParameter(3, 305.);
  expFit->FixParameter(4, 0.089);
  g->Fit("expFit", "ME");
  expFit->Draw("SAME");

  gStyle->SetOptFit(0);
  TLegend* l = new TLegend(.55,.38,.89,.89);
  l->AddEntry(expFit, TString::Format("#chi^{2}/ndf  (%.2f/%i)", expFit->GetChisquare(), expFit->GetNDF()), "l");
  l->AddEntry((TObject*)0, TString::Format("I_{0}  (%.2f #pm %.2f)pixel", expFit->GetParameter(0), expFit->GetParError(0)), "");
  l->AddEntry((TObject*)0, TString::Format("f_{long}  (%.2f #pm %.2f)", expFit->GetParameter(2), expFit->GetParError(2)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{long} [fix]  %.1f cm", expFit->GetParameter(3)), "");
  l->AddEntry((TObject*)0, TString::Format("#lambda_{short}  (%.1f #pm %.1f)cm", expFit->GetParameter(1), expFit->GetParError(1)), "");
  l->AddEntry((TObject*)0, TString::Format("r_{Fresnel} [fix]  %.3f", expFit->GetParameter(4)), "");
  l->Draw("SAME");

  return expFit;
}

TF1* singleExpFit_Fresnel(TGraphErrors* g)
{
  TF1* expFit = new TF1("expFit", "[0]*exp(-x/[1])+[0]*[2]*exp(-(242.35-x)/[1])", 0., 230.);
  expFit->SetParameters( 25., 300., 0.089);
  expFit->SetParNames( "I_{0}", "#Lambda_{tech}", "r");
  expFit->SetParLimits(1, 1., 1000.);
  expFit->FixParameter(2, 0.089);
  g->Fit("expFit", "ME");
  expFit->Draw("SAME");
  return expFit;
}

int main(int argc, char** argv) {
  TApplication application("analysis", &argc, argv);
  //TH1::AddDirectory(false);
  Helpers::setRootStyle();

  /* extracting median, error on median, position from files */

  QString filepath = ("/home/iwanicki/sw/opticalCon/data/4TSAACFIM00387/noMirror/"); // TODO GET DATA AND SET PATHES
  QStringList files = QDir(filepath).entryList(QStringList() << "*.root");
  const int nFiles = files.count();
  qDebug() << endl << "using" << nFiles << " files.. ";

  Double_t* mean = new Double_t[nFiles];
  Double_t* median = new Double_t[nFiles];
  Double_t* meanError = new Double_t[nFiles];
  Double_t* medianError = new Double_t[nFiles];
  Double_t* position = new Double_t[nFiles];

  int count = 0;
  foreach (QString file, files)
  {
    qDebug() << "\n  ..processing "<< file;

    QString fullPath = filepath + file;
    auto hist = Helpers::extractFromFile<TH1D>(fullPath, "cluster signal value 0 HistogramDisplay").front();
    mean[count] = meanLightYield(hist);
    median[count] = medianLightYield(hist);
    position[count] = file.left(3).toDouble();
    meanError[count] = meanLightYieldError(hist);
    medianError[count] = medianLightYieldError(hist);

    count++;
  }

  /* distance to SiPM and error calculation */

  Double_t* positionError = new Double_t[nFiles];

  Double_t posError = TMath::Sqrt((3*0.1*0.1 + 0.01*0.01 + 0.02*0.02)/12.); // divide by sqrt(12)?
  for (int i=0; i<nFiles; ++i) {
    position[i] += (0.3 + 8.7 - 6.5 + 0.4);
    positionError[i] = posError;
  }

  /* serial measurements implementation for reference */
  /*
  Double_t* serialPos = new Double_t[4];
  serialPos[0]=218.6; serialPos[1]=138.45; serialPos[2]=77.; serialPos[3]=16.25;

  files = QDir(filepath).entryList(QStringList() << "s2*.root");
  int nSerial = files.count();
  qDebug() << endl << "using" << nSerial << "serial QA files.. ";

  Double_t* newSerialMean = new Double_t[4];
  Double_t* newSerialMeanError = new Double_t[4];
  Double_t* newSerialMedian = new Double_t[4];
  Double_t* newSerialMedianError = new Double_t[4];

  count = 0;
  foreach (QString file, files)
  {
    qDebug() << "\n  ..processing "<< file;

    QString fullPath = filepath + file;
    auto hist = Helpers::extractFromFile<TH1D>(fullPath, "cluster signal value 0 HistogramDisplay").front();
    newSerialMean[count] = meanLightYield(hist);
    newSerialMedian[count] = medianLightYield(hist);
    newSerialMeanError[count] = meanLightYieldError(hist);
    newSerialMedianError[count] = medianLightYieldError(hist);

    count++;
  }
  */

  /*
  filepath = ("/home/dami/master/opticalCon/data/4TSAACFIM00387/serial_old/noMirror/");
  files = QDir(filepath).entryList(QStringList() << "*.root");
  nSerial = files.count();
  qDebug() << endl << "using" << nSerial << "serial(old) QA files.. ";

  Double_t* oldSerialMean = new Double_t[4];
  Double_t* oldSerialMeanError = new Double_t[4];
  Double_t* oldSerialMedian = new Double_t[4];
  Double_t* oldSerialMedianError = new Double_t[4];

  count = 0;
  foreach (QString file, files)
  {
    qDebug() << "\n  ..processing "<< file;

    QString fullPath = filepath + file;
    auto hist = Helpers::extractFromFile<TH1D>(fullPath, "cluster signal value 0 HistogramDisplay").front();
    oldSerialMean[count] = meanLightYield(hist);
    oldSerialMedian[count] = medianLightYield(hist);
    oldSerialMeanError[count] = meanLightYieldError(hist);
    oldSerialMedianError[count] = medianLightYieldError(hist);
    qDebug() << "mean " << oldSerialMean[count]  << "; median " << oldSerialMedian[count] << "; mean error " << oldSerialMeanError[count] << "; median error " << oldSerialMedianError[count];

    count++;
  }
  */

  /* draw and fit */

  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->cd();
  TGraphErrors* g1 = new TGraphErrors(nFiles, position, mean, positionError, meanError);
  //TGraphErrors* g1old = new TGraphErrors(4, serialPos, oldSerialMean, positionError, oldSerialMeanError);
  //TGraphErrors* g1new = new TGraphErrors(4, serialPos, newSerialMean, positionError, newSerialMeanError);
  g1->SetTitle("Mean light yield");
  g1->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g1->GetXaxis()->SetTitleOffset(1.0);
  g1->GetXaxis()->SetTitleSize(35);
  g1->GetXaxis()->SetLabelSize(30);
  g1->GetYaxis()->SetTitleOffset(0.9);
  g1->GetYaxis()->SetTitleSize(35);
  g1->GetYaxis()->SetLabelSize(30);
  g1->GetYaxis()->SetTitle("mean channel signal value [pixel]");
  g1->GetYaxis()->SetRangeUser(10.,26.);
  g1->GetXaxis()->SetLimits(0.,230.);
  g1->Draw("ape");
  /*
  g1old->SetMarkerStyle(22);
  g1old->SetMarkerSize(1.5);
  g1old->SetMarkerColor(kBlue);
  g1old->Draw("same p");
  g1new->SetMarkerStyle(22);
  g1new->SetMarkerSize(1.5);
  g1new->SetMarkerColor(kGreen);
  g1new->Draw("same p");
  */

  TLegend* l1 = new TLegend(.12,.15,.4,.33);
  l1->AddEntry(g1 ,"moveable PMT 2018-feb", "pe");
  //l1->AddEntry(g1new ,"Serial PMTs 2017-nov", "p");
  //l1->AddEntry(g1old ,"Serial PMTs 2017-jan", "p");
  l1->Draw("SAME");

  //TF1* expFitMean = doubleExpFit(g1);
  //TF1* expFitMean = doubleExpFit_Mirror(g1);
  TF1* expFitMean = doubleExpFit_Fresnel(g1);
  //TF1* expFitMean = singleExpFit_Fresnel(g1);
  //TF1* expFitMean = doubleExpFit_Fresnel_fixLong(g1);

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->cd();
  TGraphErrors* g2 = new TGraphErrors(nFiles, position, median, positionError, medianError);
  //TGraphErrors* g2old = new TGraphErrors(4, serialPos, oldSerialMedian, positionError, oldSerialMedianError);
  //TGraphErrors* g2new = new TGraphErrors(4, serialPos, newSerialMedian, positionError, newSerialMedianError);
  g2->SetTitle("Median light yield");
  g2->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g2->GetXaxis()->SetTitleOffset(1.0);
  g2->GetXaxis()->SetTitleSize(35);
  g2->GetXaxis()->SetLabelSize(30);
  g2->GetYaxis()->SetTitleOffset(0.9);
  g2->GetYaxis()->SetTitleSize(35);
  g2->GetYaxis()->SetLabelSize(30);
  g2->GetYaxis()->SetTitle("median channel signal value [pixel]");
  g2->GetYaxis()->SetRangeUser(10., 26.);
  g2->GetXaxis()->SetLimits(0., 230.);
  g2->Draw("ape");
  /*
  g2old->SetMarkerStyle(22);
  g2old->SetMarkerColor(kBlue);
  g2old->SetMarkerSize(1.5);
  g2old->Draw("same p");
  g2new->SetMarkerStyle(22);
  g2new->SetMarkerColor(kGreen);
  g2new->SetMarkerSize(1.5);
  g2new->Draw("same p");
  */

  TLegend* l2 = new TLegend(.12,.15,.4,.33);
  l2->AddEntry(g1 ,"moveable PMT 2018-feb", "pe");
  //l2->AddEntry(g1new ,"Serial PMTs 2017-nov", "p");
  //l2->AddEntry(g1old ,"Serial PMTs 2017-jan", "p");
  l2->Draw("SAME");

  //TF1* expFitMedian = doubleExpFit(g2);
  //TF1* expFitMedian = doubleExpFit_Mirror(g2);
  TF1* expFitMedian = doubleExpFit_Fresnel(g2);
  //TF1* expFitMedian = singleExpFit_Fresnel(g2);
  //TF1* expFitMedian = doubleExpFit_Fresnel_fixLong(g2);

  /* residuals */
  Double_t* zeros = new Double_t[nFiles];
  for (int i=0; i<nFiles; ++i) {
    zeros[i] = 0;
    mean[i] -= expFitMean->Eval(position[i]);
    median[i] -= expFitMedian->Eval(position[i]);
    double posERR_mean = TMath::Abs( expFitMean->Eval(position[i]-positionError[i]) - expFitMean->Eval(position[i]+positionError[i]))/2.;
    double posERR_median = TMath::Abs( expFitMedian->Eval(position[i]-positionError[i]) - expFitMedian->Eval(position[i]+positionError[i]))/2.;
    meanError[i] = TMath::Sqrt(meanError[i]*meanError[i]+ posERR_mean*posERR_mean);
    medianError[i] = TMath::Sqrt(medianError[i]*medianError[i]+ posERR_median*posERR_median);
  }

  /*
  for (int i=0; i<4; ++i) {
    double posERR_mean = TMath::Abs( expFitMean->Eval(serialPos[i]-positionError[i]) - expFitMean->Eval(serialPos[i]+positionError[i]))/2.;
    double posERR_median = TMath::Abs( expFitMedian->Eval(serialPos[i]-positionError[i]) - expFitMedian->Eval(serialPos[i]+positionError[i]))/2.;

    newSerialMean[i] -= expFitMean->Eval(serialPos[i]);
    newSerialMedian[i] -= expFitMedian->Eval(serialPos[i]);
    newSerialMeanError[i] = TMath::Sqrt(meanError[i]*meanError[i]+ posERR_mean*posERR_mean);
    newSerialMedianError[i] = TMath::Sqrt(medianError[i]*medianError[i]+ posERR_median*posERR_median);

    oldSerialMean[i] -= expFitMean->Eval(serialPos[i]);
    oldSerialMedian[i] -= expFitMedian->Eval(serialPos[i]);
    oldSerialMeanError[i] = TMath::Sqrt(meanError[i]*meanError[i]+ posERR_mean*posERR_mean);
    oldSerialMedianError[i] = TMath::Sqrt(medianError[i]*medianError[i]+ posERR_median*posERR_median);
  }
  */

  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->cd();
  TGraphErrors* g1Res = new TGraphErrors(nFiles, position, mean, zeros, meanError);
  //TGraphErrors* g1oldRes = new TGraphErrors(4, serialPos, oldSerialMean, zeros, oldSerialMeanError);
  //TGraphErrors* g1newRes = new TGraphErrors(4, serialPos, newSerialMean, zeros, newSerialMeanError);
  TF1* zeroLine = new TF1("zeroLine", "0.", 0., 300.);
  zeroLine->SetLineStyle(4);
  g1Res->SetTitle("Residual plot (fit on mean values)");
  g1Res->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g1Res->GetXaxis()->SetTitleOffset(1.0);
  g1Res->GetXaxis()->SetTitleSize(35);
  g1Res->GetXaxis()->SetLabelSize(30);
  g1Res->GetYaxis()->SetTitleOffset(0.9);
  g1Res->GetYaxis()->SetTitleSize(35);
  g1Res->GetYaxis()->SetLabelSize(30);
  g1Res->GetYaxis()->SetTitle("(y - f(x))  [pixel]");
  g1Res->GetXaxis()->SetLimits(10., 230.);
  g1Res->GetYaxis()->SetRangeUser(-0.3,0.3);
  g1Res->Draw("ape");
  /*
  g1oldRes->SetMarkerStyle(22);
  g1oldRes->SetMarkerColor(kBlue);
  g1oldRes->SetMarkerSize(1.5);
  g1oldRes->Draw("same p");
  g1newRes->SetMarkerStyle(22);
  g1newRes->SetMarkerColor(kGreen);
  g1newRes->SetMarkerSize(1.5);
  g1newRes->Draw("same p");
  zeroLine->Draw("SAME");
  */

  TCanvas* c4 = new TCanvas("c4", "c4");
  c4->cd();
  TGraphErrors* g2Res = new TGraphErrors(nFiles, position, median, zeros, medianError);
  //TGraphErrors* g2oldRes = new TGraphErrors(4, serialPos, oldSerialMedian, zeros, oldSerialMedianError);
  //TGraphErrors* g2newRes = new TGraphErrors(4, serialPos, newSerialMedian, zeros, newSerialMedianError);
  g2Res->SetTitle("Residual plot (fit on median values)");
  g2Res->GetXaxis()->SetTitle("distance to SiPM [cm]");
  g2Res->GetXaxis()->SetTitleOffset(1.0);
  g2Res->GetXaxis()->SetTitleSize(35);
  g2Res->GetXaxis()->SetLabelSize(30);
  g2Res->GetYaxis()->SetTitleOffset(0.9);
  g2Res->GetYaxis()->SetTitleSize(35);
  g2Res->GetYaxis()->SetLabelSize(30);
  g2Res->GetYaxis()->SetTitle("(y - f(x))  [pixel]");
  g2Res->GetXaxis()->SetLimits(10., 230.);
  g2Res->GetYaxis()->SetRangeUser(-0.3,0.3);
  g2Res->Draw("ape");
  /*
  g2oldRes->SetMarkerStyle(22);
  g2oldRes->SetMarkerColor(kBlue);
  g2oldRes->SetMarkerSize(1.5);
  g2oldRes->Draw("samep");
  g2newRes->SetMarkerStyle(22);
  g2newRes->SetMarkerColor(kGreen);
  g2newRes->SetMarkerSize(1.5);
  g2newRes->Draw("samep");
  zeroLine->Draw("SAME");
  */

  application.Run();
  return 0;
}
