#ifndef extrapolationUtils_Yield_pp_h
#define extrapolationUtils_Yield_pp_h


// Custom headers
#include "../Macros/fitUtils.h"


double yield_PtRap_pp(const double& pT, const double& rap, const double& sqrtS, const double& mass, const int& charge)
{
  // check input parameters
  if (pT<0) throw runtime_error("[ERROR] Negative pT value for pp yield!");
  if (mass<=0.05) throw runtime_error("[ERROR] Invalid mass value for pp yield!");
  if (sqrtS<=0) throw runtime_error("[ERROR] Invalid beam energy value for pp yield!");
  if (charge==0) throw runtime_error("[ERROR] Can not describe neutral particles for pp yield!");

  // extract the mass dependent parameters
  double C(0), n(0), dNdy(0), beamEneRef(0);
  if (sqrtS>=9) {
    beamEneRef = 13;
    if (mass>0.) {
      // source: Deuteron spectra in pp 13 TeV, https://alice-figure.web.cern.ch/node/11817  ,  https://alice-figure.web.cern.ch/node/11380
      if (charge>0) { C = 0.154183; n = 3.918055; dNdy = 0.00025987494; }
      else          { C = 0.109075; n = 3.008592; dNdy = 0.00028896250; }
    }
    // compute Helium3 and Triton based on model expectation for ratio to Deuteron
    if (mass>2.5) {
      // source: Suppression of light nuclei production in collisions of small systems at the Large Hadron Collider, https://doi.org/10.1016/j.physletb.2019.03.033
      const auto dNdeta = 6.01*TMath::Power(sqrtS/7., 0.206); // using as reference the averaged dN/deta at 7 TeV
      const auto Tk = 80.6 + 83.*TMath::Power(1. + 2.33*dNdeta/67.3,  -1./2.33);
      const auto R = 197.3 * TMath::Power(4.18125*dNdeta, 1./3.) / TMath::Sqrt(938.272*Tk);
      const auto dOverP = 4.E-3/TMath::Power(1. + 1.6*1.6/R/R, 1.5);
      // use two-body coalenscence model, as it extrapolates better Triton at 7 TeV
      const auto He3OverP = 7.1E-6/TMath::Power(1. + 1.15*1.15/R/R, 1.5)/TMath::Power(1. + 1.6*1.6/R/R, 1.5); // Two-body coalescence
      const auto H3OverHe3 = TMath::Power(1. + 1.15*1.15/R/R, 1.5)/TMath::Power(1. + 1.039*1.039/R/R, 1.5); // Two-body coalescence
      //const auto He3OverP = 7.1E-6/TMath::Power(1. + 1.24*1.24/R/R, 3.); // Three-body coalescence
      //const auto H3OverHe3 = TMath::Power(1. + 1.24*1.24/R/R, 3.)/TMath::Power(1. + 1.12*1.12/R/R, 3.); // Three-body coalescence
      // case: Helium3
      dNdy *= He3OverP / dOverP;
      if (std::abs(charge)==1) {
	// case: Triton
	dNdy *= H3OverHe3;
      }
      if (mass>3.0) {
	// source: Extrapolated He4/He3 ratio from light nuclei dN/dy in pp 7 TeV, Phys. Rev. C 97, 024615 (2018)
	dNdy *= 0.000933191617864478;
      }
    }
  }
  else {
    beamEneRef = 7;
    if (mass>0.) {
      // source: Deuteron spectra in pp 7 TeV, Phys. Rev. C 97, 024615 (2018)
      if (charge>0) { C = 0.219976; n = 6.805574; dNdy = 0.00020761550; }
      else          { C = 0.221778; n = 7.540039; dNdy = 0.00019701998; }
    }
    if (mass>2.5) {
      // source: Triton spectra in pp 7 TeV, Phys. Rev. C 97, 024615 (2018)
      if (charge>0) { C = 0.130642; n = 6.245758; dNdy = 2.0272614e-07; }
      else          { C = 0.230485; n = 4.817156; dNdy = 2.0004223e-07; }
    }
    if (std::abs(charge)==2) {
      if (mass>2.5) {
	// source: Helium3 spectra in pp 7 TeV, Phys. Rev. C 97, 024615 (2018)
	if (charge>0) { C = 0.180322; n = 12.242740; dNdy = 2.8067547e-07; }
	else          { C = 0.155575; n = 3.613376;  dNdy = 1.6578813e-07; }
      }
      if (mass>3.0) {
	// source: Extrapolated He4/He3 ratio from light nuclei dN/dy in pp 7 TeV, Phys. Rev. C 97, 024615 (2018)
	dNdy *= 0.000933191617864478;
      }
    }
  }

  // compute pT spectrum
  static TF1* fPt(0);
  if (!fPt) fPt = new TF1("fPtTsallisLevyYield", tsallisLevyYield, 0., 100., 5);
  const auto res_pT = evalTF(*fPt, {1., mass, C, n}, pT);

  // compute y spectrum using data from -0.5 < y < 0.5
  // Assume a Landau energy dependent Gaussian rapidity distribution, based on https://arxiv.org/pdf/1604.02651.pdf
  const auto width_ref = TMath::Sqrt(TMath::Log(beamEneRef/0.000938/2));
  const auto yInt0p5 = std::sqrt(2*TMath::Pi()) * width_ref * TMath::Erf(std::sqrt(2)/width_ref/4); // from dNch/dy vs s
  const auto beamEneScale = TMath::Power(sqrtS, 2*0.103)/TMath::Power(beamEneRef, 2*0.103); // from dNch/deta vs s
  
  const auto width = TMath::Sqrt(TMath::Log(sqrtS/0.000938/2));
  const auto res_y = beamEneScale * dNdy * TMath::Gaus(rap, 0., width) / yInt0p5;

  return res_y * res_pT;
};

double yield_PtRap_pp(const double* x, const double* p)
{
  // input variables
  const auto& t = p[3];
  const auto& pT  = t==0 ? x[0] : (t==1 ? p[4] : x[0]);
  const auto& rap = t==0 ? x[1] : (t==1 ? x[0] : p[4]);
  if (t<0 || t>2) throw runtime_error("[ERROR] Wrong type value for pp yield3D!");
  
  // input parameters
  const auto& sqrtS  = p[0];
  const auto& mass   = p[1];
  const auto& charge = p[2];

  return yield_PtRap_pp(pT, rap, sqrtS, mass, charge);
};

double yield_PtEta_pp(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& sqrtS  = p[0];
  const auto& mass   = p[1];
  const auto& charge = p[2];
  const auto& eta    = p[3];

  const auto rap = computeRapidity(pT, eta, mass);
  const auto mT = TMath::Sqrt(mass*mass + pT*pT);
  const auto arg = mass/mT/TMath::CosH(rap);
  const auto jac = TMath::Sqrt(1. - arg*arg);

  return jac*yield_PtRap_pp(pT, rap, sqrtS, mass, charge);
};

double yield_Eta_pp(const double* x, const double* p)
{
  // input variables
  const auto& eta = x[0];
  
  // input parameters
  const auto& sqrtS  = p[0];
  const auto& mass   = p[1];
  const auto& charge = p[2];
  
  static TF1* f1D(0);
  if (!f1D) f1D = new TF1("f1DYield_Eta", yield_PtEta_pp, 0., 100., 4);
  
  const auto res = evalTF(*f1D, {sqrtS, mass, charge, eta}, 0., 100.);

  return res;
};

double spectra_PtRap_pp(const double * x, const double* p)
{
  const auto& pT = p[3]==1 ? p[4] : x[0];
  return yield_PtRap_pp(x, p) / (2*TMath::Pi()*pT);
};

double dist_Inc1D_pp(TF1& f1D, const UChar_t& id, const double* x, const double* p)
{
  // input parameters
  const auto& scale   = p[0];
  const auto& type    = p[1]; // function of:   (0) pT,     (1) rapidity,     (2) pT (integrate |y|)
  const auto& sqrtS   = p[2];
  const auto& mass    = p[3];
  const auto& charge  = p[4];
  const auto& varMin  = p[5];
  const auto& varMax  = p[6];
  
  // check input parameters
  if (type<0 || type>2) throw runtime_error("[ERROR] Wrong type value for pp dist 1D!");
  if (scale<=0) throw runtime_error("[ERROR] Wrong scale value for pp dist 1D!");
  if (varMin>=varMax) throw runtime_error("[ERROR] Wrong var range values for pp dist 1D!");
  
  double res(0);
  if (type<2) {
    res += evalTF(f1D, {sqrtS, mass, charge, type+1, x[0]}, varMin, varMax, id);
  }
  else if (type==2 && varMin==0) {
    res += evalTF(f1D, {sqrtS, mass, charge, 1, x[0]}, -varMax, varMax, id);
  }
  else if (type==2 && varMin>0) {
    res += evalTF(f1D, {sqrtS, mass, charge, 1, x[0]}, -varMax, -varMin, id);
    res += evalTF(f1D, {sqrtS, mass, charge, 1, x[0]}, varMin, varMax, id);
  }

  return scale * res;
};

double yield_Inc1D_pp(const double * x, const double* p)
{
  // input parameters
  const auto& varMax = p[6];
  
  static TF1* f1D(0);
  if (!f1D) f1D = new TF1("f1DYield", yield_PtRap_pp, -varMax, varMax, 5);

  return dist_Inc1D_pp(*f1D, 0, x, p);
};

double yield_Diff1D_pp(const double * x, const double* p)
{
  // input parameters
  const auto& type   = p[1]; // function of:   (0) pT,     (1) rapidity,     (3) pT (integrate |y|)
  const auto& varMin = p[5];
  const auto& varMax = p[6];
  
  auto binWidth = (varMax - varMin);
  if (type==2) binWidth *= 2;

  return yield_Inc1D_pp(x, p) / binWidth;
};

double spectra_Inc1D_pp(const double * x, const double* p)
{
  // input parameters
  const auto& varMax = p[6];
  
  static TF1* f1D(0);
  if (!f1D) f1D = new TF1("f1DSpectra", spectra_PtRap_pp, -varMax, varMax, 5);

  return dist_Inc1D_pp(*f1D, 1, x, p);
};

double spectra_Diff1D_pp(const double * x, const double* p)
{
  // input parameters
  const auto& type   = p[1]; // function of:   (0) pT,     (1) rapidity,     (3) pT (integrate |y|)
  const auto& varMin = p[5];
  const auto& varMax = p[6];
  
  auto binWidth = (varMax - varMin);
  if (type==2) binWidth *= 2;

  return spectra_Inc1D_pp(x, p) / binWidth;
};

double dist_Inc_pp(TF2& f2D, const UChar_t& id, const std::array<double, 4> range, const double& sqrtS, const double& mass, const double& charge, const bool& useAbsRap)
{
  // check input parameters
  if (range[0]>=range[1]) throw runtime_error("[ERROR] Wrong var1 range values for pp dist inclusive!");
  if (range[2]>=range[3]) throw runtime_error("[ERROR] Wrong var2 range values for pp dist inclusive!");
  
  double res(0);
  if (useAbsRap && range[2]==0) {
    res += evalTF(f2D, {sqrtS, mass, charge, 0}, range[0], range[1], -range[3], range[3], id);
  }
  else if (useAbsRap && range[2]>0) {
    res += evalTF(f2D, {sqrtS, mass, charge, 0}, range[0], range[1], -range[3], -range[2], id);
    res += evalTF(f2D, {sqrtS, mass, charge, 0}, range[0], range[1], range[2], range[3], id);
  }
  else {
    res += evalTF(f2D, {sqrtS, mass, charge, 0}, range[0], range[1], range[2], range[3], id);
  }

  return res;
};

double yield_Inc_pp(const std::array<double, 4> range, const double& sqrtS, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  static TF2* f2D(0);
  if (!f2D) f2D = new TF2("f2DYield", yield_PtRap_pp, range[0], range[1], -range[3], range[3], 4);

  return dist_Inc_pp(*f2D, 0, range, sqrtS, mass, charge, useAbsRap);
};

double yield_Diff_pp(const std::array<double, 4> range, const double& sqrtS, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  auto binWidth = (range[3] - range[2]) * (range[1] - range[0]);
  if (useAbsRap) binWidth *= 2;

  return yield_Inc_pp(range, sqrtS, mass, charge, useAbsRap) / binWidth;
};

double spectra_Inc_pp(const std::array<double, 4> range, const double& sqrtS, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  static TF2* f2D(0);
  if (!f2D) f2D = new TF2("f2DSpectra", spectra_PtRap_pp, range[0], range[1], -range[3], range[3], 4);

  return dist_Inc_pp(*f2D, 1, range, sqrtS, mass, charge, useAbsRap);
};

double spectra_Diff_pp(const std::array<double, 4> range, const double& sqrtS, const double& mass, const int& charge, const bool& useAbsRap=true)
{
  auto binWidth = (range[3] - range[2]) * (range[1] - range[0]);
  if (useAbsRap) binWidth *= 2;

  return spectra_Inc_pp(range, sqrtS, mass, charge, useAbsRap) / binWidth;
};


#endif // #ifndef extrapolationUtils_Yield_pp_h
