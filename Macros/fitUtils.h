#ifndef fitUtils_h
#define fitUtils_h

// ROOT headers
#include "TF3.h"
#include "TF2.h"
#include "TF1.h"
#include "TMath.h"
// C++ headers
#include <string>


typedef std::map<std::vector<double>, double> CACHE_TYPE;
const auto EPS = 1.E-9;


double computeRapidity(const double& pT, const double& eta, const double& mass)
{
  const auto pTCosH = pT*TMath::CosH(eta), pTSinH = pT*TMath::SinH(eta);
  return TMath::Log((TMath::Sqrt(mass*mass + pTCosH*pTCosH) + pTSinH)/TMath::Sqrt(mass*mass + pT*pT));
};


std::string getStr(const std::string& s, const std::vector<double>& v) {
  auto res = s + " >> ";
  for (const auto& n: v) res += Form(" %g ,", n);
  return res;
};


double evalTF(TF1& f, const std::vector<double>& pars, const double& var)
{
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  fVars.push_back(var);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Eval(var);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};

double evalTF(TF1& f, const std::vector<double>& pars, const double& vMin, const double& vMax, const UChar_t& id=0)
{
  if (vMin>=vMax) throw runtime_error("[ERROR] Wrong var range values for evalTF!");
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  for (auto& v: {vMin, vMax}) fVars.push_back(v);
  if (id!=0) fVars.push_back(id);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Integral(vMin, vMax, EPS);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};

double evalTF(TF2& f, const std::vector<double>& pars, const double& var1, const double& var2)
{
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  for (auto& v: {var1, var2}) fVars.push_back(v);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Eval(var1, var2);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};

double evalTF(TF2& f, const std::vector<double>& pars, const double& v1Min, const double& v1Max, const double& v2Min, const double& v2Max, const UChar_t& id=0)
{
  if (v1Min>=v1Max) throw runtime_error("[ERROR] Wrong var1 range values for evalTF!");
  if (v2Min>=v2Max) throw runtime_error("[ERROR] Wrong var2 range values for evalTF!");
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  for (auto& v: {v1Min, v1Max, v2Min, v2Max}) fVars.push_back(v);
  if (id!=0) fVars.push_back(id);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Integral(v1Min, v1Max, v2Min, v2Max, EPS);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};

double evalTF(TF3& f, const std::vector<double>& pars, const double& var1, const double& var2, const double& var3)
{
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  for (auto& v: {var1, var2, var3}) fVars.push_back(v);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Eval(var1, var2, var3);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};

double evalTF(TF3& f, const std::vector<double>& pars, const double& v1Min, const double& v1Max, const double& v2Min, const double& v2Max, const double& v3Min, const double& v3Max, const UChar_t& id=0)
{
  if (v1Min>=v1Max) throw runtime_error("[ERROR] Wrong var1 range values for evalTF!");
  if (v2Min>=v2Max) throw runtime_error("[ERROR] Wrong var2 range values for evalTF!");
  if (v3Min>=v3Max) throw runtime_error("[ERROR] Wrong var3 range values for evalTF!");
  static std::unique_ptr<CACHE_TYPE> result;
  if (!result) result.reset(new CACHE_TYPE());
  auto fVars = pars;
  for (auto& v: {v1Min, v1Max, v2Min, v2Max, v3Min, v3Max}) fVars.push_back(v);
  if (id!=0) fVars.push_back(id);
  if (result->find(fVars)==result->end()) {
    f.SetParameters(pars.data());
    (*result)[fVars] = f.Integral(v1Min, v1Max, v2Min, v2Max, v3Min, v3Max, EPS);
    if (std::isnan((*result)[fVars])) throw std::runtime_error(getStr("[ERROR] NAN results for", fVars));
  }
  return (*result)[fVars];
};


double blastWaveV2Integrand(const double& pT, const double& phi_s, const double& mass, const double& T_kin, const double& rho_0, const double& rho_a, const double& s_2, const double& n)
{
  // Based on https://arxiv.org/pdf/hep-ph/0101136.pdf
  // check parameters
  if (pT<0) throw runtime_error("[ERROR] Negative pT value for blastWaveV2Integrand!");
  if (mass<=0.05) throw runtime_error("[ERROR] Invalid mass value for blastWaveV2Integrand!");
  if (T_kin<=0) throw runtime_error("[ERROR] Null T_kin value for blastWaveV2Integrand!");

  // transverse mass
  const auto mT = TMath::Sqrt(mass*mass + pT*pT);

  // rho
  const auto rho = rho_0 + rho_a*TMath::Cos(2*phi_s);

  const auto arg1 = std::min(pT * TMath::SinH(rho) / T_kin , 700.); // avoid FPE
  const auto arg2 = mT * TMath::CosH(rho) / T_kin;

  const auto C = (n==1. ? TMath::Cos(2*phi_s) : 1.);
  const auto F = 1. + 2*s_2*TMath::Cos(2*phi_s);
  
  return C * TMath::BesselI(2, arg1) * TMath::BesselK1(arg2) * F;
};

double blastWaveV2_Phi_Integrand(double* x, double* p)
{
  // input variables
  const auto& phi_s = x[0]; // azimuthal angle

  // input parameters
  const auto& pT    = p[0];
  const auto& mass  = p[1];
  const auto& T_kin = p[2];
  const auto& rho_0 = p[3];
  const auto& rho_a = p[4];
  const auto& s_2   = p[5];
  const auto& n     = p[6];

  return blastWaveV2Integrand(pT, phi_s, mass, T_kin, rho_0, rho_a, s_2, n);
};

double blastWaveV2(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& mass  = p[0];
  const auto& T_kin = p[1];
  const auto& rho_0 = p[2];
  const auto& rho_a = p[3];
  const auto& s_2   = p[4];

  static TF1* fIntV2(0);
  if (!fIntV2) fIntV2 = new TF1("fIntBlastWaveV2_Phi_Integrand", blastWaveV2_Phi_Integrand, 0., 2*TMath::Pi(), 7);
  const auto norm = evalTF(*fIntV2, {pT, mass, T_kin, rho_0, rho_a, s_2, 0.}, 0., 2*TMath::Pi());
  const auto res = evalTF(*fIntV2, {pT, mass, T_kin, rho_0, rho_a, s_2, 1.}, 0., 2*TMath::Pi());
  
  return res / norm;
};


double blastWaveIntegrand(const double& r, const double& pT, const double& mass, const double& T_kin, const double& beta, const double& n)
{
  // check parameters
  if (r<0 || r>1) throw runtime_error("[ERROR] Invalid r value for blastWaveIntegrand!");
  if (pT<0) throw runtime_error("[ERROR] Negative pT value for blastWaveIntegrand!");
  if (mass<=0.05) throw runtime_error("[ERROR] Invalid mass value for blastWaveIntegrand!");
  if (T_kin<=0) throw runtime_error("[ERROR] Invalid T_kin value for blastWaveIntegrand!");

  // transverse mass
  const auto mT = TMath::Sqrt(mass*mass + pT*pT);

  // beta and rho (keep within reasonable limits)
  const auto beta_T = std::min(beta * TMath::Power(r, n), 0.99999999999);
  const auto rho = TMath::ATanH(beta_T);

  const auto arg1 = std::min(pT * TMath::SinH(rho) / T_kin , 700.); // avoid FPE
  const auto arg2 = mT * TMath::CosH(rho) / T_kin;

  return r * mT * TMath::BesselI0(arg1) * TMath::BesselK1(arg2);
};

double blastWave_R_Integrand(double* x, double* p)
{
  // input variables
  const auto& r = x[0]; // r/R

  // input parameters
  const auto& pT    = p[0];
  const auto& mass  = p[1];
  const auto& T_kin = p[2];
  const auto& beta  = p[3];
  const auto& n     = p[4];

  return blastWaveIntegrand(r, pT, mass, T_kin, beta, n);
};

double blastWaveYieldFunc(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& mass  = p[0];
  const auto& T_kin = p[1];
  const auto& beta  = p[2];
  const auto& n     = p[3];

  static TF1* fInt(0);
  if (!fInt) fInt = new TF1("fIntBlastWave_R_Integrand", blastWave_R_Integrand, 0., 1., 5);
  const auto res = evalTF(*fInt, {pT, mass, T_kin, beta, n}, 0., 1.);
  
  return (2*TMath::Pi()*pT) * res;
};

double blastWaveYield(const double* x, const double* p)
{
  // implementation of BGBW (dN/dpt)

  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& dNdy  = p[0];
  const auto& mass  = p[1];
  const auto& T_kin = p[2];
  const auto& beta  = p[3];
  const auto& n     = p[4];

  static TF1* fInt(0);
  if (!fInt) fInt = new TF1("fBlastWaveYieldFunc", blastWaveYieldFunc, 0., 100., 4);
  const auto norm = evalTF(*fInt, {mass, T_kin, beta, n}, 0., 100.);
  const auto res = evalTF(*fInt, {mass, T_kin, beta, n}, pT);
  
  return dNdy * res / norm;
};

double blastWaveSpectra(const double* x, const double* p)
{
  return blastWaveYield(x, p) / (2*TMath::Pi()*x[0]);
};


double tsallisBlastWaveIntegrand(const double& r, const double& pT, const double& rap, const double& phi,
				 const double& mass, const double& T_kin, const double& beta_s, const double& q)
{
  // check parameters
  if (r<0 || r>1) throw runtime_error("[ERROR] Invalid r value for tsallisBlastWaveIntegrand!");
  if (pT<0) throw runtime_error("[ERROR] Negative pT value for tsallisBlastWaveIntegrand!");
  if (T_kin==0) throw runtime_error("[ERROR] Null T_kin value for tsallisBlastWaveIntegrand!");
  if (q==1) throw runtime_error("[ERROR] Unit q value for tsallisBlastWaveIntegrand!");

  // transverse mass
  const auto mT = TMath::Sqrt(mass*mass + pT*pT);

  // rho
  const auto rho = TMath::ATanH(beta_s * r);

  const auto d = mT*TMath::CosH(rho)*TMath::CosH(rap) - pT*TMath::SinH(rho)*TMath::Cos(phi);
  const auto arg = 1. + (q-1.) * d / T_kin;
  const auto fct = TMath::Power(arg, 1./(1. - q));
  
  return r * mT * TMath::CosH(rap) * fct;
};

double tsallisBlastWave_RRapPhi_Integrand(double* x, double* p)
{
  // input variables
  const auto& r   = x[0]; // r/R
  const auto& rap = x[1];
  const auto& phi = x[2]; // azimuthal angle
  
  // input parameters
  const auto& mass   = p[0];
  const auto& T_kin  = p[1];
  const auto& beta_s = p[2];
  const auto& q      = p[3];
  const auto& pT     = p[4];

  return tsallisBlastWaveIntegrand(r, pT, rap, phi, mass, T_kin, beta_s, q);
};

double tsallisBlastWaveSpectra_Pt(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& dNdy   = p[0]; // defined within -0.5 < y < 0.5
  const auto& mass   = p[1];
  const auto& T_kin  = p[2];
  const auto& beta_s = p[3];
  const auto& q      = p[4];

  static TF3* fIntTsa = 0;
  if (!fIntTsa) {
    fIntTsa = new TF3("fIntTsallisBlastWave_RRapPhi_Integrand", tsallisBlastWave_RRapPhi_Integrand, 0., 1., -0.5, 0.5, -TMath::Pi(), TMath::Pi(), 5);
    fIntTsa->SetNpx(10000);
    fIntTsa->SetNpy(10000);
    fIntTsa->SetNpz(10000);
  }
    
  const auto res = evalTF(*fIntTsa, {mass, T_kin, beta_s, q, pT}, 0., 1., -0.5, 0.5, -TMath::Pi(), TMath::Pi());
  
  return dNdy * res;
};

double tsallisBlastWaveYield_Pt(const double* x, const double* p)
{
  return (2*TMath::Pi()*x[0]) * tsallisBlastWaveSpectra_Pt(x, p);
};


double tsallisLevy(const double& pT, const double& mass, const double& C, const double& n)
{
  // Based on https://arxiv.org/pdf/2003.03184.pdf
  // check parameters
  if (pT<0) throw runtime_error("[ERROR] Negative pT value for tsallisLevy!");
  if (mass<=0.05) throw runtime_error("[ERROR] Invalid mass value for tsallisLevy!");
  if (C==0) throw runtime_error("[ERROR] Null C value for tsallisLevy!");
  if (n==0) throw runtime_error("[ERROR] Null n value for tsallisLevy!");

  // transverse mass
  const auto mT = TMath::Sqrt(mass*mass + pT*pT);

  const auto K = (n-1)*(n-2)/(n*C*(n*C + mass*(n-2)));
  
  const auto arg = 1. + (mT - mass)/(n*C);
  const auto fct = TMath::Power(arg, -n);

  const auto res = pT * K * fct;
  
  return res;
};

double tsallisLevyYieldFunc(const double* x, const double* p)
{
  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& mass = p[0];
  const auto& C    = p[1];
  const auto& n    = p[2];
  
  return tsallisLevy(pT, mass, C, n);
};

double tsallisLevyYield(const double* x, const double* p)
{
  // implementation of Levy-Tsallis (dN/dpt)

  // input variables
  const auto& pT = x[0];
  
  // input parameters
  const auto& dNdy = p[0];
  const auto& mass = p[1];
  const auto& C    = p[2];
  const auto& n    = p[3];

  static TF1* fInt(0);
  if (!fInt) fInt = new TF1("fIntTsallisLevyYieldFunc", tsallisLevyYieldFunc, 0., 100., 3);
  const auto norm = evalTF(*fInt, {mass, C, n}, 0., 100.);
  const auto res = evalTF(*fInt, {mass, C, n}, pT);
  
  return dNdy * res / norm;
};

double tsallisLevySpectra(const double* x, const double* p)
{
  return tsallisLevyYield(x, p) / (2*TMath::Pi()*x[0]);
};


void setFunction(TF1& f)
{
  const std::string name = f.GetName();
  if (name=="TsallisLevyYield") f.SetFunction(tsallisLevyYield);
  else if (name=="TsallisLevySpectra") f.SetFunction(tsallisLevyYield);
  else if (name=="BlastWaveV2") f.SetFunction(blastWaveV2);
  else if (name=="BlastWaveYield") f.SetFunction(blastWaveYield);
  else if (name=="BlastWaveSpectra") f.SetFunction(blastWaveSpectra);
  else if (name=="TsallisBlastWaveYield") f.SetFunction(tsallisBlastWaveYield_Pt);
  else if (name=="TsallisBlastWaveSpectra") f.SetFunction(tsallisBlastWaveSpectra_Pt);
};


TF1* getFitFunction(const std::string& model, const std::string& obs, const float& xMin, const float& xMax)
{
  TF1* f(0);
  if (model=="TsallisLevy") {
    if (obs=="yield") f = new TF1("TsallisLevyYield", tsallisLevyYield, xMin, xMax, 4);
    else if (obs=="spectra") f = new TF1("TsallisLevySpectra", tsallisLevySpectra, xMin, xMax, 4);
    f->SetParNames("dNdy", "mass", "C", "n");
  }
  else if (model=="BlastWave") {
    if (obs=="yield") f = new TF1("BlastWaveYield", blastWaveYield, xMin, xMax, 5);
    else if (obs=="spectra") f = new TF1("BlastWaveSpectra", blastWaveSpectra, xMin, xMax, 5);
    else if (obs=="v2") f = new TF1("blastWaveV2", blastWaveV2, xMin, xMax, 5);
    if (obs=="v2") f->SetParNames("mass", "T_kin", "rho_0", "rho_a", "s_2");
    else f->SetParNames("dNdy", "mass", "T_kin", "beta_s", "n");
  }
  else if (model=="TsallisBlastWave") {
    if (obs=="yield") f = new TF1("TsallisBlastWaveYield", tsallisBlastWaveYield_Pt, xMin, xMax, 5);
    else if (obs=="spectra") f = new TF1("TsallisBlastWaveSpectra", tsallisBlastWaveSpectra_Pt, xMin, xMax, 5);
    f->SetParNames("dNdy", "mass", "T", "#beta", "q");
  }
  if (!f) throw std::runtime_error("[ERROR] Type of function not found!");
  return f;
};


#endif // #ifndef fitUtils_h
