#ifndef extrapolationUtils_v2_h
#define extrapolationUtils_v2_h


// Custom headers
#include "../Macros/utility_LightNuclei.h"
#include "extrapolationUtils_Yield_PbPb.h"


double v2_PtEtaCent_PbPb(const double& pT, const double& eta, const double& cnt, const double& mass)
{
  // check input parameters
  if (cnt<0 || pT<0) throw runtime_error("[ERROR] Negative centrality or pT value for v2!");
  if (cnt>1) throw runtime_error("[ERROR] Invalid centrality value for v2!");
  if (mass<=1.5) throw runtime_error("[ERROR] Invalid mass value for v2!");
  
  // source: Deuteron v2 in PbPb 5.02 TeV, Phys. Lett. B 805 (2020) 135414
  static TH3F* hV2(0);
  if (!hV2) { hV2 = new TH3F(); *hV2 = getHist("PbPb_5p02TeV", "AllDeuteron", "v2"); }
  
  // Use deuteron v2 and scale as v2_x(pT) = v2_d(F * pT) / F where F = A_d / A_x;
  double A(2);
  if (mass>2.5) A = 3; // Helium3 and Triton
  if (mass>3.0) A = 4; // Helium4
  
  const auto bin = hV2->FindBin((2/A)*pT, 100*cnt, std::abs(eta));
  const auto res = (A/2)*hV2->GetBinContent(bin);

  return res;
};

double v2_PtRapCent_PbPb(const double& pT, const double& rap, const double& cnt, const double& mass)
{
  const auto eta = TMath::ASinH(TMath::SinH(rap) * TMath::Sqrt(1. + mass*mass/pT/pT));
  return v2_PtEtaCent_PbPb(pT, eta, cnt, mass);
};

double v2_PtRapCent_PbPb(const double* x, const double* p)
{
  // input variables
  const auto& pT  = x[0];
  const auto& rap = x[1];
  const auto& cnt = x[2];
  
  // input parameters
  const auto& mass = p[0];

  return v2_PtRapCent_PbPb(pT, rap, cnt, mass);
};

double v2Yield_PtRapCent_PbPb(const double* x, const double* p)
{
  // input variables
  const auto& t = p[3];
  const auto& pT  = t==0 ? x[0] : (t==1 ? p[4] : (t==2 ? x[0] : (t==3 ? x[0] : (t==4 ? x[0] : (t==5 ? p[4] : p[4])))));
  const auto& rap = t==0 ? x[1] : (t==1 ? x[0] : (t==2 ? p[4] : (t==3 ? x[1] : (t==4 ? p[4] : (t==5 ? x[0] : p[5])))));
  const auto& cnt = t==0 ? x[2] : (t==1 ? x[1] : (t==2 ? x[1] : (t==3 ? p[4] : (t==4 ? p[5] : (t==5 ? p[5] : x[0])))));
  if (t<0 || t>6) throw runtime_error("[ERROR] Wrong type value for v2Yield!");
  
  // input parameters
  const auto& mass    = p[0];
  const auto& charge  = p[1];
  const auto& sqrtSNN = p[2];

  const auto v2 = v2_PtRapCent_PbPb(pT, rap, cnt, mass);
  const auto yield = yield_PtRapCent_PbPb(pT, rap, cnt, sqrtSNN, mass, charge);

  return v2 * yield;
};

double v2_2D_PbPb(const double * x, const double* p)
{
  // input parameters
  const auto& type    = p[0]; // average over: (0) pT,    (1) rapidity,    (2) centrality,    (3) |rapidity|
  const auto& sqrtSNN = p[1];
  const auto& mass    = p[2];
  const auto& charge  = p[3];
  const auto& var1Min = p[4];
  const auto& var1Max = p[5];
  
  // check input parameters
  if (type<0 || type>3) throw runtime_error("[ERROR] Wrong type value for v2 2D!");
  if (var1Min>=var1Max) throw runtime_error("[ERROR] Wrong range values for v2 2D!");

  static TF1 *fV2Yield(0), *fYield(0);
  if (!fV2Yield) fV2Yield = new TF1("fV2Yield", v2Yield_PtRapCent_PbPb, -var1Max, var1Max, 6);
  if (!fYield) fYield = new TF1("fYield", yield_PtRapCent_PbPb, -var1Max, var1Max, 6);

  double res(0), norm(0);
  if (type<3) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, type+4, x[0], x[1]}, var1Min, var1Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, type+4, x[0], x[1]}, var1Min, var1Max);
  }
  else if (type==3 && var1Min==0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 5, x[0], x[1]}, -var1Max, var1Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 5, x[0], x[1]}, -var1Max, var1Max);
  }
  else if (type==3 && var1Min>0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 5, x[0], x[1]}, -var1Max, -var1Min);
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 5, x[0], x[1]}, var1Min, var1Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 5, x[0], x[1]}, -var1Max, -var1Min);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 5, x[0], x[1]}, var1Min, var1Max);
  }

  return res / norm;
};

double v2_1D_PbPb(const double* x, const double* p)
{
  // input parameters
  const auto& type    = p[0]; // function of:   (0) pT,     (1) rapidity,     (2) centrality,     (3) pT (integrate |y|),     (4) centrality (integrate |y|)
  const auto& sqrtSNN = p[1];
  const auto& mass    = p[2];
  const auto& charge  = p[3];
  const auto& var1Min = p[4];
  const auto& var1Max = p[5];
  const auto& var2Min = p[6];
  const auto& var2Max = p[7];
  
  // check input parameters
  if (type<0 || type>4) throw runtime_error("[ERROR] Wrong type value for v2 1D!");
  if (var1Min>=var1Max) throw runtime_error("[ERROR] Wrong var1 range values for v2 1D!");
  if (var2Min>=var2Max) throw runtime_error("[ERROR] Wrong var2 range values for v2 1D!");
  
  static TF2 *fV2Yield(0), *fYield(0);
  if (!fV2Yield) fV2Yield = new TF2("fV2Yield", v2Yield_PtRapCent_PbPb, -var1Max, var1Max, -var2Max, var2Max, 5);
  if (!fYield) fYield = new TF2("fYield", yield_PtRapCent_PbPb, -var1Max, var1Max, -var2Max, var2Max, 5);
  
  double res(0), norm(0);
  if (type<3) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, type+1, x[0]}, var1Min, var1Max, var2Min, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, type+1, x[0]}, var1Min, var1Max, var2Min, var2Max);
  }
  else if (type==3 && var1Min==0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 1, x[0]}, -var1Max, var1Max, var2Min, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 1, x[0]}, -var1Max, var1Max, var2Min, var2Max);
  }
  else if (type==3 && var1Min>0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 1, x[0]}, -var1Max, -var1Min, var2Min, var2Max);
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 1, x[0]}, var1Min, var1Max, var2Min, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 1, x[0]}, -var1Max, -var1Min, var2Min, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 1, x[0]}, var1Min, var1Max, var2Min, var2Max);
  }
  else if (type==4 && var2Min==0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 3, x[0]}, var1Min, var1Max, -var2Max, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, -var2Max, var2Max);
  }
  else if (type==4 && var2Min>0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 3, x[0]}, var1Min, var1Max, -var2Max, -var2Min);
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 3, x[0]}, var1Min, var1Max, var2Min, var2Max);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, -var2Max, -var2Min);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 3, x[0]}, var1Min, var1Max, var2Min, var2Max);
  }
  
  return res / norm;
};

double v2_Inc_PbPb(const std::array<double, 6> range, const double& sqrtSNN, const double& mass, const double& charge, const bool& useAbsRap=true)
{
  // check input parameters
  if (range[0]>=range[1]) throw runtime_error("[ERROR] Wrong var1 range values for v2 inclusive!");
  if (range[2]>=range[3]) throw runtime_error("[ERROR] Wrong var2 range values for v2 inclusive!");
  if (range[4]>=range[5]) throw runtime_error("[ERROR] Wrong var3 range values for v2 inclusive!");
  
  static TF3 *fV2Yield(0), *fYield(0);
  if (!fV2Yield) fV2Yield = new TF3("fV2Yield", v2Yield_PtRapCent_PbPb, range[0], range[1], -range[3], range[3], range[4], range[5], 4);
  if (!fYield) fYield = new TF3("fYield", yield_PtRapCent_PbPb, range[0], range[1], -range[3], range[3], range[4], range[5], 4);
  
  double res(0), norm(0);
  if (useAbsRap && range[2]==0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 0}, range[0], range[1], -range[3], range[3], range[4], range[5]);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 0}, range[0], range[1], -range[3], range[3], range[4], range[5]);
  }
  else if (useAbsRap && range[2]>0) {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 0}, range[0], range[1], -range[3], -range[2], range[4], range[5]);
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 0}, range[0], range[1], range[2], range[3], range[4], range[5]);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 0}, range[0], range[1], -range[3], -range[2], range[4], range[5]);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 0}, range[0], range[1], range[2], range[3], range[4], range[5]);
  }
  else {
    res += evalTF(*fV2Yield, {mass, charge, sqrtSNN, 0}, range[0], range[1], range[2], range[3], range[4], range[5]);
    norm += evalTF(*fYield, {sqrtSNN, mass, charge, 0}, range[0], range[1], range[2], range[3], range[4], range[5]);
  }

  return res / norm;
};


#endif // #ifndef extrapolationUtils_v2_h
