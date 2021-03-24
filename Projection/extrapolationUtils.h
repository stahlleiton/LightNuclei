#ifndef extrapolationUtils_h
#define extrapolationUtils_h


// Custom headers
#include "extrapolationUtils_Yield_pp.h"
#include "extrapolationUtils_Yield_PbPb.h"
#include "extrapolationUtils_v2.h"


TF1* getProjection1D(const std::string& col, const std::string& obs, const std::array<double, 2>& range, const bool& useDiff)
{
  TF1* f(0);
  if (col.rfind("PbPb", 0)==0) {
    if (obs=="yield") {
      if (useDiff) f = new TF1("Yield_Diff1D_PbPb", yield_Diff1D_PbPb, range[0], range[1], 9);
      else         f = new TF1("Yield_Inc1D_PbPb", yield_Inc1D_PbPb, range[0], range[1], 9);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max", "var2Min", "var2Max");
    }
    else if (obs=="spectra") {
      if (useDiff) f = new TF1("Spectra_Diff1D_PbPb", spectra_Diff1D_PbPb, range[0], range[1], 9);
      else         f = new TF1("Spectra_Inc1D_PbPb", spectra_Inc1D_PbPb, range[0], range[1], 9);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max", "var2Min", "var2Max");
    }
    else if (obs=="v2") {
      f = new TF1("V2_1D_PbPb", v2_1D_PbPb, range[0], range[1], 8);
      f->SetParNames("type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max", "var2Min", "var2Max");
    }
  }
  else if (col.rfind("pp", 0)==0) {
    if (obs=="yield") {
      if (useDiff) f = new TF1("Yield_Diff1D_pp", yield_Diff1D_pp, range[0], range[1], 7);
      else         f = new TF1("Yield_Inc1D_pp", yield_Inc1D_pp, range[0], range[1], 7);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "varMin", "varMax");
    }
    else if (obs=="spectra") {
      if (useDiff) f = new TF1("Spectra_Diff1D_pp", spectra_Diff1D_pp, range[0], range[1], 7);
      else         f = new TF1("Spectra_Inc1D_pp", spectra_Inc1D_pp, range[0], range[1], 7);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "varMin", "varMax");
    }
  }
  if (!f) throw std::runtime_error("[ERROR] Type of function not found for projection 1D!");
  return f;
};


TF2* getProjection2D(const std::string& col, const std::string& obs, const std::array<double, 4>& range, const bool& useDiff)
{
  TF2* f(0);
  if (col.rfind("PbPb", 0)==0) {
    if (obs=="yield") {
      if (useDiff) f = new TF2("Yield_Diff2D_PbPb", yield_Diff2D_PbPb, range[0], range[1], range[2], range[3], 7);
      else         f = new TF2("Yield_Inc2D_PbPb", yield_Inc2D_PbPb, range[0], range[1], range[2], range[3], 7);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max");
    }
    else if (obs=="spectra") {
      if (useDiff) f = new TF2("Spectra_Diff2D_PbPb", spectra_Diff2D_PbPb, range[0], range[1], range[2], range[3], 7);
      else         f = new TF2("Spectra_Inc2D_PbPb", spectra_Inc2D_PbPb, range[0], range[1], range[2], range[3], 7);
      f->SetParNames("scale", "type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max");
    }
    else if (obs=="v2") {
      f = new TF2("V2_2D_PbPb", v2_2D_PbPb, range[0], range[1], range[2], range[3], 6);
      f->SetParNames("type", "sqrt(sNN)", "mass", "charge", "var1Min", "var1Max");
    }
  }
  if (!f) throw std::runtime_error("[ERROR] Type of function not found for projection 2D!");
  return f;
};


TF2* getProjection2D(const std::string& col, const std::string& obs, const double& mass, const int& charge=0)
{
  TF2* f(0);
  if (col.rfind("pp", 0)==0) {
    const auto sqrtS = getBeamEnergy(col);
    if (obs=="yield") {
      f = new TF2("Yield_PtRap_pp", yield_PtRap_pp, 0., 100., -10., 10., 4);
      f->SetParNames("sqrt(s)", "mass", "charge");
      f->SetParameters(sqrtS, mass, charge, 0);
    }
    else if (obs=="spectra") {
      f = new TF2("Spectra_PtRap_pp", spectra_PtRap_pp, 0., 100., -10., 10., 4);
      f->SetParNames("sqrt(s)", "mass", "charge");
      f->SetParameters(sqrtS, mass, charge, 0);
    }
  }
  if (!f) throw std::runtime_error("[ERROR] Type of function not found for projection 2D!");
  return f;
};


TF3* getProjection3D(const std::string& col, const std::string& obs, const double& mass, const int& charge=0)
{
  TF3* f(0);
  if (col.rfind("PbPb", 0)==0) {
    const auto sqrtSNN = getBeamEnergy(col);
    if (obs=="yield") {
      f = new TF3("Yield_PtRapCent_PbPb", yield_PtRapCent_PbPb, 0., 100., -10., 10., 0., 1., 4);
      f->SetParNames("sqrt(sNN)", "mass", "charge");
      f->SetParameters(sqrtSNN, mass, charge, 0);
    }
    else if (obs=="spectra") {
      f = new TF3("Spectra_PtRapCent_PbPb", spectra_PtRapCent_PbPb, 0., 100., -10., 10., 0., 1., 4);
      f->SetParNames("sqrt(sNN)", "mass", "charge");
      f->SetParameters(sqrtSNN, mass, charge, 0);
    }
    else if (obs=="v2") {
      f = new TF3("V2_PtRapCent_PbPb", v2_PtRapCent_PbPb, 0., 100., -10., 10., 0., 1., 1);
      f->SetParNames("mass");
      f->SetParameter(0, mass);
    }
  }
  if (!f) throw std::runtime_error("[ERROR] Type of function not found for projection 3D!");
  return f;
};


double getProjection(const std::string& col, const std::string& obs, const std::array<double, 6>& range, const double& mass, const int& charge, const bool& useDiff, const bool& useAbsRap)
{
  double res(-99.);
  if (col.rfind("PbPb", 0)==0) {
    const auto sqrtSNN = getBeamEnergy(col);
    if (obs=="yield") {
      if (useDiff) res = yield_Diff_PbPb(range, sqrtSNN, mass, charge, useAbsRap);
      else         res = yield_Inc_PbPb(range, sqrtSNN, mass, charge, useAbsRap);
    }
    else if (obs=="spectra") {
      if (useDiff) res = spectra_Diff_PbPb(range, sqrtSNN, mass, charge, useAbsRap);
      else         res = spectra_Inc_PbPb(range, sqrtSNN, mass, charge, useAbsRap);
    }
    else if (obs=="v2") {
      res = v2_Inc_PbPb(range, sqrtSNN, mass, charge, useAbsRap);
    }
  }
  else if (col.rfind("pp", 0)==0) {
    const auto sqrtS = getBeamEnergy(col);
    std::array<double, 4> rng({range[0], range[1], range[2], range[3]});
    if (obs=="yield") {
      if (useDiff) res = yield_Diff_pp(rng, sqrtS, mass, charge, useAbsRap);
      else         res = yield_Inc_pp(rng, sqrtS, mass, charge, useAbsRap);
    }
    else if (obs=="spectra") {
      if (useDiff) res = spectra_Diff_pp(rng, sqrtS, mass, charge, useAbsRap);
      else         res = spectra_Inc_pp(rng, sqrtS, mass, charge, useAbsRap);
    }
  }
  if (res==-99.) throw std::runtime_error("[ERROR] Type of function not found for projection!");
  return res;
};


#endif // #ifndef extrapolationUtils_h
