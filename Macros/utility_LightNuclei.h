#ifndef utility_LightNuclei_h
#define utility_LightNuclei_h

// ROOT headers
#include "TGraphAsymmErrors.h"
#include "TH3F.h"
#include "TFile.h"
// C++ headers
#include <iostream>
#include <vector>
#include <map>
#include <string>


// Global constants
const std::map<std::string, std::map<std::string, std::string> > LABEL =
  {
   {"",     { {"Deuteron", "d"}, {"Triton", "t"}, {"Helium3", "{}^{3}He"}, {"Helium4", "{}^{4}He"} }},
   {"Anti", { {"Deuteron", "#bar{d}"}, {"Triton", "#bar{t}"}, {"Helium3", "{}^{3}#bar{He}"}, {"Helium4", "{}^{4}#bar{He}"} }},
   {"All",  { {"Deuteron", "d + #bar{d}"}, {"Triton", "t + #bar{t}"}, {"Helium3", "{}^{3}He+{}^{3}#bar{He}"}, {"Helium4", "{}^{4}He+{}^{4}#bar{He}"} }}
  };
const std::map<std::string, int> CHARGE = {{"Pion", 1}, {"Kaon", 1}, {"Proton", 1}, {"Deuteron", 1}, {"Triton", 1}, {"Helium3", 2}, {"Helium4", 2}};
const std::map<std::string, float> MASS = {{"Pion", 0.13957018}, {"Kaon", 0.493677}, {"Proton", 0.9382720813}, {"Deuteron", 1.87561}, {"Triton", 2.80925}, {"Helium3", 2.80923}, {"Helium4", 3.72742}};
const std::map<std::string, std::pair<std::string, std::string> > VAR({{"absRap", {"|y|", ""}}, {"pT", {"p_{T}", "GeV"}}, {"absEta", {"|#eta|", ""}}, {"cent", {"cent.", "%"}}, {"dNdeta", {"dN_{ch}/d#eta", ""}}});
const std::vector<int> COLOR({kRed, kOrange+3, kOrange, kYellow+1, kGreen+7, kGreen+4, kGreen+1, kBlue+7, kBlue+4, kBlue});
const std::vector<int> MARKER({20, 21, 22, 34});


double roundVal(const double& a, const UInt_t& n)
{
  auto f = std::pow(10., n); return std::round(a*f)/f;
};


double ratioUnc(const double& A, const double& B, const double& uncA, const double& uncB)
{
  auto res = std::sqrt((uncA*uncA/A/A) + (uncB*uncB/B/B));
    return std::abs(A/B)*res;
};



double getBeamEnergy(const std::string& col)
{
  if (col.rfind("_")==std::string::npos) return 0;
  auto s = col.substr(col.rfind("_")+1);
  auto c = s.substr(0, s.rfind("TeV"));
  if (c.find("p")!=std::string::npos) c.replace(c.find("p"), 1, ".");
  return std::stod(c);
};


TH3F getHist(const std::string& col, const std::string& par, const std::string& obs, const bool& check=true)
{
  const auto fileName = "../Data/hData_"+col+".root";
  TFile file(fileName.c_str(), "READ");
  if (!file.IsOpen() || file.IsZombie()) throw std::runtime_error("[ERROR] Failed to open file: "+fileName);
  const auto name = "hData_"+par+"_"+obs;
  const auto hist = dynamic_cast<TH3F*>(file.Get(name.c_str()));
  if (check && !hist) throw std::runtime_error("[ERROR] Failed to extract "+name+" from "+fileName);
  if (!hist) return TH3F();
  return *hist;
};


#endif // #ifndef utility_LightNuclei_h
