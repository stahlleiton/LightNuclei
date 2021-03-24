#ifndef test_Extrapolation_C
#define test_Extrapolation_C

// Custom headers
#include "../Fit/utils.h"
#include "../Projection/extrapolationUtils.h"


// Functions
void addExtrapolation(FITFUNC& fitM, const GRAPH& input);


void test_Extrapolation(const std::set<std::string>& pars={})
{
  GRAPH g_ALICE;
  //fillDataGraphs(g_ALICE, type, data_ALICE);
  fillDataGraphs(g_ALICE, {"PbPb_5p02TeV", "PbPb_2p76TeV", "pp_13TeV", "pp_7TeV"}, pars, {"yield", "spectra", "v2"});
	  
  FITFUNC f_ALICE;
  addExtrapolation(f_ALICE, g_ALICE);

  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetFillColor(kWhite);
  
  // Plot the fitted graphs
  plotGraphs(g_ALICE, "Extrapolation", f_ALICE);
};


void addExtrapolation(FITFUNC& fitM, const GRAPH& input)
{
  for (auto& col : input) {
    const auto sqrtSNN = getBeamEnergy(col.first);
    for (auto& chg : col.second) {
      for (auto& par : chg.second) {
	if (par.second.empty()) continue;
	auto& fM = fitM[col.first][chg.first][par.first];
	// Perform extrapolation
	for (auto& obs : par.second) {
	  for (auto& r : obs.second) {
	    auto bin = r.first;
	    for (auto& c : r.second) {
	      bin.add(c.first);
	      auto& fit = fM[obs.first][r.first][c.first];
	      // Define x range
	      const std::string xVar = c.second.GetTitle();
	      const auto xMin = c.second.GetX()[0]-c.second.GetEXlow()[0];
	      const auto xMax = c.second.GetX()[c.second.GetN()-1]+c.second.GetEXhigh()[c.second.GetN()-1];
	      double rapMin(0), rapMax(0);
	      if (bin.has("absRap")) { rapMin = bin.at("absRap").first; rapMax = bin.at("absRap").second; }
	      else if (bin.has("absEta")) { rapMin = bin.at("absEta").first; rapMax = bin.at("absEta").second; } // works since we assume flat eta dependence
	      // Set fit function
	      fit.reset(getProjection1D(col.first, obs.first, {xMin, xMax}, true));
	      const auto charge = (chg.first=="Anti" ? -1 : 1) * CHARGE.at(par.first);
	      if (col.first.rfind("PbPb", 0)==0) {
		if (obs.first=="v2") fit->SetParameters(3, sqrtSNN, MASS.at(par.first), charge, rapMin, rapMax, bin.at("cent").first/100., bin.at("cent").second/100.);
		else fit->SetParameters(1., 3, sqrtSNN, MASS.at(par.first), charge, rapMin, rapMax, bin.at("cent").first/100., bin.at("cent").second/100.);
	      }
	      else if (col.first.rfind("pp", 0)==0) {
		fit->SetParameters(1., 2, sqrtSNN, MASS.at(par.first), charge, rapMin, rapMax);
	      }
	    }
	  }
	}
      }
    }
  }
};


#endif // #ifndef test_Extrapolation_C
