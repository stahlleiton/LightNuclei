#ifndef fit_LightNuclei_C
#define fit_LightNuclei_C

// ROOT headers
#include "TFitResult.h"
// Custom headers
#include "utils.h"
#include "../Data/data_LightNuclei.h"


// Functions
void fitGraphs(GRAPH& input, FITFUNC& fitM, const std::string& dir);


void fit_LightNuclei(const std::set<std::string>& pars={})
{
  GRAPH g_ALICE;
  fillDataGraphs(g_ALICE, pars, data_ALICE, {"yield", "spectra"});
	  
  FITFUNC f_ALICE;
  fitGraphs(g_ALICE, f_ALICE, "Output");

  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetFillColor(kWhite);
  
  // Plot the fitted graphs
  plotGraphs(g_ALICE, "ALICE", f_ALICE);
};


bool extractFits(FITFUNC1& fM, const GRAPH1& input, const std::string& col, const std::string& chg, const std::string& par, const std::string& obs, const std::string& dir)
{
  const auto fName = dir+"/fitResults_"+col+"_"+chg+par+"_"+obs+".root";
  TFile file(fName.c_str(), "READ"); file.cd();
  if (!file.IsOpen() || file.IsZombie()) return false;
  for (const auto& r : input.at(obs)) {
    for (const auto& c : r.second) {
      const auto name = "fitFunc_"+r.first.str(1)+"_"+c.first.str(1);
      const auto fnc = dynamic_cast<TF1*>(file.Get(name.c_str()));
      if (!fnc) return false;
      fM[obs][r.first][c.first].reset(fnc);
    }
  }
  return true;
};


void storeFits(const FITFUNC1& fM, const std::string& col, const std::string& chg, const std::string& par, const std::string& obs, const std::string& dir)
{
  std::map<std::string, TF1*> fitMap;
  for (auto& r : fM.at(obs)) {
    for (auto& c : r.second) {
      if (!c.second) return;
      if (c.second->GetNumberFitPoints()>c.second->GetNumberFreeParameters() && (c.second->GetNDF()==0 || c.second->GetProb()<0.05)) return;
      const auto name = "fitFunc_"+r.first.str(1)+"_"+c.first.str(1);
      fitMap[name] = c.second.get();
    }
  }
  if (fitMap.empty()) return;
  const auto fName = dir+"/fitResults_"+col+"_"+chg+par+"_"+obs+".root";
  TFile file(fName.c_str(), "RECREATE"); file.cd();
  if (!file.IsOpen() || file.IsZombie()) return;
  for (const auto& f : fitMap) f.second->Write(f.first.c_str());
  file.Write();
  file.Close();
};


void fitGraphs(GRAPH& input, FITFUNC& fitM, const std::string& dir)
{
  gSystem->mkdir(dir.c_str(), kTRUE);
  for (auto& col : input) {
    for (auto& chg : col.second) {
      for (auto& par : chg.second) {
	if (par.second.empty()) continue;
	auto& fM = fitM[col.first][chg.first][par.first];
	for (auto& obs : par.second) {
	  // Extract fits
	  if (extractFits(fM, par.second, col.first, chg.first, par.first, obs.first, dir)) continue;
	  // Perform fits
	  fM.clear();
	  for (auto& r : obs.second) {
	    auto bin = r.first;
	    for (auto& c : r.second) {
	      bin.add(c.first);
	      auto& fit = fM[obs.first][r.first][c.first];
	      // Define fit range
	      const std::string xVar = c.second.GetTitle();
	      const auto xMin = c.second.GetX()[0]-c.second.GetEXlow()[0];
	      const auto xMax = c.second.GetX()[c.second.GetN()-1]+c.second.GetEXhigh()[c.second.GetN()-1];
	      if (xVar!="pT") continue;
	      // Set fit function
	      if (obs.first=="yield" || obs.first=="spectra") {
		if (col.first.rfind("pp_",0)==0) {
		  fit.reset(getFitFunction("TsallisLevy", obs.first, xMin, xMax));
		  fit->SetParameters(2.E-7, MASS.at(par.first), 0.23, 7.);
		  //fit->SetParLimits(0, -1.E4, 1.0); // dNdy
		  fit->FixParameter(1, MASS.at(par.first)); // mass
		  fit->SetParLimits(2, 0.01, 0.99); // C
		  fit->SetParLimits(3, 0.01, 20.); // n
		  /*
		  fit.reset(getFitFunction("TsallisBlastWave", obs.first, xMin, xMax));
		  fit->SetParameters(0.02, 0.8, 0.08, 0.2, 1.14);
		  fit->SetParLimits(1, 0.01, 2.0); // mass
		  fit->SetParLimits(2, 0.0001, 0.99); // T
		  fit->SetParLimits(3, 0.0001, 0.5); // beta
		  fit->SetParLimits(4, 1.001, 1.5); // q
		  */
		}
		else if (col.first.rfind("PbPb_",0)==0) {
		  fit.reset(getFitFunction("BlastWave", obs.first, xMin, xMax));
		  if (par.first=="Deuteron") {
		    fit->SetParameters(0.08, 0, 0.090, 0.90, 0.74);
		    fit->SetParLimits(0, 0.001, 0.150); // norm
		    fit->SetParLimits(2, 0.050, 0.220); // T_kin
		    fit->SetParLimits(3, 0.55, 0.95); // beta_s
		    fit->SetParLimits(4, 0.40, 3.3); // n
		    /*
		      fit->SetParameters(1.E1, 0, 0.3, 0.6, 0.01);
		      fit->SetParLimits(0, 0.0, 1.E5); // norm
		      fit->SetParLimits(2, 0.01, 0.99); // T_kin
		      fit->SetParLimits(3, 0.0, 0.99); // beta_s
		      fit->SetParLimits(4, 0.0, 3.0); // n
		    */
		  }
		  else if (par.first=="Helium3") {
		    fit->SetParameters(0.0002, 0, 0.075, 0.88, 0.66);
		    fit->SetParLimits(0, 0.00001, 0.00150); // norm
		    fit->SetParLimits(2, 0.050, 0.400); // T_kin
		    fit->SetParLimits(3, 0.55, 0.95); // beta_s
		    fit->SetParLimits(4, 0.40, 3.3); // n
		  }
		  else if (par.first=="Helium4") {
		    fit->SetParameters(0.0002, 0, 0.075, 0.88, 0.66);
		    fit->SetParLimits(0, 0.0000, 0.00150); // norm
		    fit->SetParLimits(2, 0.050, 0.400); // T_kin
		    fit->SetParLimits(3, 0.55, 0.95); // beta_s
		    fit->SetParLimits(4, 0.40, 3.3); // n
		  }
		  else if (par.first=="Triton") {
		    fit->SetParameters(0.0001, 0, 0.120, 0.75, 0.5);
		    fit->SetParLimits(0, 0.000001, 0.00050); // norm
		    fit->SetParLimits(2, 0.005, 0.350); // T_kin
		    fit->SetParLimits(3, 0.50, 0.94); // beta_s
		    fit->SetParLimits(4, 0.0, 3.0); // n
		    fit->FixParameter(3, 0.75);
		  }
		  fit->FixParameter(1, MASS.at(par.first));
		}
	      }
	      else if (obs.first=="v2") {
		if (col.first.rfind("PbPb_",0)==0) {
		  fit.reset(getFitFunction("BlastWave", obs.first, xMin, xMax));
		  if (par.first=="Deuteron") {
		    fit->SetParameters(0, 0.100, 0.81, 0.0310, 0.0321);
		    fit->SetParLimits(1, 0.090, 0.140); // T_kin
		    fit->SetParLimits(2, 0.50, 1.00); // rho_0
		    fit->SetParLimits(3, 0.005, 0.040); // rho_a
		    fit->SetParLimits(4, 0.000, 0.200); // s_2
		    //fit->FixParameter(4, 0.062); //s_2
		    //fit->FixParameter(1, 0.1); // T_kin
		    //fit->FixParameter(2, 0.81);
		  }
		  else if (par.first=="Helium3") {
		    fit->SetParameters(0, 0.106, 0.88, 0.64, 0.0406);
		    fit->SetParLimits(1, 0.005, 1.600); // T_kin
		    fit->SetParLimits(3, 0.005, 1.500); // rho_a
		    fit->SetParLimits(4, 0.000, 0.200); // s_2
		    //fit->FixParameter(4, 0.0); //s_2
		  }
		  fit->FixParameter(0, MASS.at(par.first));
		}
	      }
	      if (!fit) continue;
	      // Perform fits
	      const auto pars = fit->GetParameters();
	      for (size_t i=0; i<5; i++) {
		std::cout << col.first << " , " << chg.first+par.first << " , " << obs.first << " , " << r.first.str() << " , " << c.first.str() << " , " << xMin << " , " << xMax << std::endl;
		const auto& fitResult = c.second.Fit(fit.get(), "RN0S", "", xMin, xMax);
		std::cout << col.first << " , " << chg.first+par.first << " , " << obs.first << " , " << r.first.str() << " , " << c.first.str() << std::endl;
		const auto& edm = fitResult->Edm();
		const auto& chi2 = fitResult->Chi2();
		if (edm>0 && edm<1.E-2 && chi2<10.) break;
		if (edm>0 && edm<1. && chi2<1.5) break;
		fit->SetChisquare(100.);
		fit->SetParameters(pars);
	      }
	      std::cout << "FIT RESULT: " << fit->GetChisquare() << std::endl;
	      //if (fit->GetChisquare()>50) throw std::runtime_error("STOP");
	    }
	  }
	  // Store fits
	  storeFits(fM, col.first, chg.first, par.first, obs.first, dir);
	}
      }
    }
  }
};


#endif // #ifndef fit_LightNuclei_C
