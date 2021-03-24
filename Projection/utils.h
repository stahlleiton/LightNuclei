#ifndef LightNuclei_Projection_utils_h
#define LightNuclei_Projection_utils_h

// Custom headers
#include "../Macros/bin.h"
#include "../Macros/utility_LightNuclei.h"
#include "extrapolationUtils.h"
#include "projectedData_Binning.h"
// ROOT headers
#include "TFile.h"
// C++ headers
#include <iostream>
#include <vector>
#include <map>
#include <string>


// Type definition
typedef std::map< BINF , std::map< std::string , std::map< std::string , TGraphAsymmErrors > > > GRAPH1;
typedef std::map< std::string , std::map< BINF , GRAPH1 > > GRAPH2;
typedef std::map< std::string , GRAPH2 > GRAPH;
typedef std::map<std::string, GRAPH> GRAPHM;
typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, bool> > > > BOOLCONT;


// Functions


BOOLCONT extractGraph(GRAPH& gM, const std::string& type, const std::string& sel, const DATA_CONTBIN& data)
{
  BOOLCONT missing;
  bool foundAll(true);
  for (const auto& col : data) {
    for (const auto& chg : col.second) {
      for (const auto& par : chg.second) {
	for (const auto& obs : par.second) {
	  // Extract graphs in file
	  const auto fileName = "Output/GRAPH_"+type+"/GRAPH_"+type+"_"+chg.first+par.first+"_"+col.first+"_"+obs.first+"_"+sel.c_str()+".root";
	  TFile file(fileName.c_str(), "READ");
	  bool found(true);
	  if (file.IsOpen() && !file.IsZombie()) {
	    for (const auto& r : obs.second) {
	      for (const auto& c : r.second) {
		const auto name = "graph_"+col.first+"_"+chg.first+par.first+"_"+obs.first+"_"+r.first.str(1)+"_"+c.first.str(1);
		const auto& graph = dynamic_cast<TGraphAsymmErrors*>(file.Get(name.c_str()));
		if (graph) gM[col.first][obs.first][r.first][c.first][chg.first][par.first] = *graph;
		else { found = false; break; }
	      }
	      if (!found) break;
	    }
	    file.Close();
	  }
	  else found = false;
	  if (found) std::cout << "[INFO] Found all graphs in " << fileName << std::endl;
	  missing[col.first][chg.first][par.first][obs.first] = !found;
	  if (!found) foundAll = false;
	}
      }
    }
  }
  if (foundAll) missing.clear();
  return missing;
};


#endif // #ifndef LightNuclei_Projection_utils_h
