#ifndef projectedData_LightNuclei_h
#define projectedData_LightNuclei_h

#include "../Macros/bin.h"


// Bins for total yield vs |y| in PbPb
const CONTBIN totalYield_vs_AbsRap_PbPb_set1 =
  {
   {
    {{"pT", {0.0, 100.}}},
    {
     {
      {{"cent", {0., 20.}}},
      {
       { {"absRap", {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0}} }
      }
     }
    }
   }
  };

const CONTBIN totalYield_vs_AbsRap_PbPb_set2 =
  {
   {
    {{"pT", {0.0, 100.}}},
    {
     {
      {{"cent", {0., 20.}}},
      {
       { {"absRap", {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0}} }
      }
     }
    }
   }
  };


// Bins for v2 vs pT
const CONTBIN v2_vs_pT_A2 =
  {
   {
    {{"cent", {20., 50.}}},
    {
     {
      {{"absRap", {0., 1.}}},
      {
       { {"pT", {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0}} }
      }
     },
     {
      {{"absRap", {1., 2.}}},
      {
       { {"pT", {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0}} }
      }
     },
     {
      {{"absRap", {2., 3.}}},
      {
       { {"pT", {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0}} }
      }
     }
    }
   }
  };

const CONTBIN v2_vs_pT_A3 =
  {
   {
    {{"cent", {10., 50.}}},
    {
     {
      {{"absRap", {0., 1.}}},
      {
       { {"pT", {1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.5, 5.3, 6.0, 7.5}} }
      }
     },
     {
      {{"absRap", {1., 2.}}},
      {
       { {"pT", {1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.5, 5.3, 6.0, 7.5}} }
      }
     },
     {
      {{"absRap", {2., 3.}}},
      {
       { {"pT", {1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.5, 5.3, 7.5}} }
      }
     },
     {
      {{"absRap", {0., 2.}}},
      {
       { {"pT", {1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.5, 5.3, 6.0, 7.5}} }
      }
     }
    }
   }
  };

const CONTBIN v2_vs_pT_A4 =
  {
   {
    {{"cent", {10., 50.}}},
    {
     {
      {{"absRap", {0., 1.}}},
      {
       { {"pT", {1.6, 3.2, 4.0, 4.8, 5.6, 7.0, 10.0}} }
      }
     },
     {
      {{"absRap", {1., 2.}}},
      {
       { {"pT", {1.6, 3.2, 4.0, 4.8, 5.6, 7.0, 10.0}} }
      }
     },
     {
      {{"absRap", {2., 3.}}},
      {
       { {"pT", {1.6, 3.2, 5.6, 7.0, 10.}} }
      }
     },
     {
      {{"absRap", {0., 2.}}},
      {
       { {"pT", {1.6, 2.4, 3.2, 4.0, 4.8, 5.6, 7.0, 10.0}} }
      }
     }
    }
   }
  };


// Bins for total yield in pp
const CONTBIN totalYield_pp =
  {
   {
    {{"cent", {0.0, 100.}}},
    {
     {
      {{"absRap", {0., 0.5}}},
      {
       { {"pT", {0., 100.}} }
      }
     }
    }
   }
  };


// Set binning for projected data
const DATA_CONTBIN binning_ProjectedData_V4 =
  {
   {"pp_14TeV",
    {
     {"", //Particle
      {
       {"Triton",
	{
	 {"yield", totalYield_pp }
	}
       }
      }
     },
     {"Anti", //Antiparticle
      {
       {"Triton",
	{
	 {"yield", totalYield_pp }
	}
       }
      }
     }
    }
   }
  };

const DATA_CONTBIN binning_ProjectedData =
  {
   {"PbPb_5p5TeV",
    {
     {"", //Particle
      {
       {"Deuteron",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       },
       {"Helium3",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       },
       {"Helium4",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set2 }
	}
       },
       {"Triton",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       }
      }
     },
     {"Anti", //Antiparticle
      {
       {"Deuteron",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       },
       {"Helium3",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       },
       {"Helium4",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set2 }
	}
       },
       {"Triton",
	{
	 {"yield", totalYield_vs_AbsRap_PbPb_set1 }
	}
       }
      }
     },
     {"All", //Particle + antiparticle
      {
       {"Helium3",
	{
	 {"v2", v2_vs_pT_A3},
	}
       },
       {"Helium4",
	{
	 {"v2", v2_vs_pT_A4},
	}
       }
      }
     }
    }
   },
   {"pp_14TeV",
    {
     {"", //Particle
      {
       {"Deuteron",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Helium3",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Helium4",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Triton",
	{
	 {"yield", totalYield_pp }
	}
       }
      }
     },
     {"Anti", //Antiparticle
      {
       {"Deuteron",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Helium3",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Helium4",
	{
	 {"yield", totalYield_pp }
	}
       },
       {"Triton",
	{
	 {"yield", totalYield_pp }
	}
       }
      }
     }
    }
   }
  };


#endif // #ifndef projectedData_LightNuclei_h
