#include <map>
#include <string>

std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > > FUNCMAP_ =
  {
   // BTL region
   { "BTL",
     {
      { "Deuteron",
        {
         { "dEdx",
	   {
	    { "Mean", "1.96031*(1+3.02004/x/x+0.645994*log(x)/x/x+0.0806513*log(x))" },
	    { "Error", "min(0.28369*(1.+-0.340068/x+2.63488/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.17304*(1+2.67008/x/x+0.895547*log(x)/x/x+0.0802195*log(x))" },
	    { "Error", "min(0.344465*(1.+-0.166354/x+1.69837/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.34749*(1+2.43048/x/x+0.67152*log(x)/x/x+0.0655253*log(x))" },
	    { "Error", "min(0.390026*(1.+0.130539/x+0.797627/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.92314*1.92314+x*x)+-0.00829856*log(x))/x" },
	    { "Error", "min(0.536286*(1+log(x)*(-0.678549+-0.00621141*x*x)+0.034797*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Helium3",
        {
         { "dEdx",
	   {
	    { "Mean", "5.79232*(1+0.260312/x/x+0.161423*log(x)/x/x+0.00999868*log(x))" },
	    { "Error", "min(0.829191*(1.+-0.268972/x+1.0859/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "6.39809*(1+0.27072/x/x+0.106738*log(x)/x/x+0.00592375*log(x))" },
	    { "Error", "min(0.830905*(1.+-0.0952092/x+0.820337/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "6.60005*(1+0.311164/x/x+-0.0286522*log(x)/x/x+0.00138079*log(x))" },
	    { "Error", "min(0.851314*(1.+0.00120642/x+0.822314/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.4789*1.4789+x*x)+-0.00959373*log(x))/x" },
	    { "Error", "min(0.315024*(1+log(x)*(-0.64854+-0.00485013*x*x)+0.0323571*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Helium4",
        {
         { "dEdx",
	   {
	    { "Mean", "5.53002*(1+0.356349/x/x+0.548526*log(x)/x/x+0.026021*log(x))" },
	    { "Error", "min(0.818758*(1.+0.00853419/x+0.914788/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "6.10562*(1+0.379018/x/x+0.474894*log(x)/x/x+0.0218456*log(x))" },
	    { "Error", "min(0.822955*(1.+0.108531/x+0.772818/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "6.25071*(1+0.389879/x/x+0.438522*log(x)/x/x+0.0192591*log(x))" },
	    { "Error", "min(0.854559*(1.+0.0421126/x+1.06265/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.92885*1.92885+x*x)+-0.00979419*log(x))/x" },
	    { "Error", "min(0.535387*(1+log(x)*(-0.618446+-0.0039706*x*x)+0.0232663*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Kaon",
        {
         { "dEdx",
	   {
	    { "Mean", "2.36912*(1+0.141667/x/x+-0.136268*log(x)/x/x+0.0255214*log(x))" },
	    { "Error", "min(0.250634*(1.+0.021398/x+0.289117/x/x),0.332)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.66357*(1+0.15377/x/x+-0.207741*log(x)/x/x+0.0162399*log(x))" },
	    { "Error", "min(0.304427*(1.+0.173729/x+0.157061/x/x),0.414)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.74322*(1+0.180252/x/x+-0.204191*log(x)/x/x+0.014962*log(x))" },
	    { "Error", "min(0.332246*(1.+0.259724/x+0.136391/x/x),0.685)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(0.585039*0.585039+x*x)+-0.00761037*log(x))/x" },
	    { "Error", "min(0.0609854*(1+log(x)*(-0.669705+0.000157112*x*x)+0.0931937*x*x)/x/x,0.050)" },
	   }
	 },
        }
      },
      { "Pion",
        {
         { "dEdx",
	   {
	    { "Mean", "2.50849*(1+0.00189057*x)" },
	    { "Error", "min(0.251715*(1.+0.145213/(x+2.82047e-11)),0.300)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.77042*(1+0.00149014*x)" },
	    { "Error", "min(0.308329*(1.+0.161878/(x+7.85383e-11)),0.369)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.85989*(1+0.00110288*x)" },
	    { "Error", "min(0.337683*(1.+0.219922/(x+1.97231e-11)),0.476)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "sqrt(0.997774+0.02/x/x)" },
	    { "Error", "min(0.00563828*(1+-0.00491058*x),0.008)" },
	   }
	 },
        }
      },
      { "Proton",
        {
         { "dEdx",
	   {
	    { "Mean", "2.30327*(1+0.738838/x/x+-0.293608*log(x)/x/x+0.0277897*log(x))" },
	    { "Error", "min(0.230421*(1.+0.331153/x+0.680514/x/x),0.781)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.53927*(1+0.733537/x/x+-0.250319*log(x)/x/x+0.0256227*log(x))" },
	    { "Error", "min(0.288783*(1.+0.412346/x+0.277028/x/x),0.866)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.62202*(1+0.716126/x/x+-0.190874*log(x)/x/x+0.0233202*log(x))" },
	    { "Error", "min(0.316375*(1.+0.53748/x+0.156494/x/x),0.658)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.00148*1.00148+x*x)+-0.00693493*log(x))/x" },
	    { "Error", "min(0.1*(1+log(x)*(-0.302522+0.00616925*x*x)+0.034885*x*x)/x/x,0.129)" },
	   }
	 },
        }
      },
      { "Triton",
        {
         { "dEdx",
	   {
	    { "Mean", "1.62015*(1+4.66715/x/x+4.92688*log(x)/x/x+0.15069*log(x))" },
	    { "Error", "min(0.286484*(1.+-0.710248/x+5.73213/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "1.83732*(1+4.05723/x/x+4.79025*log(x)/x/x+0.140063*log(x))" },
	    { "Error", "min(0.344342*(1.+-0.329482/x+3.72654/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.02159*(1+3.79802/x/x+3.93521*log(x)/x/x+0.116199*log(x))" },
	    { "Error", "min(0.399554*(1.+-0.354045/x+2.97666/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(2.89309*2.89309+x*x)+-0.0125391*log(x))/x" },
	    { "Error", "min(1.27228*(1+log(x)*(-0.672347+-0.00553714*x*x)+0.0249916*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
     }
   },
   // ETL region
   { "ETL",
     {
      { "Deuteron",
        {
         { "dEdx",
	   {
	    { "Mean", "2.03098*(1+0.444951/x/x+0.993581*log(x)/x/x+0.0413861*log(x))" },
	    { "Error", "min(0.127074*(1.+-1.23243/x+5.20104/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.21518*(1+1.31936/x/x+0.443112*log(x)/x/x+0.0395841*log(x))" },
	    { "Error", "min(0.189601*(1.+-2.20128/x+11.4857/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.2858*(1+1.30624/x/x+0.627367*log(x)/x/x+0.0443814*log(x))" },
	    { "Error", "min(0.251782*(1.+-1.09503/x+5.45087/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.92224*1.92224+x*x)+-0.00580526*log(x))/x" },
	    { "Error", "min(0.651749*(1+log(x)*(-0.638896+-0.00434826*x*x)+0.0205403*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Helium3",
        {
         { "dEdx",
	   {
	    { "Mean", "2.9599*(1+0.168964/x/x+0.0401449*log(x)/x/x+0.000548417*log(x))" },
	    { "Error", "min(0.333049*(1.+-0.183676/x+1.75928/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "4.10713*(1+0.265627/x/x+-0.0419816*log(x)/x/x+-0.00668154*log(x))" },
	    { "Error", "min(0.838455*(1.+0.0846674/x+0.886645/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "4.54891*(1+0.269442/x/x+-0.127015*log(x)/x/x+-0.00850254*log(x))" },
	    { "Error", "min(0.703937*(1.+-0.145216/x+1.79489/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.43408*1.43408+x*x)+-0.003689*log(x))/x" },
	    { "Error", "min(0.38562*(1+log(x)*(-0.721372+-0.00704187*x*x)+0.0323578*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Helium4",
        {
         { "dEdx",
	   {
	    { "Mean", "2.95919*(1+0.246532/x/x+0.0742251*log(x)/x/x+0.000397814*log(x))" },
	    { "Error", "min(0.337129*(1.+-0.186709/x+2.23903/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "3.85918*(1+0.263763/x/x+0.549037*log(x)/x/x+0.0142678*log(x))" },
	    { "Error", "min(0.841157*(1.+-0.196699/x+1.99589/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "4.42189*(1+0.317186/x/x+0.229182*log(x)/x/x+0.000734983*log(x))" },
	    { "Error", "min(0.683652*(1.+0.170088/x+1.47462/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(1.91678*1.91678+x*x)+-0.00629957*log(x))/x" },
	    { "Error", "min(0.640444*(1+log(x)*(-0.64165+-0.00460962*x*x)+0.0212883*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Kaon",
        {
         { "dEdx",
	   {
	    { "Mean", "2.29641*(1+0.0740719/x/x+-0.0943651*log(x)/x/x+0.00210887*log(x))" },
	    { "Error", "min(0.0871965*(1.+1.74758/x+-0.7/x/x),0.204)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.52472*(1+0.0816454/x/x+-0.111076*log(x)/x/x+-0.00666966*log(x))" },
	    { "Error", "min(0.140753*(1.+1.15822/x+-0.217901/x/x),0.334)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.62626*(1+0.080205/x/x+-0.121248*log(x)/x/x+-0.00237001*log(x))" },
	    { "Error", "min(0.221507*(1.+0.354084/x+0.0016693/x/x),0.336)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(0.472675*0.472675+x*x)+-0.000504868*log(x))/x" },
	    { "Error", "min(0.057614*(1+log(x)*(-0.85823+-0.0100394*x*x)+0.0950164*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
      { "Pion",
        {
         { "dEdx",
	   {
	    { "Mean", "2.32815*(1+-0.000111344*x)" },
	    { "Error", "min(0.0860073*(1.+1.76937/(x+1.30943)),0.167)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.53153*(1+-0.000913883*x)" },
	    { "Error", "min(0.143697*(1.+1.74953/(x+2.54064)),0.233)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.64549*(1+-0.000418609*x)" },
	    { "Error", "min(0.0591433*(1.+315.713/(x+95.6565)),0.288)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "sqrt(0.999768+0.02/x/x)" },
	    { "Error", "min(0.0035*(1+0.0014066*x),0.200)" },
	   }
	 },
        }
      },
      { "Proton",
        {
         { "dEdx",
	   {
	    { "Mean", "2.2472*(1+0.353374/x/x+-0.133121*log(x)/x/x+0.00728233*log(x))" },
	    { "Error", "min(0.0878658*(1.+2.12522/x+-0.320175/x/x),0.321)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "2.50741*(1+0.448156/x/x+-0.295403*log(x)/x/x+-0.00605313*log(x))" },
	    { "Error", "min(0.136546*(1.+1.38709/x+0.783348/x/x),0.532)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.61316*(1+0.47098/x/x+-0.332115*log(x)/x/x+-0.00309796*log(x))" },
	    { "Error", "min(0.21894*(1.+0.484796/x+0.45723/x/x),0.507)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(0.939109*0.939109+x*x)+-0.000787422*log(x))/x" },
	    { "Error", "min(0.1*(1+log(x)*(-0.639383+-0.0052219*x*x)+0.0540903*x*x)/x/x,0.094)" },
	   }
	 },
        }
      },
      { "Triton",
        {
         { "dEdx",
	   {
	    { "Mean", "1.91501*(1+-0.518267/x/x+3.12393*log(x)/x/x+0.0569626*log(x))" },
	    { "Error", "min(0.131486*(1.+-2.22695/x+12.0178/x/x),2.000)" },
	   }
	 },
         { "dEdx1",
	   {
	    { "Mean", "1.9883*(1+0.976406/x/x+3.10859*log(x)/x/x+0.0734723*log(x))" },
	    { "Error", "min(0.188039*(1.+-3.06989/x+24.4769/x/x),2.000)" },
	   }
	 },
         { "dEdx2",
	   {
	    { "Mean", "2.08826*(1+1.2746/x/x+3.00164*log(x)/x/x+0.0714078*log(x))" },
	    { "Error", "min(0.233949*(1.+-0.505466/x+9.52472/x/x),2.000)" },
	   }
	 },
         { "invBeta",
	   {
	    { "Mean", "(sqrt(2.91716*2.91716+x*x)+-0.0106045*log(x))/x" },
	    { "Error", "min(1.09878*(1+log(x)*(-0.47741+-0.000598835*x*x)+0.00625027*x*x)/x/x,0.200)" },
	   }
	 },
        }
      },
     }
   },
  };
