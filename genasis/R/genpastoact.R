genpastoact<-function(x,date_start=NA,date_end=NA,input="openair",output=NA,method="auto",days=28,pollutant=NA,temp=NA,wind=NA) {

  gaps<-function(volumepsm,film,density,surface,tsp,organic,particler,samplingr,temperature,aa,bb,cc,dd,rrt,time) density*10^(0.6366*(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature)))-3.1744)*volumepsm*(1-exp(-time*samplingr/(surface*density*10^(0.6366*(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature)))-3.1744)*film)))*(1-10^(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature))+log10(organic)-11.91)*tsp/(10^(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature))+log10(organic)-11.91)*tsp+1)+10^(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature))+log10(organic)-11.91)*tsp*particler/(10^(aa+bb/(273+temperature)+rrt*(cc+dd/(273+temperature))+log10(organic)-11.91)*tsp+1))/time
    
  # Constants of GENASIS model.
  com<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"acenaphtene","acenaphtylene",NA,"alpha-HCH","anthracene","benzo(a)anthracene","benzo(a)pyrene","benzo(b)fluoranthene","benzo(ghi)perylene","benzo(k)fluoranthene","beta-HCH",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"fluoranthene","fluorene","gamma-HCH",NA,"HCB","chrysene","indeno(123cd)pyrene","naphthalene",NA,NA,NA,NA,"pp-DDD","pp-DDE","pp-DDT",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PCB 118",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PCB 101",NA,NA,NA,NA,NA,NA,NA,NA,"PCB 138",NA,NA,NA,NA,"PCB 153",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PCB 180",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PCB 28",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PCB 52",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"PeCB",NA,"phenantrene","pyrene",NA,NA,NA)
  coc<-c("","1234678HpCDD","1234678HpCDF","","123478HxCDD","123478HxCDF","1234789HpCDF","","123678HxCDD","123678HxCDF","","12378PeCDD","123789HxCDD","123789HxCDF","12378PeCDF","","","","234678HxCDF","23478PeCDF","","2378TCDD","2378TCDF","","","","","ACE","ACY","ALD","A_HCH","ANT","BAA","BAP","BBFLA","BGP","BKFLA","B_HCH","","","","","","","","","","","","","","","","","","","","","","","","","A_CHL","CIS_NONA","D_HCH","DBAHA","","DIE","A_ENDS","END","FLA","FLU","G_HCH","HEP","HCB","CHRY","IP","NAP","OP_DDT","OCDD","OCDF","OXCHL","PP_DDD","PP_DDE","PP_DDT","PBDE100","PBDE126","PBDE153","PBDE154","","PBDE17","PBDE183","PBDE28","PBDE47","PBDE66","PBDE77","PBDE85","PBDE99","PCB105","PCB114","PCB118","","PCB16_32","","PCB189","","","","","","","PCB101","PCB110","PCB128","","","","","","PCB137","PCB138","PCB141","","PCB149","PCB151","PCB153","PCB156","","PCB17","PCB170","PCB171","PCB174","","PCB177","","","PCB18","PCB180","PCB183","PCB185","PCB187","","","","","PCB194","PCB195","","","","PCB200","PCB205","PCB206","PCB207","","PCB209","","","","PCB28","","","PCB31","PCB33","","","","PCB42","PCB44","","","PCB47","","PCB49","PCB52","PCB56_60","","","PCB66","","PCB74","","","","PCB87","","PCB95","","PCB99","","","","","","","","","PCN50","","PCN52","","","","","","","","","","PCN69","","PECB","PER","PHE","PYR","RET","G_CHL","TRANS_NONA")
  ooo<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,6.803602834,1.126956976,NA,6.75576367,4.246188594,0.913565156,0.291865455,0.583312601,0.515647932,0.827066021,5.807453416,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,5.83851886,5.156868718,8.056294483,NA,5.587025685,1.665274297,0.421454386,8.53087735,NA,NA,NA,NA,10.44528388,9.729836268,6.219969211,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,13.18681319,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,7.217533937,NA,NA,NA,NA,NA,NA,NA,NA,9.529161835,NA,NA,NA,NA,8.99534983,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,12.02678745,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,6.917050691,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,6.910007316,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,7.746833482,NA,6.234607907,5.041244844,NA,NA,NA)
  too<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.741534406,-0.568389589,NA,2.037530138,1.438764359,-0.752041898,-1.591200358,-0.824783513,-0.877126571,-0.806302973,1.881805571,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.859286176,1.72859837,2.142023552,NA,1.932658516,-0.153294648,-1.104989678,2.230207438,NA,NA,NA,NA,2.229418313,2.406941437,1.834503441,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.700257785,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.32674389,NA,NA,NA,NA,NA,NA,NA,NA,2.455395218,NA,NA,NA,NA,2.431104258,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.559440658,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.228671419,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.0503888,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.974052878,NA,2.213984682,1.463321468,NA,NA,NA)
  tot<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.013834856,0.110074364,NA,-0.024876003,0.027327269,0.11527873,0.090675599,0.053973311,0.025201713,0.070641131,-0.014975085,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.006309686,-0.011352595,-0.020710324,NA,-0.023493341,0.076519079,0.02636289,-0.016529517,NA,NA,NA,NA,0.02421464,-0.019837283,-0.001468515,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.024768629,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.038905422,NA,NA,NA,NA,NA,NA,NA,NA,-0.012312868,NA,NA,NA,NA,-0.020244097,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.002339142,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.034378321,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.017119587,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.005648169,NA,-0.033691043,0.012008167,NA,NA,NA)
  two<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.767560703,-0.258588571,NA,1.389518744,1.043745548,-1.565961394,-1.541150838,-1.227208148,-1.182624766,-1.527653921,1.216074957,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.070960995,1.43223736,1.02706243,NA,2.087042827,-1.387439282,-2.347207219,2.642350323,NA,NA,NA,NA,1.527616444,1.930475344,1.319842658,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.790655768,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.943720098,NA,NA,NA,NA,NA,NA,NA,NA,2.040641435,NA,NA,NA,NA,2.310989828,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.5619772,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,1.948091489,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.965204231,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,2.030074338,NA,1.581230182,0.633079489,NA,NA,NA)
  twt<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.013521929,0.106592701,NA,-0.016793449,0.032341818,0.124973491,0.090057582,0.059677329,0.028305949,0.079340282,-0.003833491,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.003269229,-0.007751521,-0.006228672,NA,-0.02536926,0.090908678,0.036456229,-0.021537453,NA,NA,NA,NA,0.030090919,-0.014138407,0.007733518,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.006698611,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.020639028,NA,NA,NA,NA,NA,NA,NA,NA,-0.004420616,NA,NA,NA,NA,-0.018690052,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.020780451,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.030947348,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.002996872,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.004950156,NA,-0.026002463,0.022096412,NA,NA,NA)
  tww<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.006354579,-0.07567419,NA,0.157226059,0.096359125,0.197322354,-0.012407391,0.09665498,0.077427859,0.17503776,0.171941662,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.192422643,0.07233889,0.26972906,NA,-0.037683735,0.302175905,0.31578834,-0.100600139,NA,NA,NA,NA,0.174276559,0.116367725,0.123116144,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.212683942,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.331070081,NA,NA,NA,NA,NA,NA,NA,NA,0.095093767,NA,NA,NA,NA,0.029336877,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.218646952,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.068473224,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.26499566,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,-0.013746772,NA,0.154449326,0.202654132,NA,NA,NA)
  
  # Constants of GAPS model.
  aa<-c(-2.98,-2.98,-2.98,0.55,0.55,0.55,-2.98,0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55,-2.98,0.55,-2.98,0.55,0.55,0.55,0.55,0.55,-2.98,-2.98,-2.98,-2.98,-2.3,-2.22,-4.37,-3.23,-3.28,-5.32,-6.3,-6.3,-7.28,-6.3,-7.69,-5.91,-4.11,-4.84,-3.37,-1.82,-3.86,-6.48,-4.68,-5.42,-3.94,-2.39,-4.43,-7.05,-5.25,-6.07,-4.52,-2.96,-5.01,-7.63,-5.83,-6.73,-5.09,-3.54,-5.58,-8.29,-7.42,-7.45,-7.36,-3.53,-3.82,-5.9,-6.75,-4.26,-2.79,-3.61,-3.95,-2.39,-5.32,-7.28,-1.24,-5.95,-2.98,-2.98,-5.64,-3.94,-7.49,-5.63,-7.18,-8.41,-5.39,-4.62,-5.8,-3.45,-3.71,-3.54,-6.47,-7.88,-5.69,-6.22,-4.64,20.478,20.478,20.478,-2.0296,-2.0296,-2.0296,20.478,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,20.478,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-2.0296,-5.3894,-4.9793,-7.0795,-5.3439,-7.0234,-5.6245,-5.1902,-5.4892,-6.6258,-7.7754,-6.0162,-6.5852,-6.3517,-6.3365,-6.0457,-5.8767,-5.5495,-7.8384,-6.5377,-7.0871,-8.2035,-7.826,NA,-6.3,-3.28,-4.26,-5.01,-8.03,-9.17)
  bb<-c(1672,1672,1672,986,986,986,1672,986,986,986,986,986,986,986,986,1672,986,1672,986,986,986,986,986,1672,1672,1672,1672,2609,2555,3709,3231,3259,4613,5263,5263,5912,5263,4937,5001,3808,4296,3320,2290,3645,5380,4188,4675,3700,2670,4025,5760,4568,5110,4080,3049,4405,6140,4947,5544,4460,3429,4784,5127,5092,4856,5967,3428,3790,4333,4436,3909,2934,3415,3455,2914,4613,5912,1905,4590,1672,1672,4179,4185,5116,4603,5459,6077,5131,4931,5298,3803,4672,3889,5068,5576,4936,5331,4757,-5194,-5194,-5194,1310.7,1310.7,1310.7,-5194,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,-5194,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,1310.7,3786.4,3650.6,4684.6,4000.3,4682.3,4180.6,3908.4,4058.5,4629,5044.4,4394.4,4684.3,4586.2,4614.8,4417,4452.2,4666.2,5271.4,5039.6,5002.9,5395.9,5293.6,NA,5263,3259,3909,4405,5036,5501)
  cc<-c(7e-05,7e-05,7e-05,-0.0032,-0.0032,-0.0032,7e-05,-0.0032,-0.0032,-0.0032,-0.0032,-0.0032,-0.0032,-0.0032,-0.0032,7e-05,-0.0032,7e-05,-0.0032,-0.0032,-0.0032,-0.0032,-0.0032,7e-05,7e-05,7e-05,7e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7e-05,7e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-49.35,-49.35,-49.35,-5.5305,-5.5305,-5.5305,-49.35,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-49.35,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,-5.5305,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0)
  dd<-c(0.857,0.857,0.857,1.714,1.714,1.714,0.857,1.714,1.714,1.714,1.714,1.714,1.714,1.714,1.714,0.857,1.714,0.857,1.714,1.714,1.714,1.714,1.714,0.857,0.857,0.857,0.857,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.857,0.857,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18678,18678,18678,5879.8,5879.8,5879.8,18678,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,18678,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,5879.8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0)
  rrt<-c(2379,2994,2898,2573,2781,2708,2986,2780,2788,2714,2382,2587,2802,2772,2593,2152,2322,2290,2748,2551,2181,2386,2338,1993,1985,1985,1944,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3196,3147,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.5267,0.5178,0.5114,0.5061,0.3787,0.5753,0.6106,0.5725,0.3761,0.4585,0.3406,0.4735,0.2668,0.4746,0.4961,0.5559,0.5371,0.5192,0.5262,0.5165,0.4945,0.5358,0.5404,0.5318,0.5212,0.5104,0.5036,0.5245,0.5693,0.541,0.3704,0.5951,0.5688,0.5631,0.5484,0.5662,0.5454,0.533,0.369,0.5787,0.5529,0.5589,0.55,0.3538,0.5958,0.5825,0.5802,0.6297,0.6191,0.6029,0.598,0.5997,0.5854,0.6328,0.6497,0.6227,0.6184,0.666,0.4102,0.3947,0.3933,0.3993,0.296,0.3598,0.3985,0.4054,0.3114,0.4465,0.4412,0.4361,0.4339,0.4134,0.4183,0.4258,0.4264,0.424,0.4213,0.4704,0.3362,0.4415,0.4603,0.3291,0.4558,0.5027,0.4842,0.4928,0.4903,0.4645,0.4607,0.4872,0.4778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NA,0,0,0,0,0,0)
  volumepsm<-0.000210
  film<-0.00567
  density<-21000
  surface<-0.0370
  tsp<-25
  organic<-0.45
  particler<-1
  
  # Twinning of length to cover both names and abbreviations (IDs).
  com<-c(com,coc)
  ooo<-c(ooo,ooo)
  too<-c(too,too)
  tot<-c(tot,tot)
  two<-c(two,two)
  twt<-c(twt,twt)
  tww<-c(tww,tww)
  aa<-c(aa,aa)
  bb<-c(bb,bb)
  cc<-c(cc,cc)
  dd<-c(dd,dd)
  rrt<-c(rrt,rrt)

  ## Selection of relevant compound(s).
  # If no pollutant is specified
  if (is.na(max(pollutant))) {
    # If x is defined as genasis-type data frame, load the compounds.
    if (class(x)=="data.frame"&input=="genasis") {
      compounds<-as.character(unique(x[,2]))
    } 
    
    # If x is defined as openair-type data frame, load the compounds.
    if (class(x)=="data.frame"&input=="openair") {
      compounds<-colnames(x)[-which(is.element(colnames(x),c("date","date_end","temp","wind","note")))]
    } 
    
    # If x is vector, load the compounds.
    if (class(x)!="data.frame") {
      compounds<-pollutant
    }
  } else {
    compounds<-pollutant
  }
  
  
  ## Verification of compounds.
  if (length(compounds[which(!is.element(compounds,com))])>0) {
    stop(paste0(length(unique(compounds[which(!is.element(compounds,com))]))," compounds were not recognized (",paste(unique(compounds[which(!is.element(compounds,com))]), collapse=', '),")."))
  }
  compos<-compounds
  
  # Output type selection.
  if (is.na(output)&class(x)=="data.frame") {
    output<-input
  }

  
  ## Creating column variables.
  # If x is defined as genasis-type data frame, set column vectors.
  if (class(x)=="data.frame"&input=="genasis") {
    valu      <-as.numeric(x[which(is.element(x[,2],compos)),1])
    comp      <-as.character(x[which(is.element(x[,2],compos)),2])
    date_start<-gendate(x[which(is.element(x[,2],compos)),3])
    date_end  <-gendate(x[which(is.element(x[,2],compos)),4])
    if (ncol(x)>4) {
      temp    <-as.numeric(x[which(is.element(x[,2],compos)),5])
    } else {
      temp    <-rep(NA,length(x[which(is.element(x[,2],compos)),2]))
    }
    if (ncol(x)>5) {
      wind    <-as.numeric(x[which(is.element(x[,2],compos)),6])
    } else {
      wind    <-rep(NA,length(x[which(is.element(x[,2],compos)),2]))
    }
    note      <-rep(NA,length(x[which(is.element(x[,2],compos)),2]))
    
  }
  
  # If x is defined as openair-type data frame, set column vectors.  
  if (class(x)=="data.frame"&input=="openair") {
    valu      <-c()
    comp      <-c()
    date_start<-c()
    date_end  <-c()
    temp      <-c()
    wind      <-c()
    for (compound in compos) {
      valu      <-c(valu,as.numeric(x[,compound]))
      comp      <-c(comp,as.character(rep(compound,nrow(x))))
      date_start<-as.Date(c(as.character(date_start),as.character(x[,"date"])))
      if (is.element("date_end",colnames(x))) {
        date_end<-as.Date(c(as.character(date_end),as.character(x[,"date_end"])))} else {
          date_end<-as.Date(c(as.character(date_end),as.character(x[,"date"])))
        }
      if (is.element("temp",colnames(x))) {
        temp<-c(temp,as.numeric(x[,"temp"]))} else {
          temp<-c(temp,rep(NA,nrow(x)))
        }
      if (is.element("wind",colnames(x))) {
        wind<-c(wind,as.numeric(x[,"wind"]))} else {
          wind<-c(wind,rep(NA,nrow(x)))
        }
    }
    note      <-rep(NULL,length(x))
  }
  
  # If x is defined as numeric vector, set column vectors.
  if (class(x)!="data.frame") {
    valu      <-x
    comp      <-rep(pollutant,length(x))
    date_start<-date_start
    date_end  <-date_end
    suppressWarnings(if (length(temp)==1&is.na(temp)) {
      temp    <-rep(NA,length(x))
    })
    suppressWarnings(if (length(wind)==1&is.na(wind)) { 
      wind    <-rep(NA,length(x))
    })
    note      <-rep(NULL,length(x))
    
    if (length(date_start[which(is.na(date_start))])==length(date_start)) {
      date_start<-rep(as.Date(Sys.Date()-days),length(x))
      note<-rep(paste0("Start/end date was set to actual (-",days," days), because not known."),length(x))
    }
    if(length(date_end)==1&is.na(date_end[1])) { # Set length of date_end to desired one.
      date_end<-rep(NA,length(date_start))
    }
    date_end<-as.Date(date_end)
    if(is.na(max(date_end))) { # If there are some NAs in date_end, will be replaced.
      date_end[which(is.na(date_end))]<-gendate(date_start[which(is.na(date_end))])+days 
    }
    
    if (length(temp)!=length(x)) {warning("The length of \"temp\" and \"x\" differ.")}
    if (length(wind)!=length(x)) {warning("The length of \"wind\" and \"x\" differ.")}
  }
  
  rnote<-note
  
  ## Passive to active recalculation.
  rvalu<-rep(NA,length(valu))
  for (i in 1:length(valu)) {
    
    period<-as.numeric(gendate(date_end[i])-gendate(date_start[i]))
    if (is.na(period)) {
      period<-days
    }
    
    # Different methods of recalculation:
    if (!is.na(method)&method=="const") {    
      rvalu[i]<-(valu[i]/period)/ooo[which(com==comp[i])[1]]
    }
    
    if (!is.na(method)&method=="temp") {    
      rvalu[i]<-(valu[i]/period)/(exp(too[which(com==comp[i])[1]]+temp[i]*tot[which(com==comp[i])[1]]))  
      }
    
    if (!is.na(method)&method=="tempwind") {    
      rvalu[i]<-(valu[i]/period)/(exp(two[which(com==comp[i])[1]]+temp[i]*twt[which(com==comp[i])[1]]+wind[i]*tww[which(com==comp[i])[1]]))    
    }

    if (is.na(method)|method=="auto") {    
      if (is.na(temp[i])) {
        rvalu[i]<-(valu[i]/period)/ooo[which(com==comp[i])[1]]
        rnote[i]<-paste0(note[i]," Multiplied by constant.")
      }
      if (!is.na(temp[i])&is.na(wind[i])) {
        rvalu[i]<-(valu[i]/period)/(exp(too[which(com==comp[i])[1]]+temp[i]*tot[which(com==comp[i])[1]]))
        rnote[i]<-paste0(note[i]," Temperature dependent.")
      }
      if (!is.na(temp[i])&!is.na(wind[i])) {
        rvalu[i]<-(valu[i]/period)/(exp(two[which(com==comp[i])[1]]+temp[i]*twt[which(com==comp[i])[1]]+wind[i]*tww[which(com==comp[i])[1]]))
        rnote[i]<-paste0(note[i]," Temperature and wind dependent.")
      }     
    }
    
    if (!is.na(method)&method=="gaps4") {
      rvalu[i]<-valu[i]/gaps(volumepsm,film,density,surface,tsp,organic,particler,4,temp[i],aa[which(com==comp[i])[1]],bb[which(com==comp[i])[1]],cc[which(com==comp[i])[1]],dd[which(com==comp[i])[1]],rrt[which(com==comp[i])[1]],period)/period
      rnote[i]<-paste0(note[i]," Gaps model with sampling rate 4.")
    }
    
    if (!is.na(method)&method=="gaps7") {
      rvalu[i]<-valu[i]/gaps(volumepsm,film,density,surface,tsp,organic,particler,7,temp[i],aa[which(com==comp[i])[1]],bb[which(com==comp[i])[1]],cc[which(com==comp[i])[1]],dd[which(com==comp[i])[1]],rrt[which(com==comp[i])[1]],period)/period
      rnote[i]<-paste0(note[i]," Gaps model with sampling rate 7.")
    }
    }

 
  ## Results assignment.
  rcomp<-comp
  rdate_start<-date_start
  rdate_end<-date_end
  rtemp<-temp
  rwind<-wind
  
  
  ## Generation of desired data type.
  # Output in genasis mode.
  if (!is.na(output)&output=="genasis") {
    result<-data.frame(rvalu,rcomp,as.Date(rdate_start,origin="1970-01-01"),as.Date(rdate_end,origin="1970-01-01"),rtemp,rwind,rnote)
    colnames(result)<-c("valu","comp","date_start","date_end","temp","wind","note")
  }
  
  # Output in openair mode.
  if (!is.na(output)&output=="openair") {
    unidates<-as.data.frame(matrix(NA,0,2))
    for (a in unique(rdate_start)) {
      for (b in unique(rdate_end[which(rdate_start==a)])) {
        unidates<-rbind(unidates,as.Date(c(a,b),origin="1970-01-01"))
      }
    }
    
    unidates[,1]<-as.Date(unidates[,1],origin="1970-01-01")
    unidates[,2]<-as.Date(unidates[,2],origin="1970-01-01")
    
    result<-as.data.frame(matrix(NA,nrow(unidates),length(compos)+5))
    colnames(result)<-c("date","date_end","temp","wind","note",as.character(compos))
    
    
    for (i in 1:nrow(unidates)) {
      result[i,1]<-unidates[i,1]
      result[i,2]<-unidates[i,2]
      result[i,3]<-rtemp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2])][1]
      result[i,4]<-rwind[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2])][1]     
      if (paste(rcomp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&!is.na(rnote))], collapse=', ')=="") {
        result[i,5]<-NA
      } else {
        result[i,5]<-paste0("Model found for: ",paste(rcomp[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&!is.na(rnote))], collapse=', '),".")
      }
      for (compound in compos) {
        result[i,compound]<-rvalu[which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound)][1]
      }
      if (length(which(rdate_start==unidates[i,1]&rdate_end==unidates[i,2]&rcomp==compound))>1) {
        warning(paste0("There was more than 1 record with the same start and end date for ",compound,", only the first record was used."))
        
      }
    }
    result[,1]<-as.Date(result[,1],origin="1970-01-01")
    result[,2]<-as.Date(result[,2],origin="1970-01-01")
  }
  
  # Output in vector-mode.
  if (class(x)!="data.frame"&is.na(output)) {
    result<-rvalu 
  }  
  return(res=result)
}