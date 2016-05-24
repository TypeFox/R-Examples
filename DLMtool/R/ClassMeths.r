
# Define Classes and Set Methods for DLMtool package
# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# DLMtool Classes:
# DLM_data - An object for storing data for analysis using data-limited methods

# Stock - An operating model component that specifies the parameters of the 
#   population dynamics model
# Fleet - The component of the operating model that controls fishing dynamics
# Observation - An operating model component that controls the observation model
# OM - An object containing all the parameters needed to control the MSE which 
#   is  built from components: Stock, Fleet and Observation objects.

# MSE - A Management Strategy Evaluation object that contains information about 
#  simulation conditions and performance of data-limited methods

# DLM_fease - An object for storing information about what data are available 
#  or might be available

# DLM_general - An object for storing general toolkit data. The data are stored
#   in the right format in the slot data

# lmmodel - An object for storing fitted linear model objects in this case the 
#   relationship between M, age-at-maturity and the von B. K parameter.

# Create DLM_data class 
setClass("DLM_data", representation(Name = "character", Year = "vector", 
    Cat = "matrix", Ind = "matrix", Rec = "matrix", t = "vector", 
    AvC = "vector", Dt = "vector", Mort = "vector", FMSY_M = "vector", 
    BMSY_B0 = "vector", Cref = "vector", Bref = "vector", Iref = "vector", 
    L50 = "vector", L95 = "vector", LFC = "vector", LFS = "vector", 
    CAA = "array", Dep = "vector", Abun = "vector", vbK = "vector", 
    vbLinf = "vector", vbt0 = "vector", wla = "vector", wlb = "vector", 
    steep = "vector", CV_Cat = "vector", CV_Dt = "vector", CV_AvC = "vector", 
    CV_Ind = "vector", CV_Mort = "vector", CV_FMSY_M = "vector", 
    CV_BMSY_B0 = "vector", CV_Cref = "vector", CV_Bref = "vector", 
    CV_Iref = "vector", CV_Rec = "vector", CV_Dep = "vector", 
    CV_Abun = "vector", CV_vbK = "vector", CV_vbLinf = "vector", 
    CV_vbt0 = "vector", CV_L50 = "vector", CV_LFC = "vector", 
    CV_LFS = "vector", CV_wla = "vector", CV_wlb = "vector", 
    CV_steep = "vector", sigmaL = "vector", MaxAge = "vector", 
    Units = "character", Ref = "numeric", Ref_type = "character", 
    Log = "list", params = "list", PosMPs = "vector", MPs = "vector", 
    OM = "data.frame", Obs = "data.frame", TAC = "array", TACbias = "array", 
    Sense = "array", CAL_bins = "numeric", CAL = "array", MPrec = "vector", 
    ML = "array", Lbar = "array", Lc = "array", LHYear = "numeric",
	Misc = "list")) 

# initialize DLM_data
setMethod("initialize", "DLM_data", function(.Object,stock="nada"){
#.Object
  #}) 
  #.Object<-new('DLM_data')
  # run an error check here
  if(file.exists(stock)){
    dat <- read.csv(stock,header=F,colClasses="character") # read 1st sheet
    dname<-dat[,1]
    dat<-dat[,2:ncol(dat)]

    .Object@Name<-dat[match("Name", dname),1]
    .Object@Year <-as.numeric(dat[match("Year",dname),dat[match("Year",dname),]!=""])
    .Object@Cat <-matrix(as.numeric(dat[match("Catch",dname),dat[match("Catch",dname),]!=""]),nrow=1)
    .Object@Ind<-matrix(as.numeric(dat[match("Abundance index",dname),1:length(.Object@Year)]),nrow=1)
    .Object@Rec<-matrix(as.numeric(dat[match("Recruitment",dname),1:length(.Object@Year)]),nrow=1)
    .Object@t<-as.numeric(dat[match("Duration t",dname),1])
    .Object@AvC<-as.numeric(dat[match("Average catch over time t",dname),1])
    .Object@Dt<-as.numeric(dat[match("Depletion over time t",dname),1])
    .Object@Mort<-as.numeric(dat[match("M",dname),1])
    .Object@FMSY_M<-as.numeric(dat[match("FMSY/M",dname),1])
    .Object@BMSY_B0<-as.numeric(dat[match("BMSY/B0",dname),1])
    .Object@Cref<-as.numeric(dat[match("Cref",dname),1])
    .Object@Bref<-as.numeric(dat[match("Bref",dname),1])
    .Object@Iref<-as.numeric(dat[match("Iref",dname),1])
    .Object@L50<-as.numeric(dat[match("Length at 50% maturity",dname),1])
    .Object@L95<-as.numeric(dat[match("Length at 95% maturity",dname),1])
    .Object@LFC<-as.numeric(dat[match("Length at first capture",dname),1])
    .Object@LFS<-as.numeric(dat[match("Length at full selection",dname),1])
    .Object@Dep<-as.numeric(dat[match("Current stock depletion",dname),1])
    .Object@Abun<-as.numeric(dat[match("Current stock abundance",dname),1])
    .Object@vbK<-as.numeric(dat[match("Von Bertalanffy K parameter", dname),1])
    .Object@vbLinf<-as.numeric(dat[match("Von Bertalanffy Linf parameter", dname),1])
    .Object@vbt0<-as.numeric(dat[match("Von Bertalanffy t0 parameter", dname),1])
    .Object@wla<-as.numeric(dat[match("Length-weight parameter a", dname),1])
    .Object@wlb<-as.numeric(dat[match("Length-weight parameter b", dname),1])
    .Object@steep<-as.numeric(dat[match("Steepness", dname),1])
    .Object@sigmaL<-as.numeric(dat[match("Sigma length composition", dname),1])
  
    .Object@CV_Cat<-as.numeric(dat[match("CV Catch", dname),1])
    .Object@CV_Dt<-as.numeric(dat[match("CV Depletion over time t", dname),1])
    .Object@CV_AvC<-as.numeric(dat[match("CV Average catch over time t", dname),1])
    .Object@CV_Ind<-as.numeric(dat[match("CV Abundance index", dname),1])
    .Object@CV_Mort<-as.numeric(dat[match("CV M", dname),1])
    .Object@CV_Rec<-as.numeric(dat[match("CV Rec", dname),1])
    .Object@CV_FMSY_M<-as.numeric(dat[match("CV FMSY/M", dname),1])
    .Object@CV_BMSY_B0<-as.numeric(dat[match("CV BMSY/B0", dname),1])
    .Object@CV_Cref<-as.numeric(dat[match("CV Cref", dname),1])
    .Object@CV_Bref<-as.numeric(dat[match("CV Bref", dname),1])
    .Object@CV_Iref<-as.numeric(dat[match("CV Iref", dname),1])
    .Object@CV_Dep<-as.numeric(dat[match("CV current stock depletion", dname),1])
    .Object@CV_Abun<-as.numeric(dat[match("CV current stock abundance", dname),1])
    .Object@CV_vbK<-as.numeric(dat[match("CV von B. K parameter", dname),1])
    .Object@CV_vbLinf<-as.numeric(dat[match("CV von B. Linf parameter", dname),1])
    .Object@CV_vbt0<-as.numeric(dat[match("CV von B. t0 parameter", dname),1])
    .Object@CV_L50<-as.numeric(dat[match("CV Length at 50% maturity", dname),1])
    .Object@CV_LFC<-as.numeric(dat[match("CV Length at first capture", dname),1])
    .Object@CV_LFS<-as.numeric(dat[match("CV Length at full selection", dname),1])
    .Object@CV_wla<-as.numeric(dat[match("CV Length-weight parameter a", dname),1])
    .Object@CV_wlb<-as.numeric(dat[match("CV Length-weight parameter b", dname),1])
    .Object@CV_steep<-as.numeric(dat[match("CV Steepness", dname),1])
    .Object@MaxAge<-as.numeric(dat[match("Maximum age", dname),1])
    .Object@MPrec<-as.numeric(dat[match("MPrec", dname),1])
    
    if(length(grep("CAL",dname))>1){
      CAL_bins<-as.numeric(dat[match("CAL_bins",dname),dat[match("CAL_bins",dname),]!=""])
      nCAL<-length(CAL_bins)-1
      .Object@CAL_bins<-CAL_bins
      CALdat<-grep("CAL ",dname)
      .Object@CAL<-array(as.numeric(as.matrix(dat[CALdat,1:nCAL])),dim=c(1,length(CALdat),nCAL))
    }

    CAAy<-grep("CAA",dname)[1:length(grep("CAA",dname))]
    CAAa<-sum(dat[CAAy[1],]!="")
    if(!is.na(CAAa)){
      .Object@CAA<-array(as.numeric(as.matrix(dat[CAAy,1:CAAa])),dim=c(1,length(CAAy),CAAa))
    }
    
    .Object@ML<-matrix(as.numeric(dat[match("Mean length",dname),1:length(.Object@Year)]),nrow=1)
    .Object@Lbar<-matrix(as.numeric(dat[match("Mean length Lc",dname),1:length(.Object@Year)]),nrow=1)
    .Object@Lc<-matrix(as.numeric(dat[match("Modal length",dname),1:length(.Object@Year)]),nrow=1)
    
    .Object@LHYear<-as.numeric(dat[match("LHYear",dname),1])
    .Object@Units<-dat[match("Units", dname),1]
    .Object@Ref<-as.numeric(dat[match("Reference TAC",dname),1])
    .Object@Ref_type<-dat[match("Reference TAC type",dname),1]
    .Object@Log[[1]]<-paste("Created:", Sys.time())
    .Object@params<-new('list')
    .Object@OM<-data.frame(NA)
    .Object@Obs<-data.frame(NA)
    .Object@TAC<-array(NA,dim=c(1,1,1))
    .Object@TACbias<-array(NA,dim=c(1,1,1))
    .Object@Sense<-array(NA,dim=c(1,1,1))
    .Object@PosMPs<-NA
    .Object@MPs<-NA

  }else{
    if(stock!="MSE"){
      if(!is.na(stock))print("Couldn't find specified csv file, blank DLM object created")
    }
  }
  # Default values -------------------------------------------------------------
  if(NAor0(.Object@CV_Cat)).Object@CV_Cat<-0.2
  if(NAor0(.Object@CV_Dt)).Object@CV_Dt<-0.25
  if(NAor0(.Object@CV_AvC)).Object@CV_AvC<-0.2
  if(NAor0(.Object@CV_Ind)).Object@CV_Ind<-0.2
  if(NAor0(.Object@CV_Mort)).Object@CV_Mort<-0.2
  if(NAor0(.Object@CV_FMSY_M)).Object@CV_FMSY_M<-0.2
  if(NAor0(.Object@CV_BMSY_B0)).Object@CV_BMSY_B0<-0.045
  if(NAor0(.Object@CV_Cref)).Object@CV_Cref<-0.2
  if(NAor0(.Object@CV_Bref)).Object@CV_Bref<-0.2
  if(NAor0(.Object@CV_Iref)).Object@CV_Iref<-0.2
  if(NAor0(.Object@CV_Rec)).Object@CV_Rec<-0.2
  if(NAor0(.Object@CV_Dep)).Object@CV_Dep<-0.25
  if(NAor0(.Object@CV_Abun)).Object@CV_Abun<-0.25
  if(NAor0(.Object@CV_vbK)).Object@CV_vbK<-0.1
  if(NAor0(.Object@CV_vbLinf)).Object@CV_vbLinf<-0.1
  if(NAor0(.Object@CV_vbt0)).Object@CV_vbt0<-0.1
  if(NAor0(.Object@CV_L50)).Object@CV_L50<-0.1
  if(NAor0(.Object@CV_LFC)).Object@CV_LFC<-0.2
  if(NAor0(.Object@CV_LFS)).Object@CV_LFS<-0.2
  if(NAor0(.Object@CV_wla)).Object@CV_wla<-0.1
  if(NAor0(.Object@CV_wlb)).Object@CV_wlb<-0.1
  if(NAor0(.Object@CV_steep)).Object@CV_steep<-0.2
  if(length(.Object@sigmaL)==0).Object@sigmaL<-0.2
  if(length(.Object@CAA)==0).Object@CAA<-array(NA,c(1,1,1))
  if(length(.Object@CAL)==0).Object@CAL<-array(NA,c(1,1,1))
  if(length(.Object@CAL_bins)==0).Object@CAL_bins<-1
  if(length(.Object@TAC)==0).Object@TAC<-array(1,c(1,1))
  if(length(.Object@TACbias)==0).Object@TACbias<-array(1,c(1,1))
  if(length(.Object@Sense)==0).Object@Sense<-array(1,c(1,1))
  if(length(.Object@ML)==0).Object@ML<-array(NA,c(1,1))
  if(length(.Object@Lbar)==0).Object@Lbar<-array(NA,c(1,1))
  if(length(.Object@Lc)==0).Object@Lc<-array(NA,c(1,1))
  
  
  .Object
})

# Create Stock class
setClass("Stock",representation(Name="character",maxage="numeric",R0="numeric", 
  M="numeric", Msd="numeric",Mgrad="numeric",h="numeric",SRrel="numeric",  
  Linf="numeric",K="numeric",t0="numeric", Ksd="numeric",Kgrad="numeric", 
  Linfsd="numeric",Linfgrad="numeric",recgrad="numeric", a="numeric", 
  b="numeric",D="numeric",Perr="numeric", Period="numeric", Amplitude="numeric",
  Size_area_1="numeric", 
  Frac_area_1="numeric", Prob_staying="numeric",AC="numeric", 
  L50="numeric", L50_95="numeric",Source="character"))
# initialize Stock
setMethod("initialize", "Stock", function(.Object,file=NA){

  if (!is.na(file)) {
    if (file.exists(file)) {
      dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      dname<-dat[,1]
      dat<-dat[,2:ncol(dat)]
      
      .Object@Name<-dat[match("Name", dname),1]
      .Object@maxage<-as.numeric(dat[match("maxage",dname),1])
      .Object@R0<-as.numeric(dat[match("R0",dname),1])
      .Object@M<-as.numeric(dat[match("M",dname),1:2])
      .Object@Msd<-as.numeric(dat[match("Msd",dname),1:2])
      .Object@Mgrad<-as.numeric(dat[match("Mgrad",dname),1:2])
      .Object@h<-as.numeric(dat[match("h",dname),1:2])
      .Object@SRrel<-as.numeric(dat[match("SRrel",dname),1])
      .Object@Linf<-as.numeric(dat[match("Linf",dname),1:2])
      .Object@K<-as.numeric(dat[match("K",dname),1:2])
      .Object@t0<-as.numeric(dat[match("t0",dname),1:2])
      .Object@Ksd<-as.numeric(dat[match("Ksd",dname),1:2])
      .Object@Kgrad<-as.numeric(dat[match("Kgrad",dname),1:2])
      .Object@Linfsd<-as.numeric(dat[match("Linfsd",dname),1:2])
      .Object@Linfgrad<-as.numeric(dat[match("Linfgrad",dname),1:2])
      .Object@recgrad<-as.numeric(dat[match("recgrad",dname),1:2])
      .Object@a<-as.numeric(dat[match("a",dname),1])
      .Object@b<-as.numeric(dat[match("b",dname),1])
      .Object@D<-as.numeric(dat[match("D",dname),1:2])
      .Object@Perr<-as.numeric(dat[match("Perr",dname),1:2])
	  .Object@Period<-as.numeric(dat[match("Period",dname),1:2])
	  .Object@Amplitude<-as.numeric(dat[match("Amplitude",dname),1:2])
      .Object@AC<-as.numeric(dat[match("AC",dname),1:2])
      .Object@Size_area_1<-as.numeric(dat[match("Size_area_1",dname),1:2])
      .Object@Frac_area_1<-as.numeric(dat[match("Frac_area_1",dname),1:2])
      .Object@Prob_staying<-as.numeric(dat[match("Prob_staying",dname),1:2])
      .Object@L50<-as.numeric(dat[match("L50", dname),1:2])
      .Object@L50_95<-as.numeric(dat[match("L50_95", dname),1:2])
      .Object@Source<-dat[match("Source", dname),1]
	} else {
	  message("File doesn't exist")
	}
   }	
  .Object

})

# Create Fleet class 
setClass("Fleet",slots=c(Name="character",nyears="numeric", Spat_targ="numeric",
  Fsd="numeric", qinc="numeric",qcv="numeric", 
  EffYears="numeric", EffLower="numeric", EffUpper="numeric",
  SelYears="numeric", AbsSelYears="numeric", L5="numeric", LFS="numeric",
  Vmaxlen="numeric", L5Lower="numeric", L5Upper="numeric", LFSLower="numeric",
  LFSUpper="numeric",  VmaxLower="numeric", VmaxUpper="numeric"))
# initialize Fleet 
setMethod("initialize", "Fleet", function(.Object,file=NA){
  if (!is.na(file)) {
    if (file.exists(file)) {
      dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      dname<-dat[,1]
      dat<-dat[,2:ncol(dat)]
      
      .Object@Name<-dat[match("Name", dname),1]
      .Object@nyears <-as.numeric(dat[match("nyears",dname),1])
      .Object@Spat_targ<-as.numeric(dat[match("Spat_targ",dname),1:2])
      .Object@Fsd<-as.numeric(dat[match("Fsd",dname),1:2])
      # .Object@Fgrad<-as.numeric(dat[match("Fgrad",dname),1:2])
      nEffYears <- ncol(dat[match("EffYears",dname),])
	  oldw <- getOption("warn")
	  options(warn=-1)
	  chk <- as.numeric(dat[match("EffYears",dname),1:nEffYears])
	  options(warn = oldw)
	  ind <- which(!is.na(chk))
	  nEffYears <- length(ind)
      .Object@EffYears <-as.numeric(dat[match("EffYears",dname),1:nEffYears])
      .Object@EffLower <-as.numeric(dat[match("EffLower",dname),1:nEffYears])
      .Object@EffUpper <-as.numeric(dat[match("EffUpper",dname),1:nEffYears])
      
      .Object@qinc<-as.numeric(dat[match("qinc",dname),1:2])
      .Object@qcv<-as.numeric(dat[match("qcv",dname),1:2])
	  
	  chkName <- match("SelYears",dname) # Check if vector of selectivity years exists
	  if (is.finite(chkName)) {
	    nSelYears <- ncol(dat[match("SelYears",dname),])
		oldw <- getOption("warn")
		options(warn=-1)
		chk <- as.numeric(dat[match("SelYears",dname),1:nSelYears])
		options(warn = oldw)
		ind <- which(is.finite(chk))
	    nSelYears <- length(ind)
	    chk <- length(ind)
	    if (is.finite(chk) &  chk > 0) { # parameters for selectivity years exists 
		  .Object@SelYears <- as.numeric(dat[match("SelYears",dname),1:nSelYears])
		  .Object@L5Lower <- as.numeric(dat[match("L5Lower",dname),1:nSelYears])
		  .Object@L5Upper <- as.numeric(dat[match("L5Upper",dname),1:nSelYears])
		  .Object@LFSLower <- as.numeric(dat[match("LFSLower",dname),1:nSelYears])
		  .Object@LFSUpper <- as.numeric(dat[match("LFSUpper",dname),1:nSelYears])
		  .Object@VmaxLower <- as.numeric(dat[match("VmaxLower",dname),1:nSelYears])
		  .Object@VmaxUpper <- as.numeric(dat[match("VmaxUpper",dname),1:nSelYears]) 
	    }
	  }
      # These are ignored in MSE if L5Lower etc are set 
	  .Object@L5 <- as.numeric(dat[match("L5",dname),1:2])
      .Object@LFS <- as.numeric(dat[match("LFS",dname),1:2])
      .Object@Vmaxlen <-as.numeric(dat[match("Vmaxlen",dname),1:2])
	  
    } else {
	  message("File doesn't exist")
	}
   }	
 	
  .Object
})


# Create Observation class 
setClass("Observation",representation(Name="character",LenMcv="numeric",
  Cobs="numeric",Cbiascv="numeric",CAA_nsamp="numeric",CAA_ESS="numeric",
  CAL_nsamp="numeric",CAL_ESS="numeric",CALcv="numeric",
  Iobs="numeric",Mcv="numeric",Kcv="numeric",t0cv="numeric",Linfcv="numeric",
  LFCcv="numeric",LFScv="numeric",B0cv="numeric",
  FMSYcv="numeric",FMSY_Mcv="numeric",BMSY_B0cv="numeric",
  rcv="numeric", Dbiascv="numeric",Dcv="numeric",
  Btbias="numeric",Btcv="numeric",Fcurbiascv="numeric",Fcurcv="numeric",
  hcv="numeric",Icv="numeric",maxagecv="numeric",Reccv="numeric",
  Irefcv="numeric",Crefcv="numeric",Brefcv="numeric",beta="numeric"))
# initialize Observation 
setMethod("initialize", "Observation", function(.Object,file=NA){
  if (!is.na(file)) {
    if (file.exists(file)) {
      dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
      dname<-dat[,1]
      dat<-dat[,2:ncol(dat)]
      
      .Object@Name<-dat[match("Name", dname),1]
      .Object@LenMcv<-as.numeric(dat[match("LenMcv",dname),1])
      .Object@Cobs<-as.numeric(dat[match("Cobs",dname),1:2])
      .Object@Cbiascv<-as.numeric(dat[match("Cbiascv",dname),1])
      .Object@CAA_nsamp<-as.numeric(dat[match("CAA_nsamp",dname),1:2])
      .Object@CAA_ESS<-as.numeric(dat[match("CAA_ESS",dname),1:2])
      .Object@CAL_nsamp<-as.numeric(dat[match("CAA_nsamp",dname),1:2])
      .Object@CAL_ESS<-as.numeric(dat[match("CAA_ESS",dname),1:2])
      .Object@CALcv<-as.numeric(dat[match("CALcv",dname),1:2])
      .Object@Iobs<-as.numeric(dat[match("Iobs",dname),1:2])
      .Object@Mcv<-as.numeric(dat[match("Mcv",dname),1])
      .Object@Kcv<-as.numeric(dat[match("Kcv",dname),1])
      .Object@t0cv<-as.numeric(dat[match("t0cv",dname),1])
      .Object@Linfcv<-as.numeric(dat[match("Linfcv",dname),1])
      .Object@LFCcv<-as.numeric(dat[match("LFCcv",dname),1])
      .Object@LFScv<-as.numeric(dat[match("LFScv",dname),1])
      .Object@B0cv<-as.numeric(dat[match("B0cv",dname),1])
      .Object@FMSYcv<-as.numeric(dat[match("FMSYcv",dname),1])
      .Object@FMSY_Mcv<-as.numeric(dat[match("FMSY_Mcv",dname),1])
      .Object@BMSY_B0cv<-as.numeric(dat[match("BMSY_B0cv",dname),1])
      .Object@rcv<-as.numeric(dat[match("rcv",dname),1])
      .Object@Dbiascv<-as.numeric(dat[match("Dbiascv",dname),1])
      .Object@Dcv<-as.numeric(dat[match("Dcv",dname),1:2])
      .Object@Btbias<-as.numeric(dat[match("Btbias",dname),1:2])
      .Object@Btcv<-as.numeric(dat[match("Btcv",dname),1:2])
      .Object@Fcurbiascv<-as.numeric(dat[match("Fcurbiascv",dname),1])
      .Object@Fcurcv<-as.numeric(dat[match("Fcurcv",dname),1:2])
      .Object@hcv<-as.numeric(dat[match("hcv",dname),1])
      .Object@Icv<-as.numeric(dat[match("Icv",dname),1])
      .Object@maxagecv<-as.numeric(dat[match("maxagecv",dname),1])
      .Object@Reccv<-as.numeric(dat[match("Reccv",dname),1:2])
      .Object@Irefcv<-as.numeric(dat[match("Irefcv",dname),1])
      .Object@Crefcv<-as.numeric(dat[match("Crefcv",dname),1])
      .Object@Brefcv<-as.numeric(dat[match("Brefcv",dname),1])
      .Object@beta<-as.numeric(dat[match("beta",dname),1:2])
	} else {
	  message("File doesn't exist")
	}
  }	
  .Object

})

# Create OM class
setClass("OM",representation(Name="character",nyears="numeric",maxage="numeric",
  R0="numeric",M="numeric", Msd="numeric",Mgrad="numeric",h="numeric", 
  SRrel="numeric",Linf="numeric",K="numeric",t0="numeric", Ksd="numeric", 
  Kgrad="numeric",Linfsd="numeric",Linfgrad="numeric",recgrad="numeric",
  a="numeric",b="numeric",D="numeric", Size_area_1="numeric", 
  Frac_area_1="numeric",Prob_staying="numeric", Source="character", 
  L50="numeric", L50_95="numeric", SelYears="numeric", AbsSelYears="numeric",
  L5="numeric", LFS="numeric",  Vmaxlen="numeric", 
  L5Lower="numeric", L5Upper="numeric", LFSLower="numeric",
  LFSUpper="numeric",  VmaxLower="numeric", VmaxUpper="numeric",
  beta="numeric", 
  Spat_targ="numeric", Fsd="numeric", Period="numeric", Amplitude="numeric",
  EffYears="numeric", EffLower="numeric", EffUpper="numeric", 
  # Fgrad="numeric", 
  qinc="numeric",qcv="numeric",AC="numeric", Cobs="numeric",Cbiascv="numeric",
  CAA_nsamp="numeric",CAA_ESS="numeric", CAL_nsamp="numeric",CAL_ESS="numeric", 
  CALcv="numeric", Iobs="numeric",Perr="numeric", Mcv="numeric",Kcv="numeric", 
  t0cv="numeric",Linfcv="numeric", LFCcv="numeric", LFScv="numeric", 
  B0cv="numeric",FMSYcv="numeric",FMSY_Mcv="numeric",BMSY_B0cv="numeric",
  LenMcv="numeric",rcv="numeric", Dbiascv="numeric",Dcv="numeric", 
  Btbias="numeric",Btcv="numeric", Fcurbiascv="numeric",Fcurcv="numeric", 
  hcv="numeric", Icv="numeric",maxagecv="numeric", Reccv="numeric", 
  Irefcv="numeric",Crefcv="numeric",Brefcv="numeric"))
# initialize OM 
setMethod("initialize", "OM", function(.Object,Stock,Fleet,Observation){
  if(class(Stock)!='Stock')print(paste('Could not build operating model:',deparse(substitute(Stock)),'not of class Stock'))
  if(class(Fleet)!='Fleet')print(paste('Could not build operating model:',deparse(substitute(Fleet)),'not of class Fleet'))
  if(class(Observation)!='Observation')print(paste('Could not build operating model:',deparse(substitute(Observation)),'not of class Observation'))
  if(class(Stock)!='Stock'|class(Fleet)!='Fleet'|class(Observation)!='Observation')stop()
   
  .Object@Name<-paste("Stock:",Stock@Name,"  Fleet:",Fleet@Name,"  Observation model:",Observation@Name,sep="")
  # Now copy the values for stock, fleet and observation slots to same slots in the Sim object
  Sslots<-slotNames(Stock)
  for(i in 2:length(Sslots)) {
    tt <- .hasSlot(Stock,Sslots[i]) # For back-compatibility
    if (tt) slot(.Object,Sslots[i])<-slot(Stock,Sslots[i])
  }	
  Fslots<-slotNames(Fleet)
  for(i in 2:length(Fslots)) {
    tt <- .hasSlot(Fleet,Fslots[i])
	if (tt)  slot(.Object,Fslots[i])<-slot(Fleet,Fslots[i])
  }	
  Oslots<-slotNames(Observation)
  for(i in 2:length(Oslots)) {
    tt <- .hasSlot(Observation,Oslots[i])   
    if (tt) slot(.Object,Oslots[i])<-slot(Observation,Oslots[i])
  } 
  .Object
})

# Create MSE class 
setClass("MSE",representation(Name="character",nyears="numeric", 
  proyears="numeric",nMPs="numeric",MPs="character", nsim="numeric",
  OM="data.frame",Obs="data.frame",B_BMSY="array", F_FMSY="array",
  B="array",FM="array",C="array",TAC="array",SSB_hist="array",
  CB_hist="array",FM_hist="array"))

setMethod("initialize", "MSE", function(.Object,Name,nyears,proyears,nMPs,MPs,
                                                nsim,OMtable,Obs,B_BMSYa,F_FMSYa,Ba,FMa,Ca,TACa,SSB_hist,CB_hist,FM_hist){
  .Object@Name<-Name
  .Object@nyears <-nyears
  .Object@proyears<-proyears
  .Object@nMPs<-nMPs
  .Object@MPs<-MPs
  .Object@nsim<-nsim
  .Object@OM<-OMtable
  .Object@Obs<-Obs
  .Object@B_BMSY<-B_BMSYa
  .Object@F_FMSY<-F_FMSYa
  .Object@B<-Ba
  .Object@FM<-FMa
  .Object@C<-Ca
  .Object@TAC<-TACa
  .Object@SSB_hist<-SSB_hist
  .Object@CB_hist<-CB_hist
  .Object@FM_hist<-FM_hist
  .Object
})
  

# Create DLM_fease class 
setClass("DLM_fease",representation(Name="character",Case="character",Catch="numeric",
         Index="numeric",Natural_mortality_rate="numeric",Maturity_at_length="numeric",
         Growth="numeric",Length_weight_conversion="numeric",Fleet_selectivity="numeric",
         Catch_at_length="numeric",Catch_at_age="numeric",Recruitment_index="numeric",
         Stock_recruitment_relationship="numeric",Target_catch="numeric",Target_biomass="numeric",
         Target_index="numeric",Abundance="numeric"))
# initialize DLM_fease 
setMethod("initialize", "DLM_fease", function(.Object,file="nada",ncases=1){
  # run an error check here
  if(file.exists(file)){
    dat <- read.csv(file,header=F,colClasses="character") # read 1st sheet
    nr<-nrow(dat)
    ncases=ncol(dat)-1
    dname<-dat[,1]
    if(ncases==1)dat<-array(dat[,2:ncol(dat)],dim=c(nr,ncases))
    if(ncases>1)dat<-dat[,2:ncol(dat)]
    .Object@Name<- dat[match("Name", dname),1]
    .Object@Case<-as.character(dat[match("Case",dname),1:ncases])
    .Object@Catch<-as.numeric(dat[match("Catch",dname),1:ncases])
    .Object@Index<-as.numeric(dat[match("Index",dname),1:ncases])
    .Object@Natural_mortality_rate<-as.numeric(dat[match("Natural_mortality_rate",dname),1:ncases])
    .Object@Maturity_at_length<-as.numeric(dat[match("Maturity_at_length",dname),1:ncases])
    .Object@Growth<-as.numeric(dat[match("Growth",dname),1:ncases])
    .Object@Length_weight_conversion<-as.numeric(dat[match("Length_weight_conversion",dname),1:ncases])
    .Object@Fleet_selectivity<-as.numeric(dat[match("Fleet_selectivity",dname),1:ncases])
    .Object@Catch_at_length<-as.numeric(dat[match("Catch_at_length",dname),1:ncases])
    .Object@Catch_at_age<-as.numeric(dat[match("Catch_at_age",dname),1:ncases])
    .Object@Recruitment_index<-as.numeric(dat[match("Recruitment_index",dname),1:ncases])
    .Object@Stock_recruitment_relationship<-as.numeric(dat[match("Stock_recruitment_relationship",dname),1:ncases])
    .Object@Target_catch<-as.numeric(dat[match("Target_catch",dname),1:ncases])
    .Object@Target_biomass<-as.numeric(dat[match("Target_biomass",dname),1:ncases])
    .Object@Target_index<-as.numeric(dat[match("Target_index",dname),1:ncases])
    .Object@Abundance<-as.numeric(dat[match("Abundance",dname),1:ncases])
  }else{
    .Object@Name<-"Blank DLM_Fease"
    .Object@Case<-"Case 1"
    .Object@Catch<-0
    .Object@Index<-0
    .Object@Natural_mortality_rate<-0
    .Object@Maturity_at_length<-0
    .Object@Growth<-0
    .Object@Length_weight_conversion<-0
    .Object@Fleet_selectivity<-0
    .Object@Catch_at_length<-0
    .Object@Catch_at_age<-0
    .Object@Recruitment_index<-0
    .Object@Stock_recruitment_relationship<-0
    .Object@Target_catch<-0
    .Object@Target_biomass<-0
    .Object@Target_index<-0
    .Object@Abundance<-0
  }  
  .Object
  
})


# Create DLM_general object
setClass("DLM_general",representation(Name="character",data="list"))
# initialize DLM_general
setMethod("initialize", "DLM_general", function(.Object){
  .Object
})

# Create lmmodel class
setClass("lmmodel",representation(Name="character",models="list"))
# initialize lmmodel 
setMethod("initialize", "lmmodel", function(.Object,Name,models){
  .Object@Name<-Name
  .Object@models<-models
  .Object
})

# Define generic plot method for DLM data objects
setMethod("plot",
  signature(x = "DLM_data"),
  function(x,funcs=NA,maxlines=6,perc=0.5,xlims=NA){
    
    DLM_data<-x
    cols<-rep(c('black','red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey'),4)
    ltys<-rep(1:4,each=13)
    
    if(is.na(funcs[1]))funcs<-DLM_data@MPs

    nMPs<-length(funcs)
    nplots<-ceiling(nMPs/maxlines)
    maxl<-ceiling(nMPs/nplots)
    mbyp <- split(1:nMPs, ceiling(1:nMPs/maxl))   # assign methods to plots

    if(is.na(xlims[1])|length(xlims)!=2){
      xlims<-quantile(DLM_data@TAC,c(0.005,0.95),na.rm=T)
      if(xlims[1]<0)xlims[1]<-0
    }
    if(!NAor0(DLM_data@Ref)){
      if(xlims[1]>DLM_data@Ref)xlims[1]<-max(0,0.98*DLM_data@Ref)
      if(xlims[2]<DLM_data@Ref)xlims[2]<-1.02*DLM_data@Ref
    }
    ylims<-c(0,1)

    #for(m in 1:nMPs){
     # if(sum(!is.na(DLM_data@TAC[m,,1]))>2){
       # dens<-density(DLM_data@TAC[m,,1],na.rm=T)
        #print(quantile(dens$y,0.99,na.rm=T))
      #  if(quantile(dens$y,0.9,na.rm=T)>ylims[2])ylims[2]<-quantile(dens$y,0.90,na.rm=T)
      #}
    #}

    #dev.new2(width=10,height=0.5+7*nplots)
    par(mfrow=c(ceiling(nplots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.35,0.35,0.35,0.05))

    for(p in 1:nplots){
      m<-mbyp[[p]][1]
      plot(NA,NA,xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col="white",lwd=3,type="l")
      abline(h=0)
      if(!NAor0(DLM_data@Ref)){
        abline(v=DLM_data@Ref,col="light grey",lwd=2)
        if(!NAor0(DLM_data@Ref_type[1]))legend('right',DLM_data@Ref_type,text.col="grey",bty='n')
      }
      #plot(density(DLM@TAC[m,,1],from=0,na.rm=T),xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col=coly[m],lty=ltyy[m],type="l")

      if(!is.na(perc[1]))abline(v=quantile(DLM_data@TAC[m,,1],p=perc,na.rm=T),col=cols[m],lty=ltys[m])
      #if(length(mbyp[[p]])>0){
        for(ll in 1:length(mbyp[[p]])){
          m<-mbyp[[p]][ll]
          if(sum(!is.na(DLM_data@TAC[m,,1]))>10){  # only plot if there are sufficient non-NA TAC samples
            x<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$x
            y<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$y
            y<-y/max(y)
            lines(x,y,col=cols[ll])
          }else{
            print(paste("Method ",funcs[m]," produced too many NA TAC values for plotting densities",sep=""))
          }
          if(!is.na(perc[1]))abline(v=quantile(DLM_data@TAC[m,,1],p=perc,na.rm=T),col=cols[ll],lty=2)
        }
      #}
      cind<-1:length(mbyp[[p]])
      legend('topright',funcs[mbyp[[p]]],text.col=cols[cind],col=cols[cind],lty=1,bty='n',cex=0.75)
    }

    mtext(paste("TAC (",DLM_data@Units,")",sep=""),1,outer=T,line=0.5)
    mtext(paste("Standardized relative frequency",sep=""),2,outer=T,line=0.5)
    mtext(paste("TAC calculation for ",DLM_data@Name,sep=""),3,outer=T,line=0.5)
})

# Define generic summary method for DLM data objects
setMethod("summary",
          signature(object = "DLM_data"),
          function(object){
  
  scols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
  
  #dev.new2(width=8,height=4.5)
  par(mai=c(0.35,0.9,0.2,0.01),c(0.3,0,0,0))
  layout(matrix(c(1,2,1,2,1,2,3,3,3,3),nrow=2))
  plot(object@Year,object@Cat[1,],col="blue",type="l",xlab="Year",ylab=paste("Catch (",object@Units,")",sep=""),ylim=c(0,max(object@Cat[1,],na.rm=T)))
  plot(object@Year,object@Ind[1,],col="orange",type="l",xlab="Year",ylab="Relative abundance",ylim=c(0,max(object@Ind[1,],na.rm=T)))
  
  slots<-c("Dep","Mort","FMSY_M","Dt","BMSY_B0","vbK")
  namey<-c("Stock depletion", "Natural Mortality rate","Ratio of FMSY to M","Depletion over time t","BMSY relative to unfished","Von B. k parameter")
  slotsCV<-c("CV_Dep","CV_Mort","CV_FMSY_M","CV_Dt","CV_BMSY_B0","CV_vbK")
  
  ind<-rep(TRUE,length(slotsCV))
  for(i in 1:length(slotsCV))if(NAor0(attr(object,slots[i]))|NAor0(attr(object,slotsCV[i])))ind[i]<-FALSE
  slots<-slots[ind]
  slotsCV<-slotsCV[ind]
  nrep<-150
  xstore<-array(NA,c(length(slots),nrep))
  ystore<-array(NA,c(length(slots),nrep))
  
 
  for(i in 1:length(slots)){
    mu<-attr(object,slots[i])
    cv<-attr(object,slotsCV[i])
    xstore[i,]<-qlnorm(seq(0,1,length.out=nrep),mconv(mu,cv),sdconv(mu,cv))
    ystore[i,]<-dlnorm(xstore[i,],mconv(mu,cv),sdconv(mu,cv))
  }
  
  plot(xstore[1,],ystore[1,],type="l",xlim=c(0,1.2),ylim=c(0,quantile(ystore,0.97)),xlab="",ylab="Relative frequency",col=scols[1])
  if(length(slots)>1){
    for(i in 2:length(slots)) lines(xstore[i,],ystore[i,],col=scols[i])
  }
  legend('topright',legend=namey[ind],text.col=scols[1:length(slots)],bty='n')
  mtext(paste("Data summary for",deparse(substitute(DLM_data)),sep=" "),3,font=2,line=0.25,outer=T)

})

# Define generic summary method for MSE objects 
setMethod("summary",
          signature(object = "MSE"),
          function(object){            

    MSEobj<-object      
    nm<-MSEobj@nMPs
    nsim<-MSEobj@nsim
    proyears<-MSEobj@proyears
    
    Yd<-P10<-P50<-P100<-POF<-LTY<-STY<-VY<-array(NA,c(nm,nsim))
    
    yind<-max(MSEobj@proyears-4,1):MSEobj@proyears
    RefYd<-MSEobj@OM$RefY
    yend<-max(MSEobj@proyears-9,1):MSEobj@proyears
    ystart<-1:10
    y1<-1:(MSEobj@proyears-1)
    y2<-2:MSEobj@proyears
    
    for(m in 1:nm){
      Yd[m,]<-round(apply(MSEobj@C[,m,yind],1,mean,na.rm=T)/RefYd*100,1)
      POF[m,]<-round(apply(MSEobj@F_FMSY[,m,]>1,1,sum,na.rm=T)/proyears*100,1)
      P10[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.1,1,sum,na.rm=T)/proyears*100,1)
      P50[m,]<-round(apply(MSEobj@B_BMSY[,m,]<0.5,1,sum,na.rm=T)/proyears*100,1)
      P100[m,]<-round(apply(MSEobj@B_BMSY[,m,]<1,1,sum,na.rm=T)/proyears*100,1)
      LTY[m]<-round(sum(MSEobj@C[,m,yend]/RefYd>0.5)/(MSEobj@nsim*length(yend))*100,1)
      STY[m]<-round(sum(MSEobj@C[,m,ystart]/RefYd>0.5)/(MSEobj@nsim*length(ystart))*100,1)
      AAVY<-apply(((MSEobj@C[,m,y1]-MSEobj@C[,m,y2])^2)^0.5,1,mean)/apply(MSEobj@C[,m,y2],1,mean)
      VY[m]<-round(sum(AAVY<0.1)/MSEobj@nsim*100,1)
    }
    nr<-2
    out<-cbind(MSEobj@MPs,round(apply(Yd,1,mean,na.rm=T),nr),round(apply(Yd,1,sd,na.rm=T),nr),
                             round(apply(POF,1,mean,na.rm=T),nr),round(apply(POF,1,sd,na.rm=T),nr),
                             round(apply(P10,1,mean,na.rm=T),nr),round(apply(P10,1,sd,na.rm=T),nr),
                             round(apply(P50,1,mean,na.rm=T),nr),round(apply(P50,1,sd,na.rm=T),nr),
                             round(apply(P100,1,mean,na.rm=T),nr),round(apply(P100,1,sd,na.rm=T),nr),
                             round(apply(LTY,1,mean,na.rm=T),nr),
                             round(apply(STY,1,mean,na.rm=T),nr),
                             round(apply(VY,1,mean,na.rm=T),nr))
    out<-as.data.frame(out)
    names(out)<-c("MP","Yield","stdev","POF","stdev ","P10","stdev",
                  "P50","stdev","P100","stdev","LTY","STY","VY")
    out[,1]<-as.character(out[,1])
    for(i in 2:ncol(out))out[,i]<-as.numeric(as.character(out[,i]))
    out
  })


# Plotting code for MSE object
setMethod("plot",
  signature(x = "MSE"),
  function(x){
  MSEobj<-x
  Pplot(MSEobj)
  Kplot(MSEobj)
  Tplot(MSEobj)
})


