fbp  <- function(input=NULL,output="Primary"){       
  .Deprecated(new="fbp", 
              package="cffdrs", 
              msg="The 'fwi.fbp' package and contained functions are being deprecated and replaced by the 'cffdrs' package, please update your code to use the 'cffdrs' package.",
              old="fbp")  
# The choice of output include: "Primarry" (default), "Secondary", and "All".
### detach the dataset when it is attached from last fail run
if(!is.na(charmatch("input",search()))) {detach(input)}
output<-toupper(output)
### Here is the main function:
if(is.null(input)){                                                                                                                     # default input data is NULL
    FUELTYPE     <- "C2"; ACCEL <- 0; DJ <- 180;D0 <- 0; ELV <- 0; BUIEFF <- 1; HR <- 1; FFMC <- 90;
    ISI          <- 0;BUI <- 60; WS <- 10; WD <- 0; GS  <- 0; ASPECT <- 0;PC <- 50; PDF <- 35; CC <- 80;
    GFL          <- 0.35; CBH <- 3; CFL <- 1; LAT <- 55; LONG <- -120; FMC <- 0; THETA <- 0
    input        <- as.data.frame(cbind(ACCEL,DJ,D0,ELV,BUIEFF,HR,FFMC,ISI,BUI,WS,WD,GS,ASPECT,PC,PDF,CC,GFL,CBH,CFL,LAT,LONG,FMC,THETA))
    input        <- cbind(FUELTYPE,input)
    input[,"FUELTYPE"]      <- as.character(input[,"FUELTYPE"])
    } else {
    names(input)<-toupper(names(input))
    ID<-input$ID;FUELTYPE<-toupper(input$FUELTYPE);FFMC<-input$FFMC;BUI<-input$BUI;WS<-input$WS;WD<-input$WD;FMC<-input$FMC;
    GS<-input$GS;LAT<-input$LAT;LONG<-input$LONG;ELV<-input$ELV;DJ<-input$DJ;D0<-input$D0;SD<-input$SD;SH<-input$SH;
    HR<-input$HR;PC<-input$PC;PDF<-input$PDF;GFL<-input$GFL;CC<-input$CC;THETA<-input$THETA;ACCEL<-input$ACCEL;
    ASPECT<-input$ASPECT;BUIEFF<-input$BUIEFF;CBH<-input$CBH;CFL<-input$CFL;ISI<-input$ISI
    n0 <- nrow(input)
    if(!exists("FUELTYPE")|is.null(FUELTYPE)) FUELTYPE<-rep("C2",n0);if(!exists("FFMC")|is.null(FFMC)) FFMC<-rep(90,n0);if(!exists("BUI")|is.null(BUI)) BUI<-rep(60,n0)
    if(!exists("WS")|is.null(WS)) WS<-rep(10,n0);if(!exists("WD")|is.null(WD)) WD<-rep(0,n0);if(!exists("FMC")|is.null(FMC)) FMC<-rep(0,n0);if(!exists("GS")|is.null(GS)) GS<-rep(0,n0)
    if(!exists("LAT")|is.null(LAT)) LAT<-rep(55,n0);if(!exists("LONG")|is.null(LONG)) LONG<-rep(-120,n0);if(!exists("ELV")|is.null(ELV)) ELV<-rep(0,n0)
    if(!exists("SD")|is.null(SD)) SD<-rep(0,n0);if(!exists("SH")|is.null(SH)) SH<-rep(0,n0);if(!exists("DJ")|is.null(DJ)) DJ<-rep(180,n0);if(!exists("D0")|is.null(D0)) D0<-rep(0,n0)
    if(!exists("HR")|is.null(HR)) HR<-rep(1,n0);if(!exists("PC")|is.null(PC)) PC<-rep(50,n0);if(!exists("PDF")|is.null(PDF)) PDF<-rep(35,n0);if(!exists("GFL")|is.null(GFL)) GFL<-rep(0.35,n0)
    if(!exists("CC")|is.null(CC)) CC<-rep(80,n0);if(!exists("THETA")|is.null(THETA)) THETA<-rep(0,n0);if(!exists("ACCEL")|is.null(ACCEL)) ACCEL<-rep(0,n0)
    if(!exists("ASPECT")|is.null(ASPECT)) ASPECT<-rep(0,n0);if(!exists("BUIEFF")|is.null(BUIEFF)) BUIEFF<-rep(1,n0);if(!exists("CBH")|is.null(CBH)) CBH<-rep(3,n0)
    if(!exists("CFL")|is.null(CFL)) CFL<-rep(1,n0);if(!exists("ISI")|is.null(ISI)) ISI<-rep(0,n0)
    
    ### Data cleaning up
    WD     <- WD * pi/180
    THETA  <- THETA * pi/180
    ASPECT <- ifelse(is.na(ASPECT),0,ASPECT)
    ASPECT <- ifelse(ASPECT < 0,ASPECT+360,ASPECT)
    ASPECT <- ASPECT * pi/180

    ACCEL  <- ifelse(is.na(ACCEL)|ACCEL < 0,0,ACCEL)                                              # line (no accelleration effect) */
    if (length(ACCEL[!ACCEL %in% c(0,1)])>0) warning("Input variable Accel is out of range, will be assigned to 1")
    ACCEL  <- ifelse(!ACCEL %in% c(0,1),1,ACCEL)
    
    DJ     <- ifelse(DJ < 0 | DJ > 366,0,DJ)
    DJ     <- ifelse(is.na(DJ),180,DJ)
    D0     <- ifelse(is.na(D0)|D0 < 0 | D0 > 366,0,D0)
    ELV    <- ifelse(ELV < 0 | ELV > 10000,0,ELV)
    ELV    <- ifelse(is.na(ELV),0,ELV)
    BUIEFF <- ifelse(BUIEFF < 0,0,1)
    BUIEFF <- ifelse(is.na(BUIEFF),1,BUIEFF)
    HR     <- ifelse(HR < 0,-HR,HR)                                                                                                    # Originally "T"
    HR     <- ifelse(HR > 366*24,24,HR)
    HR     <- ifelse(is.na(HR),0,HR)
    FFMC   <- ifelse(FFMC < 0 | FFMC > 101.,0,FFMC)
    FFMC   <- ifelse(is.na(FFMC),90,FFMC)
    ISI    <- ifelse(is.na(ISI)|ISI < 0 | ISI > 300,0,ISI)
    BUI    <- ifelse(BUI < 0 | BUI > 1000,0,BUI)
    BUI    <- ifelse(is.na(BUI),60,BUI)
    WS     <- ifelse(WS < 0 | WS > 300,0,WS)
    WS     <- ifelse(is.na(WS),10,WS)
    WD     <- ifelse(is.na(WD)|WD < -2*pi | WD > 2*pi,0,WD)
    GS     <- ifelse(is.na(GS)|GS < 0 | GS > 200,0,GS)
    GS     <- ifelse(ASPECT < -2*pi | ASPECT > 2*pi,0,GS)
    PC     <- ifelse(is.na(PC)|PC < 0 | PC > 100,50,PC)
    PDF    <- ifelse(is.na(PDF)|PDF < 0 | PDF > 100,35.0,PDF)
    CC     <- ifelse(CC <= 0 | CC > 100,95,CC)                                                                                         # originally "c"
    CC     <- ifelse(is.na(CC),80,CC)                                                                                         # originally "c"
    GFL    <- ifelse(is.na(GFL)|GFL <= 0 | GFL > 100,0.35,GFL)                                                                                    # changed from 0.3 to 0.35, pg 6 - 2009 */
    LAT    <- ifelse(LAT < -90 | LAT > 90,0,LAT)
    LAT    <- ifelse(is.na(LAT),55,LAT)
    LONG    <- ifelse(LONG < -180 | LONG > 360,0,LONG)
    LONG    <- ifelse(is.na(LONG),-120,LONG)
    THETA  <- ifelse(is.na(THETA)|THETA  < -2*pi | THETA > 2*pi,0,THETA)
    SD     <- ifelse(SD < 0 | SD > 100000, -999,SD)
    SD     <- ifelse(is.na(SD), 0,SD)
    SH     <- ifelse(SH < 0 | SH > 100, -999,SH)
    SH     <- ifelse(is.na(SH),0,SH)
  }

FUELTYPE<-sub("-","",FUELTYPE)
FUELTYPE<-sub(" ","",FUELTYPE)
# Convert time from hours to minutes */
HR     <- HR*60.
# Corrections to reorient WAZ, SAZ */
WAZ    <- WD + pi
WAZ    <- ifelse(WAZ > 2*pi,WAZ-2*pi,WAZ)
# nb: BMW's data set appears to have ASPECT not SAZ */

SAZ    <- ASPECT + pi
SAZ    <- ifelse(SAZ > 2*pi,SAZ-2*pi,SAZ)

# Make LONG positive for the Western Hemisphere */
LONG<-ifelse(LONG < 0, -LONG,LONG)

## /* Enter FBP Calculations */
## Initialize the output variables.
SFC    <- TFC<-HFI<-CFB<-ROS<-rep(0,length(LONG))                                                                                       # value 0 means non-fuel or ffmc == 0.
RAZ    <- rep(-999,length(LONG))                                                                                                        # value -999 means non-fuel or ffmc == 0. Nobody cares about the wind direction in such area.
if (output=="SECONDARY"|output=="ALL"|output =="S"|output =="A"){
  FROS <- BROS<-TROS<-HROSt<-FROSt<-BROSt<-TROSt<-FCFB<-
  BCFB <- TCFB<-FFI<-BFI<-TFI<-FTFC<-BTFC<-TTFC<-rep(0,length(LONG))
  TI   <- FTI<-BTI<-TTI<-LB<-WSV<- rep(-999,length(LONG))
  }

#/* presently, we do not accept a zero CBH; use near zero if necessary */
CBHs   <- c(2,3,8,4,18,7,10,0,6,6,6,6,0,0,0,0,0)
names(CBHs)<-c("C1","C2","C3","C4","C5","C6","C7","D1","M1","M2","M3","M4","S1","S2","S3","O1A","O1B")

CBH    <- ifelse(CBH <= 0 | CBH > 50 | is.na(CBH), ifelse(FUELTYPE %in% c("C6")&SD>0&SH>0,
        -11.2 + 1.06*SH + 0.00170*SD,                                                                                                  #/* 91 */
        CBHs[FUELTYPE]),CBH)
CBH    <- ifelse(CBH <0,0.0000001,CBH)

#/* presently, we do not accept a zero CFL,; use near zero if necessary */
CFLs   <- c(0.75,0.80,1.15,1.20,1.20,1.80,0.5,0,0.80,0.80,0.80,0.80,0,0,0,0,0)
names(CFLs)<-c("C1","C2","C3","C4","C5","C6","C7","D1","M1","M2","M3","M4","S1","S2","S3","O1A","O1B")
CFL    <- ifelse(CFL <= 0|CFL>2.0|is.na(CFL),CFLs[FUELTYPE],CFL)

FMC    <- ifelse(FMC <= 0 | FMC > 120 | is.na(FMC),.FMCcalc(LAT, LONG, ELV, DJ, D0),FMC)
FMC    <- ifelse(FUELTYPE %in% c("D1","S1","S2","S3","O1A","O1B"),0,FMC)
SFC    <- .SFCcalc(FUELTYPE, FFMC, BUI, PC, GFL)
BUI    <- ifelse(BUIEFF !=1,0,BUI)
WSV0   <- .Slopecalc(FUELTYPE, FFMC, BUI, WS, WAZ, GS,SAZ, FMC, SFC, PC, PDF, CC, CBH,ISI,output="WSV")                                     #/* This turns off BUI effect */
WSV    <- ifelse(GS > 0 & FFMC > 0,WSV0,WS)
RAZ0   <- .Slopecalc(FUELTYPE, FFMC, BUI, WS, WAZ, GS,SAZ, FMC, SFC, PC, PDF, CC, CBH,ISI,output="RAZ")
RAZ    <- ifelse(GS > 0 & FFMC > 0,RAZ0,WAZ)
  ISI    <- ifelse(ISI > 0,ISI,.ISIcalc(FFMC, WSV))
ROS    <- ifelse(FUELTYPE %in% c("C6"),.C6calc(FUELTYPE,ISI,BUI,FMC,SFC,CBH,option="ROS"),.ROScalc(FUELTYPE,ISI,BUI,FMC,SFC,PC,PDF,CC,CBH))
CFB    <- ifelse(FUELTYPE %in% c("C6"),.C6calc(FUELTYPE,ISI,BUI,FMC,SFC,CBH,option="CFB"),ifelse(CFL>0,.CFBcalc(FUELTYPE,FMC,SFC,ROS,CBH),0))

#CFB    <- ifelse(CFL==0,0,.CFBcalc(FUELTYPE, FMC, SFC, ROS, CBH))

TFC    <- .TFCcalc(FUELTYPE, CFL, CFB, SFC, PC, PDF)
HFI    <- .FIcalc(TFC,ROS)
CFB    <- ifelse(HR < 0,-CFB,CFB)
RAZ    <- RAZ * 180/pi
RAZ    <- ifelse(RAZ==360,0,RAZ)                                                                                                       # 360 degree is the same as 0, FBP use 0 instead (Wotton etal 2009)
FD     <- rep("I",length(CFB))
FD     <- ifelse(CFB<0.10,"S",FD)
FD     <- ifelse(CFB>=0.90,"C",FD)
CFC    <- .TFCcalc(FUELTYPE, CFL, CFB, SFC, PC, PDF,option="CFC")
if (output=="SECONDARY"|output=="ALL"|output =="S"|output =="A"){
    SF     <- ifelse (GS >= 70,10,exp(3.533 * (GS/100)^1.2))                                                                         # /* 39 */
    CSI    <- .CFBcalc(FUELTYPE,FMC,SFC,ROS,CBH,option="CSI")
    RSO    <- .CFBcalc(FUELTYPE,FMC,SFC,ROS,CBH,option="RSO")
    BE     <- .BEcalc(FUELTYPE,BUI)
    LB     <- .LBcalc(FUELTYPE, WSV)
    LBt    <- ifelse(ACCEL == 0,LB,.LBtcalc(FUELTYPE,LB,HR,CFB))
    BROS   <- .BROScalc(FUELTYPE,FFMC,BUI,WSV,FMC,SFC,PC,PDF,CC,CBH)
    FROS   <- .FROScalc(ROS, BROS, LB)
   #/* TROS is the rate of spread towards angle THETA */
    E      <- sqrt(1-1/LB/LB)                                                                                                          #/* eccentricity */
    TROS   <- ROS * (1-E)/(1-E*cos(THETA - RAZ))                                                                                       #/* note: this is the old method using the focus as the ignition point */
  #//   TROS <- ROSthetacalc(ROS, FROS, BROS, THETA)                                                                                   #MARC: what is this?
    ROSt   <- ifelse(ACCEL==0,ROS,.ROStcalc(FUELTYPE, ROS, HR, CFB))
    BROSt  <- ifelse(ACCEL==0,BROS,.ROStcalc(FUELTYPE, BROS, HR, CFB))
    FROSt  <- ifelse(ACCEL==0,FROS,.FROScalc(ROSt, BROSt, LBt))
    TROSt  <- ifelse(ACCEL==0,TROS,ROSt * (1.-sqrt(1.-1./LBt/LBt))/(1.-sqrt(1.-1./LBt/LBt)*cos(THETA - RAZ)))                          #/* note: this is the old method using the focus as the ignition point */
    FCFB   <- ifelse(CFL==0,0,ifelse(FUELTYPE %in% c("C6"),0,.CFBcalc(FUELTYPE, FMC, SFC, FROS, CBH)))
    BCFB   <- ifelse(CFL==0,0,ifelse(FUELTYPE %in% c("C6"),0,.CFBcalc(FUELTYPE, FMC, SFC, BROS, CBH)))
    TCFB   <- ifelse(CFL==0,0,ifelse(FUELTYPE %in% c("C6"),0,.CFBcalc(FUELTYPE, FMC, SFC, TROS, CBH)))
    FTFC   <- .TFCcalc(FUELTYPE, CFL, FCFB, SFC, PC, PDF)
    BTFC   <- .TFCcalc(FUELTYPE, CFL, BCFB, SFC, PC, PDF)
    TTFC   <- .TFCcalc(FUELTYPE, CFL, TCFB, SFC, PC, PDF)
  #/* equilibrium values */
    FFI    <- .FIcalc(FTFC, FROS)
    BFI    <- .FIcalc(BTFC, BROS)
    TFI    <- .FIcalc(TTFC, TROS)

  #/* For now... */
    HROSt  <- ifelse(HR < 0,-ROSt,ROSt)
    FROSt  <- ifelse(HR < 0,-FROSt,FROSt)
    BROSt  <- ifelse(HR < 0,-BROSt,BROSt)
    TROSt  <- ifelse(HR < 0,-TROSt,TROSt)

      a1<-0.115-(18.8*CFB^2.5*exp(-8*CFB))
      TI<-log(ifelse(1-RSO/ROS > 0,1-RSO/ROS,1))/(-a1)
      
      a2<-0.115-(18.8*FCFB^2.5*exp(-8*FCFB))
      FTI<-log(ifelse(1-RSO/FROS > 0,1-RSO/FROS,1))/(-a2)
      
      a3<-0.115-(18.8*BCFB^2.5*exp(-8*BCFB))
      BTI<-log(ifelse(1-RSO/BROS > 0,1-RSO/BROS,1))/(-a3)
      
      a4<-0.115-(18.8*TCFB^2.5*exp(-8*TCFB))
      TTI<-log(ifelse(1-RSO/TROS > 0,1-RSO/TROS,1))/(-a4)
      
      DH <-ifelse(ACCEL==1,.DISTtcalc(FUELTYPE, ROS, HR, CFB),ROS*HR  )
      DB <- ifelse(ACCEL==1,.DISTtcalc(FUELTYPE, BROS, HR, CFB),BROS*HR  )
      DF <- ifelse(ACCEL==1, (DH+DB)/(LBt*2), (DH+DB)/(LB*2)   )
  }
  if (exists("ID")) ID<-ID else ID <- row.names(input)  
  if (output == "PRIMARY"|output == "P"){
    FBP    <- data.frame(ID,CFB,CFC,FD,HFI,RAZ,ROS,SFC,TFC)
    FBP[,c(2:3,5:ncol(FBP))]<-apply(FBP[,c(2:3,5:ncol(FBP))],2,function(.x) ifelse(FUELTYPE %in% c("WA","NF"),0,.x))
    FBP[,"FD"]<-as.character(FBP[,"FD"])
    FBP[,"FD"]<-ifelse(FUELTYPE %in% c("WA","NF"),"NA",FBP[,"FD"])
    FBP} else
  if (output == "SECONDARY"|output == "S"){
    FBP    <- data.frame(ID,BE,SF,ISI,FFMC,FMC,D0,RSO,CSI,FROS,BROS,HROSt,FROSt,BROSt,FCFB,BCFB,FFI,BFI,FTFC,BTFC,TI,FTI,BTI,LB,LBt,WSV,DH,DB,DF,TROS,TROSt,TCFB,TFI,TTFC,TTI)
    FBP[,2:ncol(FBP)]<-apply(FBP[,2:ncol(FBP)],2,function(.x) ifelse(FUELTYPE %in% c("WA","NF"),0,.x))
    FBP} else
  if (output == "ALL"|output == "A") {
    FBP    <- data.frame(ID,CFB,CFC,FD,HFI,RAZ,ROS,SFC,TFC,BE,SF,ISI,FFMC,FMC,D0,RSO,CSI,FROS,BROS,HROSt,FROSt,BROSt,FCFB,BCFB,FFI,BFI,FTFC,BTFC,TI,FTI,BTI,LB,LBt,WSV,DH,DB,DF,TROS,TROSt,TCFB,TFI,TTFC,TTI)
    FBP[,c(2:3,5:ncol(FBP))]<-apply(FBP[,c(2:3,5:ncol(FBP))],2,function(.x) ifelse(FUELTYPE %in% c("WA","NF"),0,.x))
    FBP[,"FD"]<-as.character(FBP[,"FD"])
    FBP[,"FD"]<-ifelse(FUELTYPE %in% c("WA","NF"),"NA",FBP[,"FD"])
    FBP} 

}

