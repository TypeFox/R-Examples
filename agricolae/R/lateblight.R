lateblight <-
function(WS, Cultivar,ApplSys,InocDate, LGR,IniSpor,SR,IE, LP, InMicCol,
MatTime=c('EARLYSEASON','MIDSEASON','LATESEASON'),...)
{
MatTime <- toupper(MatTime)
MatTime <- match.arg(MatTime)
Cultivar <- toupper(Cultivar)
ApplSys <- toupper(ApplSys)
MatTime <- toupper(MatTime)
Wfile<-WS$Wfile
Sfile<-WS$Sfile
EmergDate<-WS$EmergDate
EndEpidDate<-WS$EndEpidDate
# Start model lateblight
SeLenght <- nrow(Wfile)
# Obtaining nInocDay
nInocDayFile <-subset(Wfile, Wfile[,1] <= InocDate)
nInocDay     <- nrow( nInocDayFile )
# End obtaining nInocDay
# Obtaining date after emergence
nday<-as.numeric(Wfile[,1]-EmergDate)
ModelFile<-data.frame(EmergDate,Date=Wfile[,1],nday,InocDate,EndEpidDate)
# End obtaining date after emergence
# Initialization and declaration of variables
  ActAr <- 0
  Age <- 0
  Area <- 0
  Catch <- 0
  ColFct <- 0
  ConstLfWet <- 3
  DenFct <- 0
  DirGrRt <- 0
  Divisor <- 0
  ExpCol <- 0
  fTRH <- 0
  FunEff <- 1
  GrmZoo <- 0
  GrSprg <- 0
  HoRHSpor <- 6
  HotFct <- 0
  HumFct <- 0
  IncrFct <- 0
  IndGrRt <- 0
  IE <- IE
  InfSpRt <- 0
  InfZoRt <- 0
  IniLesDia <- 0.002
  lagBlPr <- 0
  LeafFct <- 0
  LfWetFct <- 0
  LGR <- LGR
  LGRRt <- 0
  MaxDirGr <- 0.29
  MaxIndGr <- 0.81
  MaxInfSp <- 0.10
  MaxInfZo <- 0.01
  MaxZooGr <- 0.96
  n5DayPer <- 0
  nday <- 0
  nDay2 <- 0
  nDay3 <- 0
  rLGR <- 0
  rSR <- 0
  SporFct <- 0
  SSAfrac <- 0
  Tdirect <- 0
  Tindir <- 0
  TmpFct <- 0
  TotSpor <- 0
  TZooGrm <- 0
  Zoo <- 0
  ZooGrRt <- 0
  ZooSpr <- 6.567
  LesAr<-rep(0,16)
  Early <- c(0.006, 0.06, 0.19, 0.45, 0.81,1.30, 1.88, 2.20, 2.42, 2.50)
  Mid <- c(0.006,0.06,0.19,0.45,0.81,1.30,1.88,2.56,2.99,3.28,3.47,3.53)
  Late <- c(0.006,0.06,0.19,0.45,0.81,1.30,1.88,2.76,3.29,3.74,4.16,4.38,4.50)
  NumCol <- rep(0,16)
  LesDia <- rep(0,16)
  ArLP <- rep(0,SeLenght)
  SimSeverity <- rep(0,SeLenght)
  BlPr <- rep(0,SeLenght)
  AUDPC <- rep(0,SeLenght)
  MicCol <- rep(0,SeLenght)
  AttchSp <- rep(0,SeLenght)
  LatPer <- rep(0,SeLenght)
  LesExInc <- rep(0,SeLenght)
  LAI <- rep(0.01,SeLenght)
  LPAux <- 0
  rLP <- rep(0,SeLenght)
  InvrLP <- rep(0,SeLenght)
  Defol <- rep(0,SeLenght)
  InArea <-  rep(0,SeLenght)
  TotCol <-  rep(0,SeLenght)
  Sprg <- rep(0,SeLenght)
  SSA <-  rep(0,SeLenght)
# End initialization and declaration of variables
# Start simulation
Rainfall <- Wfile$Rainfall
Tmp      <- Wfile$Tmp
HumidTmp <- Wfile$humidtmp
HumidHrs <- Wfile$HumidHrs
for( nday in 2:(SeLenght)) {
# set work.Wfile (Obs = &SeLenght);
if( !is.na( Rainfall[nday]) ){
Rainfall[nday] <- Rainfall[nday] / 10  
# To go from rainfall in mm to rainfall in cm(for the fungicide sub-model)
}
# Temperature-dependent latent period;
LPAux <- LP
if( !is.na( Tmp[nday]) ){
  rLP[nday] <- 4.571693 - 0.310652 * Tmp[nday]+ 0.006755 * Tmp[nday]^2
  InvrLP[nday] <- 1 / rLP[nday]
}
  for (nDay2 in 1:(nday - 1)) {
  ArLP[nDay2] <- ArLP[nDay2] + InvrLP[nday];
  }
nDay3 <- 0
for (nDay2 in 1:(nday - 1)) {
  if (ArLP[nDay2] > (LPAux - 1)) nDay3 <- nDay2
  }
  if (nDay3 > 0) {
      SSAfrac <- LPAux - ArLP[nDay3]
      LatPer[nday] <- nday - nDay3
      }
  if (nDay3 <= 0) {
      SSAfrac <- 1
      if(InvrLP[nday]!=0)LatPer[nday] <- round(LPAux/InvrLP[nday],0)
    }
if( !is.na( LatPer[nday]) ){
    if (LatPer[nday] < 2) LatPer[nday] <- 2
    if (LatPer[nday] > 8) LatPer[nday] <- 8
    if (InvrLP[nday] < 0) {
  SSAfrac <- 1
    LatPer[nday] <- round(LPAux,0)
  }
}
# End of the temperature-dependent latent period
# Environment sub-model
    LfWetHrs <- HumidHrs[nday] + ConstLfWet
    LfWetHrs <- max (0.0, LfWetHrs)
    LfWetHrs <- min (24.0, LfWetHrs)
    LfWetTmp <- (Tmp[nday] + HumidTmp[nday]) / 2.0
# End environment sub-model
# Pathogen sub-model
# Infection
    Sprg[nday] <- Sprg[nday-1] + IniSpor

    if (nday > nInocDay) IniSpor <- 0
# Direct germination of sporangia
if( !is.na( LfWetTmp ) ){
    if ((LfWetTmp>=14.0) & (LfWetTmp <= 29.0)) {
    TDirect <- 2.1 +LfWetTmp *(-.48 + LfWetTmp*(.031 + 
    LfWetTmp*(-4.5E-04 -4.4E-06*LfWetTmp)))
    DirGrRt <- Tdirect * MaxDirGr
    GrSprg <- Sprg[nday] * DirGrRt
    }
    else  GrSprg <- 0.0
# Formation of zoospores
    if ((LfWetTmp >= 0.0) & (LfWetTmp <= 24.0)) {
      TIndir <- 1.314E-02 + LfWetTmp * (6.98E-02 + LfWetTmp*(1.815E-02 + 
      LfWetTmp * (-2.23E-03 + LfWetTmp * (7.48E-05 - 7.68E-07 * LfWetTmp))))
      IndGrRt <- TIndir * MaxIndGr * ZooSpr
      Zoo <- Sprg[nday] * IndGrRt
    }
    else Zoo <- 0.0
# Germination of zoospores
    if (( LfWetTmp >= 0.0) & ( LfWetTmp <= 28)) {
      TZooGrm <- 4.94E-03 + LfWetTmp*(.364 + LfWetTmp*(-5.813E-02 + LfWetTmp*
                (4.52E-03 + LfWetTmp*(-1.617E-04 + 2.068E-06 * LfWetTmp))))
      ZooGrRt <- TZooGrm * MaxZooGr
      GrmZoo <- Zoo * ZooGrRt
    }
    else GrmZoo <- 0.0
}
# Calculation of LfWetFct
if( !is.na( HumidTmp[nday] ) ){
    TmpFct <- min (1.0, 0.03 * HumidTmp[nday] + 0.4)
    if ((LfWetHrs * TmpFct) > 8.0 ) {
      LfWetFct <- min (0.8 + 0.0125 * (LfWetHrs * TmpFct - 8.0), 1.0)
    }
    else   LfWetFct <- max (0.0, 0.16 * (LfWetHrs * TmpFct - 3.0))

# Infection;
    if (( HumidTmp[nday] >= 1.5) & (HumidTmp[nday] <= 30.0) ) {
      fTRH <- -.4846 + HumidTmp[nday] * (.3735 + HumidTmp[nday]*(-3.556E-02 + 
      HumidTmp[nday] * (1.509E-03 - 2.405E-05 * HumidTmp[nday])))
      InfSpRt <- MaxInfSp * fTRH * LfWetFct * FunEff * IE
      InfZoRt <- MaxInfZo * fTRH * LfWetFct * FunEff * IE
      MicCol[nday] <- GrSprg * InfSpRt + GrmZoo * InfZoRt
    }
    else MicCol[nday] <- 0.0
}
# End of infection
    if (nday == nInocDay) {
        MicCol[nday] <- InMicCol
        NumCol[1] <- MicCol[nday]
    }
# Lesion expansion
if( !is.na( Tmp[nday] ) ){
    if ( (Tmp[nday] >= 8) & (Tmp[nday] <= 32.94)) {
    rLGR<-0.126636-0.038378*Tmp[nday]+0.008593*Tmp[nday]^2-0.000229*Tmp[nday]^3
    }
    if ((Tmp[nday] >= 4.09) & (Tmp[nday] < 8))rLGR<-0.064246*Tmp[nday]-0.262690
    if (Tmp[nday] > 32.94) rLGR <- 0.0
    if (Tmp[nday] < 4.09) rLGR <- 0.0

    LesExInc[nday] <- rLGR * LGR
}
# Effect of disease density
    if (nday > nInocDay ) {
     for( Age in seq(16,2,-1)) NumCol[Age] <- NumCol[Age - 1]
    }
    if (nday > nInocDay) {
      NumCol[1] <- MicCol[nday]
    }
# New latent periods
if( !is.na( LatPer[nday] ) ){
         if (LatPer[nday] == 2) ExpCol <- sum( NumCol[2:15] )
         if (LatPer[nday] == 3) ExpCol <- sum( NumCol[3:15] )
         if (LatPer[nday] == 4) ExpCol <- sum( NumCol[4:15] )
         if (LatPer[nday] == 5) ExpCol <- sum( NumCol[5:15] )
         if (LatPer[nday] == 6) ExpCol <- sum( NumCol[6:15] )
         if (LatPer[nday] == 7) ExpCol <- sum( NumCol[7:15] )
         if (LatPer[nday] == 8) ExpCol <- sum( NumCol[8:15] )
if(ExpCol>0.0)ColFct <-min(1.0, 2.7 - 0.32*log10(ExpCol)) else ColFct <- 0.0
}
# Calibrating the effect of BlPr on lesion growth
LeafFct <- ifelse(BlPr[nday] < 0.82, exp (-2.1 * BlPr[nday-1]), 
                  exp (-2.1 * 0.9 * BlPr[nday-1]))
DenFct  <- max (0.0, ColFct * LeafFct)
# Update inactive lesion area by adding in class 15
# Calculate lesion diameter
if( !is.na( LatPer[nday] ) ){
    for( Age  in seq(16,LatPer[nday],-1)) {
      IncrFct <- max (0.0, min (1.0, 1.72 - (0.12 * (Age - 1))))
      LGRRt <- LesExInc[nday] * IncrFct * DenFct
      LesDia[Age] <- LesDia[Age-1] + LGRRt
    }
    LesDia[LatPer[nday]] <- IniLesDia * DenFct
# Calculate lesion area;
    for( Age in LatPer[nday]:16) {
      LesAr[Age] <- NumCol[Age] * ((LesDia[Age]^2) * 0.785398163397448)
    }
    ActAr <- sum(LesAr[1:15])
    InArea[nday] <-InArea[nday-1] + LesAr[16]
    Defol[nday] <- ActAr + InArea[nday]
}
# End lesion expansion
# Sporulation
if( !is.na( HumidHrs[nday] ) & !is.na(HumidTmp[nday] )){
    HumFct <- min(max (((HumidHrs[nday] - HoRHSpor) / (24 - HoRHSpor)), 0), 1)
    SSA[nday] <- SSAfrac * LesAr[LatPer[nday]]
    for (Age in (LatPer[nday] + 1):15) {
      SSA[nday]<-SSA[nday] + LesAr[Age]
    }
    if ((HumidTmp[nday] >= 7.26) & (HumidTmp[nday] <= 25.29))rSR <- -4.796667 + 
    1.021578*HumidTmp[nday]-0.056437*HumidTmp[nday]^2+0.000931*HumidTmp[nday]^3
    if (HumidTmp[nday] > 25.29) rSR <- 0.0
    if (HumidTmp[nday] < 7.26) rSR <- 0.0
    if (rSR > 0) {
      SporFct <- max(0, (.36 * SR * rSR * HumFct))
      AttchSp[nday] <- SporFct * SSA[nday]
    }
    if (rSR <= 0) {
      SporFct <- 0
      AttchSp[nday] <- 0
    }

    TotSpor <- TotSpor + AttchSp[nday]
}
# End sporulation
# Dispersal
if( !is.na( Tmp[nday] ) ){
    if (Tmp[nday] > 21) HotFct <- exp (6 * (- 0.0256 * Tmp[nday] + .537))
    else HotFct <- 1
    Catch <- .025 * (1 - BlPr[nday]) * HotFct
    Sprg[nday] <- Catch * AttchSp[nday]
}
# End dispersal
# End pathogen sub-model
# Disease sub-model
    BlPr[nday] <- ifelse(BlPr[nday]>0.975, 1,max((Defol[nday]/LAI[nday-1]),0))
    if (BlPr[nday] > 1 ) BlPr[nday] <- 1
    SimSeverity[nday] <- BlPr[nday] * 100
    lagBlPr <- BlPr[nday-1]
    Area <- ((BlPr[nday] + lagBlPr) / 2)
    AUDPC[nday] <- AUDPC[nday-1] + Area
    TotCol[nday] <- TotCol[nday-1] + MicCol[nday]
# End disease sub-model
# Host growth sub-model
    n5DayPer <- nday %/% 5 + 1
    Divisor <- n5DayPer - 1
    if (MatTime=="EARLYSEASON")LAI[nday]<-Early[n5DayPer]+(((Early[n5DayPer + 1]
        - Early[n5DayPer]) / 5) * (nday-5*Divisor))
    if ( (MatTime == "EARLYSEASON")&(nday > 44))LAI[nday] <- 2.5
    if ( MatTime=="MIDSEASON" )LAI[nday]<- Mid[n5DayPer] + (((Mid[n5DayPer + 1]
        - Mid[n5DayPer]) / 5) * (nday - 5 * Divisor))
    if ( (MatTime == "MIDSEASON") &(nday > 52)) LAI[nday] <- 3.5
    if (MatTime=="LATESEASON" )LAI[nday]<-Late[n5DayPer]+ (((Late[n5DayPer + 1]
       - Late[n5DayPer]) / 5) * (nday - 5 * Divisor))
    if ( (MatTime == "LATESEASON") & (nday > 59) ) LAI[nday] <- 4.5
# End host growth sub-model
}
# End simulation
Ofile <- data.frame(Date=ModelFile$Date,nday=ModelFile$nday,
MicCol,SimSeverity,LAI, LatPer,LesExInc,AttchSp , AUDPC = AUDPC * 100,
rLP,InvrLP,BlPr,Defol)
# graphing disease progress curves
Gfile <-subset(Ofile,((Ofile$ApplSys == ApplSys )& ( Ofile$Cultivar == Cultivar)))
MinObs <- Sfile$MeanSeverity - Sfile$StDevSeverity
MaxObs <- Sfile$MeanSeverity + Sfile$StDevSeverity
Gfile <- data.frame(Sfile[,c(-1,-2)],MinObs,MaxObs)
plot(Ofile$nday,Ofile$SimSeverity,...)
points(Gfile$nday,Gfile$MeanSeverity,col="blue",cex=1,pch=19)
npoints <- nrow( Gfile )
for ( i in 1:npoints){
segments(Gfile[i,2],Gfile[i,5],Gfile[i,2],Gfile[i,6],lwd=1.5,col="blue")
}
legend("topleft",c("Disease progress curves","Weather-Severity"),
title=MatTime,lty=1,pch=c(3,19),col=c("black","blue"),cex=1)
# End graphing disease progress curves
return( list(Ofile = Ofile, Gfile=Gfile ))
}
