Chi2testMixtures <- function(Data,Means,SDs,Weights,IsLogDistribution = Means*0,PlotIt = 0,UpperLimit=max(Data,na.rm=TRUE),VarName='Data'){
# V=Chi2testMixtures(Data,Means,SDs,Weights,IsLogDistribution,PlotIt,UpperLimit)
# V$Pvalue V$BinCenters V$ObsNrInBin V$ExpectedNrInBin
#  chi-square test of Data vs a given Gauss Mixture Model
#  komfortabler Plot mit allen Mess und Fehlerwerten sowie der cdf fuer die entsprechende Chi-Quadrat funktion
#
# INPUT
# Data(1:N)                        data points
# Means(1:L,1)                     Means of Gaussians,  L ==  Number of Gaussians
# SDs(1:L,1)                     estimated Gaussian Kernels = standard deviations
# Weights(1:L,1)                   relative number of points in Gaussians (prior probabilities): 
#
# OPTIONAL
# IsLogDistribution(1:L,1) oder 0  if IsLogDistribution(i)==1, then mixture is lognormal
#                                  default == 0*(1:L)'
# PlotIt                           do a Plot of the compared cdf's and the KS-test distribution (Diff)
# UpperLimit                       test only for Data <= UpperLimit, Default = max(Data) i.e all Data.
# VarName                          Variable Name for Plots
# OUTP
# Pvalue                           Pvalue of a suiting chi-square , Pvalue ==0 if Pvalue <0.001
# BinCenters                       bin centers
# ObsNrInBin                       nr of data in bin
# ExpectedNrInBin                  nr of data that should be in bin according to GMM
# Chi2Value                        the  TestStatistic i.e.:
#                                  sum((ObsNrInBin(Ind)-ExpectedNrInBin(Ind)).^2./ExpectedNrInBin(Ind)) with
#                                  Ind = find(ExpectedNrInBin>=10);

# uses  histopt_hlp,CDFMixtures,randomLogMix_hlp,PlotMixtures       matlab's histc(..)
# Autor: RG 06/15
#1.Editor: MT 08/2015 plotten, Fehlerabfang bei kleinen Datensaetzen
  
par.defaults <- par(no.readonly=TRUE)
if(length(IsLogDistribution) == 1 && IsLogDistribution == 0) 
  IsLogDistribution = Means*0
NumberOfData=length(Data) #s.Zeile 128, Zeile 110
DataSmall=80 #Abwann ist ein Datensatz klein
#####################################################
# mit histopt-Algorithmus die Bins bestimmen, anzahl= OptimalNrOfBins
histopt_hlp <- function(Data){
# histopt(Data);
# Histogram optimaler Binanzahl. Bins besitzen alle die gleiche Groess?e.
#
# INPUT
# Data[d,1]	   		Vektor der zu zeichneten Variable
# OUTPUT
# nrOfBins			Anzahl der Bins
# nrInBins[1,d] 	Startpunkt jedes Bins als Vektor
# binMids			Mitte des Bins ??
# Letzte Ergaenzung MT, Autor unbekannt
# Ergaenzung RG, Autor unbekannt

  Data[is.infinite(Data)] = NA #MT: Korrektur, bereinigung von Inf
  optNrOfBins<-OptimalNoBins(Data)
  optNrOfBins = min(100,optNrOfBins) #RG: Aus Matlab uebernommen
# temp<-hist(Data,breaks=optNrOfBins)
 #print(optNrOfBins[1])
 #temp <- hist(Data, breaks=optNrOfBins[1], col="blue", border="light blue", main=Title)
 
 minData <- min(Data,na.rm = TRUE)
 maxData <- max(Data,na.rm = TRUE)
 i <- maxData-minData
 optBreaks <- seq(minData, maxData, i/optNrOfBins) # bins haben alle die gleiche groesse
 temp <- hist(Data, breaks=optBreaks, plot=FALSE)
 
 #box();
 Breaks <- temp$breaks
 nB <- length(Breaks)
 y <- temp$counts;
invisible(list(nrOfBins=length(Breaks)-1, nrInBins=y, binMids=temp$mids))
 }
#####################################################
ho <-histopt_hlp(Data)
OptimalNrOfBins <- ho$nrOfBins
ObsNrInBin <- ho$nrInBins
BinCenters <- ho$binMids

# BinLimits errechnen und mit histc die ObsNrInBin nachrechnen;
binwidth = diff(BinCenters)                        # differences between adjacent BinCenters
binwidth = c(binwidth,binwidth[length(binwidth)])
xx = cbind(BinCenters-binwidth/2,BinCenters+binwidth/2)
xx[1] = min(xx[1],which.min(Data))
xx[length(xx)] = max(xx[length(xx)],which.max(Data))
# in xx[,1] stehen die unteren Grenzen der Bins, in xx(:,2) die oberen Grenzen
# Shift bins so the interval is ( ] instead of [ ).
xx = Re(xx)
#xx =xx +eps(xx); #####eps!
BinLimits = c(xx[1,1],xx[,2])                         # dies sind jetzt die Bin Grenzen

# # Nachrechnen: fuer diese Bin Grenzen mit histc schauen wieviele in die jeweilgen bins fallen
# nn = histc(full(real(Data)),BinLimits);              % matlab' schnelle counting Funktion benutzen
# # Combine last bin with next to last bin
# nn(end-1) = nn(end-1)+nn(end);
# nn = nn(1:end-1);
# [nn,ObsNrInBin]
# AllOK = sum(nn-ObsNrInBin)==0  % wenn die mit histopt gerechnete anz der mit histc ger. uebeeinstimmt

# jetzt die erwartete Anzahl berechnen
CDFGaussMixture = CDFMixtures(BinLimits,Means,SDs,Weights,IsLogDistribution)$CDFGaussMixture # cdf(GMM)

AnzData =length(Data)
ExpectedNrInBinCDF = CDFGaussMixture*AnzData
ExpectedNrInBin = diff(ExpectedNrInBinCDF)


# jetzt den Vergleich Anstellen
AnzDiff = ObsNrInBin-ExpectedNrInBin
Chi2Diffs = AnzDiff*0 # init

if(NumberOfData<DataSmall){
    Ind = which(ExpectedNrInBin>=2, arr.ind=TRUE) #Bei Kleinen Datensetzen schranke runterstellen
    warning(paste('Chi2testMixtures(): Datasize',NumberOfData,'is too small. Pvalue could be misleading'))
}else{
    Ind = which(ExpectedNrInBin>=10, arr.ind=TRUE) #  ObsNrInBin mindestens 10 sonst gelten die Werte als identisch => Diff ==0
}
Chi2Diffs[Ind] =((ObsNrInBin[Ind]-ExpectedNrInBin[Ind])^2)/ExpectedNrInBin[Ind]
Chi2Value = sum(Chi2Diffs) # dies soll angeblich Chi2 Verteilt sein


# Die Chi2 Funktion via Monte-Carlo errechnen
#AnzData =length(Data);
AnzBins = length(BinCenters);
AnzRepetitions = 1000;
if(AnzBins<100)  AnzRepetitions = 2000
if(AnzBins<10)   AnzRepetitions = 5000


RandGMMDataDiff = matrix(0,AnzBins,AnzRepetitions)
for(i in 1:AnzRepetitions){
  R = RandomLogGMM(Means,SDs,Weights,IsLogDistribution,AnzData)
  #BinLimits = c(0,BinLimits,max(abs(R)))
  RandNrInBin = hist(Re(abs(R)),c(0,BinLimits,max(abs(R))),plot=F)$counts  # R's schnelle Funktion benutzen
  #BinLimits = BinLimits[2:(length(BinLimits)-1)]
  RandNrInBin[2] = RandNrInBin[1]+RandNrInBin[2]
  RandNrInBin[length(RandNrInBin)-1] = RandNrInBin[length(RandNrInBin)-1]+RandNrInBin[length(RandNrInBin)]
  RandNrInBin = RandNrInBin[3:length(RandNrInBin)-1]
  AnzDiffRand =  RandNrInBin-ExpectedNrInBin
  RandChi2Diffs = AnzDiffRand*0; # init
  if(NumberOfData<DataSmall){
    Ind = which(ExpectedNrInBin>=2, arr.ind=TRUE) #Bei Kleinen Datensetzen schranke runterstellen
  }else{
    Ind = which(ExpectedNrInBin>=10, arr.ind=TRUE) #  ObsNrInBin mindestens 10 sonst gelten die Werte als identisch => Diff ==0
  }
  RandChi2Diffs[Ind] = ( (RandNrInBin[Ind]-ExpectedNrInBin[Ind])^2)/ExpectedNrInBin[Ind]
  RandGMMDataDiff[,i] = RandChi2Diffs;
}

AllDiff =  colSums(RandGMMDataDiff);                # die Verteilung aller Differenzen
#[AllDiffCDF,AllDiffKernels] =  ecdfUnique(AllDiff); 
dummy <- ecdf(AllDiff) # cdf(Diff)
AllDiffCDF <- c(0,get("y", envir = environment(dummy)))
AllDiffKernels <- c(knots(dummy)[1],knots(dummy))#CDFuniq

# upper Limit Berucksichtigen
ClipInd = which(BinCenters<UpperLimit,arr.ind=TRUE)
BinCenters = BinCenters[ClipInd]
Chi2Value = sum(Chi2Diffs[ClipInd])
ExpectedNrInBin = ExpectedNrInBin[ClipInd]
ObsNrInBin = ObsNrInBin[ClipInd]
Chi2Diffs = Chi2Diffs[ClipInd]

# CHi2statistik berechnen
if (Chi2Value-AllDiffKernels[length(AllDiffKernels)] >1){ # Summe der Abweichungen liegt zu weit rechts
  Ch2cdfValue = 1;
} else{ #durch interpolation den wert bestimmen
#matlab: 
#Ch2cdfValue = interp1([0;AllDiffKernels],[0;AllDiffCDF],Chi2Value, 'linear');  #den MaxDiff in cdf(Diff) lokalisieren
  Ch2cdfValue = approx(rbind(0,unname(AllDiffKernels)),rbind(0,unname(AllDiffCDF)),Chi2Value, 'linear')$y;  # den MaxDiff in cdf(Diff) lokalisieren
} # if (Chi2Value-AllDiffKernels(end)) >1 # der wert liegt zu weit rechts

Pvalue = Ch2cdfValue                       # P- value fuer Chi sqare -test ausrechnen
Pvalue = round(Pvalue,5)                   # runden auf  gueltige stellen

if(PlotIt ==1){

  par(mfrow = c(2,2))
  #subplot(2,2,1)
  xlim=c(min(Data,na.rm = TRUE),max(Data,na.rm = TRUE))
  ylim=c(min(min(ObsNrInBin,na.rm = TRUE),min(ExpectedNrInBin,na.rm = TRUE)),max(max(ObsNrInBin,na.rm = TRUE),max(ExpectedNrInBin,na.rm = TRUE)))
  plot(BinCenters,ObsNrInBin,type='h',xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
  axis(1,xlim=xlim,col="black",las=1) #x-Achse
  axis(2,ylim=ylim,col="black",las=1) #y-Achse
  par(new = TRUE)
  for (i in 1:length(BinLimits)){ 
    points(c(BinLimits[i],BinLimits[i]),ylim,type='l',lwd=1,col='magenta');
    par(new = TRUE);
  }
  #hold on ;  
  par(new=TRUE)
  plot(BinCenters,ExpectedNrInBin,col='red',xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE)
  title(ylab=paste0('No. ',VarName,' in bin'))
  title(xlab=paste0(VarName,' and bins'))
  title('bar = observed, red= expected(GMM)')
  #axis tight;
  #ax=axis;
  
 
  #subplot(2,2,2)   # hier die fraglichen PDFs zeichnen    
  #xlim=c(min(abs(Data),na.rm=TRUE),max(abs(Data),na.rm=TRUE))
  plot(BinCenters,Chi2Diffs,col='red',xlab="",ylab="",type='b',xlim=xlim)
  grid() 
  #ax = axis;
  #yaxis(0,max(2,ax(4)));
  for (i in 1:length(BinLimits)){ 
    points(c(BinLimits[i],BinLimits[i]),ylim,type='l',lwd=1,col='magenta')
    par(new = TRUE);
  }
  #xaxis(min(BinLimits),max(BinLimits));
  title(xlab=paste0(VarName,', differences in bins'))
  title('squared bin differences : C2V=(Exp-Obs)^2/Exp')
  title(ylab='C2V');
 
  #subplot(2,2,3)
  Ylimits=c(0,1)
  plot(AllDiffKernels,AllDiffCDF,ylim=Ylimits,xlab="",ylab="")
  abline(v=Chi2Value,col='green')
  grid()
  abline(h=Ch2cdfValue,col='green');
  title(ylab=c('cdf(chi2)'))
  title(xlab=c('bl =Chi2;gn = sum(C2V) = ', paste(Chi2Value)))
  if (Pvalue==0){
    title(c('cdf(Chi2), Pvalue< 10e-4' ));
  }else{
    title(c('cdf(Chi2), Pvalue=', paste(Pvalue) ))
  }
  
  #subplot(2,2,4)
  Xlimits = c(min(Data,na.rm=TRUE),max(Data,na.rm=TRUE))
  #PDEplot(Data,xlim=Xlimits,ylim=Ylimits,defaultAxes=FALSE)
  paretoRadius<-ParetoRadius(Data)
  pdeVal        <- ParetoDensityEstimation(Data,paretoRadius,NULL)
  paretoDensity <- pdeVal$paretoDensity
  Ylimits = c(min(paretoDensity,na.rm=TRUE),max(paretoDensity,na.rm=TRUE))
  plot(pdeVal$kernels,paretoDensity,typ='l',col="blue",xlim = Xlimits, ylim = Ylimits, xlab = VarName, ylab = '',axes=FALSE,xaxs='i',yaxs='i') 
  title(ylab='PDE')

  #hold on; 
  par(new=TRUE)
  PlotMixtures(Data,Means,SDs,Weights=Weights,IsLogDistribution=IsLogDistribution,xlim=Xlimits,ylim=Ylimits,axes=FALSE,xlab="",ylab="",SingleGausses=T,xaxs='i',yaxs='i',MixtureColor='black', SingleColor = 'green')
  for (i in 1:length(BinLimits)){ 
    points(c(BinLimits[i],BinLimits[i]),ylim,type='l',lwd=1,col='magenta');
    par(new = TRUE);
  }
  #xaxis(min(BinLimits),max(BinLimits));
  title(paste0('black=pdf(GMM),green=pdf(',VarName,')'))
  axis(1,xlim=c(0,ceiling(max(Data,na.rm=TRUE))),col="black",las=1) #x-Achse
  axis(2,ylim=c(0,2),col="black",las=1) #y-Achse
  #drawnow;
}


par(par.defaults)
return(list(Pvalue=Pvalue,BinCenters=BinCenters,ObsNrInBin=ObsNrInBin,ExpectedNrInBin=ExpectedNrInBin,Chi2Value=Chi2Value))
}
