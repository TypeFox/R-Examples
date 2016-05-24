KStestMixtures=function(Data,Means,SDs,Weights,IsLogDistribution=Means*0,PlotIt=FALSE,UpperLimit=max(Data,na.rm=TRUE)){
# res= KStestMixtures(Data,Means,SDs,Weights,IsLogDistribution,PlotIt,UpperLimit)
# Kolmogorov-Smirnov Test Data vs a given Gauss Mixture Model
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
#
# OUTPUT
# PvalueKS                         Pvalue of a suiting Kolmogorov Smirnov Test, PvalueKS==0 if PvalueKS<0.001
# DataKernels,DataCDF              such that plot(DataKernels,DataCDF) gives the cdf(Data)
# CDFGaussMixture                  such that plot(DataKernels,DataCDF) gives the cdf(Data)

# MT 2015, reimplemented from ALUs matlab version

################################
#Hilfsfunktionen
################################

################################  
################################  
#[DataCDF,DataKernels] = ecdfUnique(Data) 
dummy <- ecdf(Data) # cdf(Diff)
DataCDF <- c(0,get("y", envir = environment(dummy)))
DataKernels <- c(knots(dummy)[1],knots(dummy))#CDFuniq# cdf(Data)

CDFGaussMixture = CDFMixtures(DataKernels,Means,SDs,Weights,IsLogDistribution)$CDFGaussMixture # cdf(GMM)


# UpperLimit bereucksichtigen, alles auf <=UpperLimit kuerzen
Data=Data[Data<=UpperLimit]
Ind = which(DataKernels<=UpperLimit, arr.ind=T)
DataKernels    = DataKernels[Ind]
DataCDF        = DataCDF[Ind]
CDFGaussMixture= CDFGaussMixture[Ind]


# KS messgroesse bestimmen
CDFdiff = abs(CDFGaussMixture-DataCDF)         # Unterschied = KS- testgroesse
MaxDiff= max(CDFdiff,na.rm=TRUE)              # wo steckt der groesste Unterschied
MaxInd=which(CDFdiff==MaxDiff, arr.ind=T)
KernelMaxDiff   = DataKernels[MaxInd]         # wo steckt der groesste Unterschied
DataCDFmaxDiff  = DataCDF[MaxInd]             # wo steckt der groesste Unterschied
GMMCDFmaxDiff   = CDFGaussMixture[MaxInd]      # wo steckt der groesste Unterschied

         
             
# Die Miller Funktion via Monte-Carlo errechnen
AnzData =length(Data)
AnzRepetitions = 1000
if(AnzData<1000) AnzRepetitions = 2000
if(AnzData<100)  AnzRepetitions = 5000
    
RandGMMDataDiff = matrix(0,AnzRepetitions,1)

for(i in c(1:AnzRepetitions)){
   R = RandomLogGMM(Means,SDs,Weights,IsLogDistribution,AnzData)
   #[RandCDF,RandKernels] = ecdfUnique(R)
   dummy <- ecdf(R) # cdf(Diff)
   RandCDF <- c(0,get("y", envir = environment(dummy)))
   RandKernels <- c(knots(dummy)[1],knots(dummy))#CDFuniq# cdf(Data)
   #
   #RandCDF =  interp1(RandKernels,RandCDF,DataKernels, 'linear') #matlab
   RandCDF = approx(RandKernels,RandCDF,DataKernels, 'linear')$y  # den MaxDiff in cdf(Diff) lokalisieren
   RandGMMDataDiff[i] = max(abs(CDFGaussMixture-RandCDF),na.rm=TRUE)
}# for i
AllDiff =  RandGMMDataDiff                     # die Verteilung aller Differenzen
#[AllDiffCDF,AllDiffKenels] =  ecdfUnique(AllDiff);  # cdf(Diff) #matlab
dummy <- ecdf(AllDiff) # cdf(AllDiff)
AllDiffCDF <- c(0,get("y", envir = environment(dummy)))
AllDiffKernels <- c(knots(dummy)[1],knots(dummy))#CDFuniq# cdf(Data)

#matlab
#MaxDiffCDFwert = interp1([0;AllDiffKenels],[0;AllDiffCDF],MaxDiff, 'linear');  # den MaxDiff in cdf(Diff) lokalisieren
MaxDiffCDFwert = approx(rbind(0,unname(AllDiffKernels)),rbind(0,unname(AllDiffCDF)),MaxDiff, 'linear')$y;  # den MaxDiff in cdf(Diff) lokalisieren

PvalueKS = MaxDiffCDFwert                        # P- value fuer KS-test ausrechnen
PvalueKS = round(PvalueKS,3)                     # runden auf 3 gueltige stellen


if(PlotIt ==1){
  par.defaults <- par(no.readonly=TRUE)
  par(mfrow = c(2,2))

  xlim=c(min(DataKernels), max(DataKernels))
  #ylim=c(min(min(DataCDF,na.rm = TRUE),min(CDFGaussMixture,na.rm = TRUE)),max(max(DataCDF,na.rm = TRUE),max(CDFGaussMixture,na.rm = TRUE)))
  ylim=c(0,1)
  plot(DataKernels,DataCDF,col='blue',xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE)
  points(DataKernels,CDFGaussMixture,col='red',type='l',lwd=2)
  
  points(KernelMaxDiff,DataCDFmaxDiff,col='green',pch=1)
  points(KernelMaxDiff,GMMCDFmaxDiff,col='green',pch=1)
  abline(v=KernelMaxDiff,col='green')
  axis(1,xlim=xlim,col="black",las=1) #x-Achse
  axis(2,ylim=ylim,col="black",las=1) #y-Achse
  title('KS-test: comparison CDF(GMM) vs CDF(Data)',xlab='Data',ylab='red =CDF(GMM), blue=CDF(Data)');
  
  #     subplot(2,2,2);   # hier die fraglichen PDFs zeichnen
  PlotMixtures(Data,Means,SDs,Weights,IsLogDistribution=IsLogDistribution,xlim=xlim,ylim=ylim,xlab='',ylab='',SingleGausses = T,SingleColor='magenta',MixtureColor='blue')
  title(paste0('max(Diff) at: ',KernelMaxDiff),xlab='Data',ylab='pdf(GMM), red= pdf(Data)')
  abline(v=KernelMaxDiff,col='green')
#     subplot(2,2,3);
  		pareto_radius2<-ParetoRadius(AllDiff) 
			pdeVal2        <- ParetoDensityEstimation(AllDiff,pareto_radius2)
			xlim=c(min(MaxDiff*0.90,pdeVal2$kernels),max(pdeVal2$kernels))
			plot(pdeVal2$kernels,pdeVal2$paretoDensity,type='l',xaxs='i',
			yaxs='i',xlab='MaxDiff, mag = max(Diff)',ylab='pdf(KS-MaxDiff)',main='KS-Distribution of MaxDiff',xlim=xlim,col='blue') 
      abline(v=MaxDiff,col='magenta')
      box()
#     subplot(2,2,4);
      xlim=c(min(MaxDiff*0.90,AllDiffKernels),max(AllDiffKernels))
      ylim=c(min(MaxDiffCDFwert*0.90,AllDiffCDF),max(AllDiffCDF))
      plot(AllDiffKernels,AllDiffCDF,type='p',ylab='cdf(KS-MaxDiff)',xlab='Diff, mag = max(Diff)',main=   paste0('Pvalue: ',PvalueKS*100,' [#]'),ylim=ylim,xlim=xlim,col='blue')
      abline(v=MaxDiff,col='magenta')
      #abline(a=c(MaxDiffCDFwert,MaxDiffCDFwert),b=c(0,MaxDiff))
      abline(h=MaxDiffCDFwert,col='magenta')
      box()
      par(par.defaults)
} #if PlotIt ==1
    
return(list(PvalueKS=PvalueKS,DataKernels=DataKernels,DataCDF=DataCDF,CDFGaussMixture=CDFGaussMixture))
}


