Bayes4Mixtures <- function(Data, Means, SDs, Weights, IsLogDistribution = 0*Means, PlotIt = FALSE, CorrectBorders = FALSE){
# V = Bayes4Mixtures(Data, Means, SDs, Weights, IsLogDistribution, PlotIt, CorrectBorders)
# INPUT
# Data(1:N)            vector of data,  may contain NaN
# Means(1:C),SDs(1:C),Weights(1:C) parameters of the Gaussians (Mean, StdDeviation, Weight)
# 
# OPTIONAL
# IsLogDistribution(1:C) 1 oder 0, gibt an ob die jeweilige Verteilung eine Lognormaverteilung ist,(default ==0*(1:C))
# PlotIt              ==TRUE Verteilungen und Posteriors werden gezeichnet (default ==0);
# CorrectBorders      ==TRUE Daten an den Grenzen werden den Randverteilungen zugeordnet
#                    (default ==0) d.h. ganz gewoehnlicher Bayes mit allen seinen Problemen
# OUTPUT
# V	List of 2
# Posteriors(1:N,1:C)          Vektor der Posteriors      korrespondierend zu Data
# NormalizationFactor(1:N)     Nenner des Bayes Theorems  korrespondierend zu Data
# 
# 
# AUTHOR: CL
# 1.Editor: MT 08/2015 : Plotten neu, Variablen vereinheitlicht
  
# Sortiere die Daten und merke die unsortiert-Reihenfolge
AnzMixtures <- length(Means)
Kernels <- unique(Data)
UNsortInd <- match(Data, Kernels)
AnzKernels <- length(Kernels)

# Berechne bedingte Wkeit p(x|ci) mit x = Daten[1:N] und ci = Klasse i, i aus 1 bis C
PDataGivenClass <- matrix(0,AnzKernels,AnzMixtures);
for(i in c(1:AnzMixtures)){
	if( IsLogDistribution[i] == 1 ){ # LogNormal
		PDataGivenClass[,i] <- symlognpdf(Kernels,Means[i],SDs[i]); # LogNormaldichte. Achtung: Muss gespiegelt werden fuer negative Werte.
	}else{
		PDataGivenClass[,i] <- dnorm(Kernels,Means[i],SDs[i]); # Gaussian
	}#end if(IsLogDistribution[i] == 1)
}#end for(i in c(1:AnzMixtures))

NormalizationFactor <- PDataGivenClass %*% Weights;  # Gewichtete Summe der Priors; 
# Zum Normalisizerungsfaktor:
# Achtung: Es soll 1.Spalte mal 1.Eintrag von Weights + 2.Spalte mal 2.Eintrag von Weights usw. gerechnet werden. 
# Dazu brauchen wir Matrixmultiplikation!
# Bei PDataGivenClass * Weights wird die 1.Zeile von P... mal 1. Wert von Weights, 2.Zeile von P mal 2. Wert usw.
# gerechnet, was nicht Sinn der Sache ist!!!

Pmax = max(NormalizationFactor);

# to prevent division error:
ZeroInd <- which(NormalizationFactor==0);
if(length(ZeroInd) > 0){
	NormalizationFactor[ZeroInd] =10^(-7)
}#end if(length(ZeroInd) > 0) 

#Posterioris nach Bayes p(c|x) mit c = Klasse (ueber-, unter- oder nicht exprimiert) und x Datensatz.
PClassGivenData <- matrix(0, AnzKernels, AnzMixtures);
for(i in c(1:AnzMixtures)){
	PClassGivenData[,i] <- PDataGivenClass[,i]*Weights[i] / NormalizationFactor
}#end for(i in c(1:AnzMixtures))

if(CorrectBorders == TRUE & (sum(IsLogDistribution)==0)){ # Randkorrektuuren anbringen
	# Daten kleiner kleinstem Modus werden diesem zugeschlagen; 
	KleinsterModus <- min(Means)
	SmallModInd <- which.min(Means)
	LowerInd <- which(Kernels<KleinsterModus);
	for(i in c(1:AnzMixtures)){
		PClassGivenData[LowerInd,i] <- 0;
	}#end for(i in c(1:AnzMixtures))
	PClassGivenData[LowerInd,SmallModInd] <- 1;
	# Daten groesser groesstem Modus werden diesem zugeschlagen; 
	GroessterModus <- max(Means)
	BigModInd <- which.max(Means)
	HigherInd <- which(Kernels>GroessterModus);
	for(i in c(1:AnzMixtures)){
		PClassGivenData[HigherInd,i] <- 0;
	}#end for(i in c(1:AnzMixtures))
	PClassGivenData[HigherInd,BigModInd] <- 1;
}#end if(CorrectBorders == TRUE & (sum(IsLogDistribution)==0))

# jetzt zurueck auf Daten
Posteriors = zeros(length(Data),AnzMixtures);
for(i in c(1:AnzMixtures)){
	Posteriors[,i] <- PClassGivenData[UNsortInd,i];  # fuer die Daten anpassen
}#end for(i in c(1:AnzMixtures))

# auch noch den Normalisierungsfaktor auf Datengroesse anpassen
Nenner <- NormalizationFactor;
NormalizationFactor <- NormalizationFactor[UNsortInd];

## MT: Neu gemacht
if(PlotIt==TRUE){
	color <- rainbow(AnzMixtures)
  xlim=c(min(Data),max(Data))
  ylim=c(0,1.05)
	plot.new()
	par(xaxs='i')
	par(yaxs='i')
	par(usr=c(xlim,ylim))
	ind=order(Data)
  if(CorrectBorders){

    for(i in 1:AnzMixtures){
      points(Data[ind], Posteriors[ind,i], col = color[i],type='l',lwd=2)
    }#end for(i in 2:AnzMixtures)   
  }else{
	for(i in 1:AnzMixtures){
		points(Data[ind], Posteriors[ind,i], col = color[i],type='l',lwd=2)
	}#end for(i in 2:AnzMixtures)
  }
}#end if(PlotIt==TRUE)
axis(1,xlim=xlim,col="black",las=1) #x-Achse
axis(2,ylim=ylim,col="black",las=1) #y-Achse
#box() #Kasten um Graphen
title(ylab='posteriori',xlab='Data')
##
res <- list(Posteriors = Posteriors, NormalizationFactor=NormalizationFactor)
return (res) 
 

# symlognpdf
#########################################################
symlognpdf <- function(Data,Means,SDs){
  #pdf = symlognpdf(Data,Means,SDs);
  # for Means>0 same as dlnorm(Data,Means,SDs); (Dichte der log-Normalverteilung)
  # for Means < 0: mirrored at y axis
  #INPUT
  #Data[1:n]  x-values
  #Means,SDs        Mean and Sdev of lognormal
  
  temp<-symlognSigmaMue(Means,SDs)
  mu<-temp$mu
  sig<-temp$sig
  if(Means>=0){
    pdfkt<-dlnorm(Data,meanlog=mu,sdlog=sig)  
  }else{
    pdfkt<-Data*0
    negDataInd<-which(Data<0)
    pdfkt[negDataInd] <- dlnorm(-Data[negDataInd],meanlog=mu,sdlog=sig)
    plot(Data,pdfkt)
  }
  return (pdfkt) 
  
  symlognSigmaMue <-  function(Means,SDs){
    
    variance<-log(SDs*SDs/(Means*Means)+1)
    sig<-sqrt(variance)
    mu<-log(abs(Means))-0.5*variance
    return (list(variance=variance,sig=sig,mu=mu)) 
    
  }
  
}
#########################################################

# zeros
#########################################################
zeros <-function (n,m=n,o=0) {
  # zeros(n)     returns an n-by-n matrix of 0s. 
  # zeros(n,1)   returns a vector of 0s 
  # zeros(n,m)   returns an n-by-m matrix of zeros
  # zeros(n,m,o) returns an 3D matrix  of zeros
  
  # ALU
  
  if (m==1) { # vector wird zurueckgegeben
    return(c(1:n)*0) ;   
  }else{      # return n-by-m matrix of ones.         
    if  (o==0){
      return(matrix(0,n,m));
    }else{   #3D matrix
      nullen = rep(0, m*n*o);  # soviele nullen wie in die 3D matrix pasen
      return(array(nullen,c(n,m,o)));
      
    } # end  if  (o==0)
  } # end if (m==1) 
} # end  function  zeros
#########################################################

}