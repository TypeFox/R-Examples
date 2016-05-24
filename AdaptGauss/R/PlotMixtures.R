PlotMixtures <- function(Data, Means, SDs, Weights = rep(1/length(Means),length(Means)), IsLogDistribution = rep(FALSE,length(Means)),SingleColor = 'blue', MixtureColor = 'red',DataColor='black',SingleGausses=FALSE,axes=TRUE,xlab, ylab,xlim, ylim,  ...){
# PlotMixtures(Data,Means,SDs,Weights,IsLogDistribution,SingleColor,MixtureColor);
# PlotMixtures(Data,Means,SDs,Weights,IsLogDistribution);
# Plot a Mixture of Gaussians

# INPUT
# Data(1:AnzCases)            this data column for which the distribution was modelled
# Means(1:L,1)                Means of Gaussians,  L ==  Number of Gaussians
# SDs(1:L,1)                estimated Gaussian Kernels = standard deviations
# OPTIONAL
# Weights(1:L,1)              relative number of points in Gaussians (prior probabilities): sum(Weights) ==1, default 1/L
# IsLogDistribution(1:L)      default ==0*(1:L), gibt an ob die jeweilige Verteilung eine Lognormalverteilung ist
#
# SingleColor                 PlotSymbol of all the single gaussians, default magenta
# MixtureColor                PlotSymbol of the mixture default black
#
# SingleGausses               Sollen die einzelnen Gauss auch gezeichnet werden, dann TRUE
# ...							            other plot arguments like xlim = c(1,10)
#
# Author: MT 08/2015, 
# 1.Editor: MT 1/2016: PDF4Mixtures als eigene Funktion ausgelagert
#Nota: Based on a Version of HeSa Feb14 (reimplemented from ALUs matlab version)      

#oldpar <- par(no.readonly = TRUE)
# labels
if(missing(xlab)){
	xlab = '' # no lable for x axis
}
if(missing(ylab)){
	ylab = 'PDE' # no lable for y axis
}
if(missing(IsLogDistribution)) IsLogDistribution = rep(FALSE,length(Means))
if(length(IsLogDistribution)!=length(Means)){
  warning(paste('Length of Means',length(Means),'does not equal length of IsLogDistribution',length(IsLogDistribution),'Generating new IsLogDistribution'))  
  IsLogDistribution = rep(FALSE,length(Means))
}
X = sort(unique(Data)) # sort ascending and make sure of uniqueness
AnzGaussians <- length(Means)
#SingleGaussian <- matrix(0,length(X),AnzGaussians)
#GaussMixture=X*0 # init
if(SingleColor  != 0){
# 	for(g in c(1:AnzGaussians)){
# 		if(IsLogDistribution[g] == TRUE){ # LogNormal 
# 			SingleGaussian[,g] <- symlognpdf(X,Means[g],SDs[g])*Weights[g] # LogNormal 
# 		}else{ # Gaussian
# 			SingleGaussian[,g] = dnorm(X,Means[g],SDs[g])*Weights[g]
# 		}# if IsLogDistribution(i) ==T  
# 		GaussMixture =  GaussMixture + SingleGaussian[,g]
# 	} # for g
	pdfV=Pdf4Mixtures(X, Means, SDs, Weights)
	SingleGaussian=pdfV$PDF
	GaussMixture=pdfV$PDFmixture
	# Limits
  GaussMixtureInd=which(GaussMixture>0.00001)
	if(missing(xlim)){ # if no limits for x-axis are comitted
    xl = max(X[GaussMixtureInd],na.rm=T)
	  xa = min(X[GaussMixtureInd], 0,na.rm=T)
		xlim = c(xa,xl)
	}
# Je plot xlim und ylim uebergeben
# Falls dies nicht geschieht, werden beide Achsen falsch skaliert!
	if(missing(ylim)){ # if no limits for y-axis are comitted
    yl <- max(GaussMixture,na.rm=T) 
	  ya <- min(GaussMixture, 0,na.rm=T)
		ylim <- c(ya,yl+0.1*yl)
     #ylim= par("yaxp")[1:2]
	}
  plot.new()
  par(xaxs='i')
  par(yaxs='i')
  par(usr=c(xlim,ylim))
  if(SingleGausses){
    	for(g in c(1:AnzGaussians)){
    		par(new = TRUE)
    		points(X, SingleGaussian[,g], col = SingleColor, type = 'l', xlim = xlim, ylim = ylim, ...)
    	}
    }
   points(X, GaussMixture, col = MixtureColor, type = 'l',xlim = xlim,...)
} else{#SingleColor  == 0
  plot(X, GaussMixture, col = MixtureColor, type = 'l', xlim, ylim,axes=FALSE, ...)
}
pareto_radius<-ParetoRadius(Data)
pdeVal        <- ParetoDensityEstimation(Data,pareto_radius)
points(pdeVal$kernels,pdeVal$paretoDensity,type='l', xlim = xlim, ylim = ylim,col=DataColor,...)   


if(axes){
  axis(1,xlim=xlim,col="black",las=1) #x-Achse
  axis(2,ylim=ylim,col="black",las=1) #y-Achse
  title(xlab=xlab,ylab=ylab)
  #box() #Kasten um Graphen
  #title(xlab=xlab,ylab=ylab)
}
#par(oldpar) # geht nicht, da sonst zB. abline( v =  0.4) nicht bei 0.4 sitzt...

#return(list(X,GaussMixture))
}# end PlotMixtures

