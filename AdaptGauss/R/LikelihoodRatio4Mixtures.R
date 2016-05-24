LikelihoodRatio4Mixtures <- function(Data,NullMixture,OneMixture,PlotIt=FALSE, LowerLimit = min(Data, na.rm=TRUE), UpperLimit = max(Data, na.rm=TRUE)){
# V = LikelihoodRatio4Mixtures(Data,NullModel,OneModel,PlotIt,LowerLimit,UpperLimit)
# Likelihood Ratio Test for 2 Gauss mixture Models
#
# INPUT
# Data[1:N]                                       data points
# NullMixture 																		=cbind(M0,S0,W0) or  cbind(M0,S0,W0,IsLog0); 
#																									the null model usually with less Gaussians than the OneMixture
# OneMixture  																		=cbind(M1,S2,W1) or  cbind(M1,S1,W1,IsLog1);
#																									the alternative model usually with more Gaussians than the OneMixture
#
# OPTIONAL
# PlotIt                           do a Plot of the compared cdf's and the KS-test distribution (Diff)
# LowerLimit                       test only for Data >= LowerLimit, Default = min(Data) i.e all Data.
# UpperLimit                       test only for Data <= UpperLimit, Default = max(Data) i.e all Data.
#
# OUTPUT
# List of 3:
# Pvalue                           the error that we make, if we accept OneMixture as the better Model over the NullMixture
# NullLogLikelihood                log likelihood of GMM Null
# OneLogLikelihood                 log likelihood of GMM One

# uses LogLikelihood4Mixtures in that  PdfForMixes
# ALU 2015
# Uebertrag von Matlab nach R: CL 02/2016
#1. Editor: MT 03/2016
	  #library(grid)# fuer multiplot
ValidDataInd = which((Data >= LowerLimit) & (Data <= UpperLimit))

## NullMixture
NullAnzMixes <- dim(NullMixture)[1]
AnzCol <- dim(NullMixture)[2]
M0 = NullMixture[,1] 
S0 = NullMixture[,2]
W0 = NullMixture[,3] 
if(AnzCol >3){ # IsLog0 gegeben;
    IsLog0 = NullMixture[,4] 
}else{
   IsLog0 = M0*0; 
}#end if IsLog0 gegeben

## OneMixture
OneAnzMixes <- dim(OneMixture)[1]
AnzCol <- dim(OneMixture)[2]
M1 = OneMixture[,1] 
S1 = OneMixture[,2] 
W1 = OneMixture[,3] 
if(AnzCol >3){ # IsLog1 gegeben;
   IsLog1 = OneMixture[,4] 
}else{
   IsLog1 = M1 *0; 
}#end if IsLog1 gegeben

# die LogLikelihood der beiden Mixturen ausrechnen
NullLL <- LogLikelihood4Mixtures(Data[ValidDataInd],M0,S0,W0,IsLog0)
NullLogLikelihood <- NullLL$LogLikelihood
NullLogPDF <- NullLL$LogPDF
NullPDF <- NullLL$PDFmixture

OneLL <- LogLikelihood4Mixtures(Data[ValidDataInd],M1,S1,W1,IsLog1)
OneLogLikelihood <- OneLL$LogLikelihood
OneLogPDF <- OneLL$LogPDF
OnePDF <-  OneLL$PDFmixture

D <- 2*(OneLogLikelihood - NullLogLikelihood)
DegreesOfFreedom <- abs(length(M1)*3-length(M0)*3)
# Wkeitsverteilung von D folgt naeherungsweise einer Chi-Quadrat-Verteilung mit DegreesOfFreedom-vielen Frei-
# heitsgraden.

Pvalue = 1- pchisq(D, DegreesOfFreedom) # cdf der Chi-sq-Verteilung an Teststatistik-Wert auswerten.

if(PlotIt ==TRUE){
   y=NULL
    x=NULL
    kernels=NULL
   paretoRadius<-ParetoRadius(Data[ValidDataInd])
   pdeVal        <- ParetoDensityEstimation(Data[ValidDataInd],paretoRadius,NULL)
   paretoDensity <- pdeVal$paretoDensity*1
   dfframe = as.data.frame(cbind(kernels = pdeVal$kernels, density = paretoDensity))
   #colnames(dfframe)=c('kernels','density')

    #  kernels=pdeVal$kernels
   pdePlot <-ggplot()+ geom_line(data = dfframe, aes(x = kernels, y = density), colour = "black")
   
   nullmodell=as.data.frame(cbind(kernels = Data[ValidDataInd],density = NullPDF[ValidDataInd]))
	 onemodell=as.data.frame(cbind(kernels = Data[ValidDataInd],density = OnePDF[ValidDataInd]))
	pdePlot <- pdePlot + geom_line(data = nullmodell, aes(x=kernels, y=density), colour="blue")
	pdePlot <- pdePlot + geom_line(data = onemodell, aes(x=kernels, y=density), colour="red")
	pdePlot <- pdePlot + ggtitle("PDE and Mixture Models")
	pdePlot <- pdePlot + ylab("blue = Nullmodel, red = Onemodel")
	pdePlot <- pdePlot + xlab("Data")
	
	nulllogmodell=as.data.frame(cbind(kernels = Data[ValidDataInd],density = NullLogPDF[ValidDataInd]))
	onelogmodell=as.data.frame(cbind(kernels = Data[ValidDataInd],density = OneLogPDF[ValidDataInd]))
	# Erstelle plot von LogPDF beider Modelle
	LogPDFplot = ggplot() + geom_line(data = nulllogmodell, aes(x=kernels, y=density), colour="blue")
	LogPDFplot =	LogPDFplot + geom_line(data = onelogmodell, aes(x=kernels, y=density), colour="red")
	LogPDFplot <- LogPDFplot + ggtitle(paste0("LogPDF, LogLikeli Null= ",NullLogLikelihood, ", One = ", OneLogLikelihood))
	LogPDFplot <- LogPDFplot + ylab("blue = Nullmodel, red = Onemodel")
	LogPDFplot <- LogPDFplot + xlab("Data")
	nuller=as.data.frame(cbind(x = c(D,D), y = c(0,pchisq(D, DegreesOfFreedom))))
	einser=as.data.frame(cbind(x = c(0,D), y = c(pchisq(D, DegreesOfFreedom),pchisq(D, DegreesOfFreedom))))
	# Erstelle plot der CDF der Chi-quadrat-Verteilung und markiere den Wert der Teststatistik
	x11 = data.frame(c(0, 10))
	#y0=c(0,pchisq(D, DegreesOfFreedom))
	#y1=c(pchisq(D, DegreesOfFreedom),pchisq(D, DegreesOfFreedom))
	ChisqCDF <- ggplot(x, aes(x11)) + stat_function(fun = pchisq, args=list(df=DegreesOfFreedom))
	ChisqCDF <- ChisqCDF + geom_line(data= nuller, aes(x,y) , colour = 'blue')
	ChisqCDF <- ChisqCDF + geom_line(data= einser, aes(x,y, color ='Myline') , colour = 'blue')
	ChisqCDF <- ChisqCDF + ggtitle(paste0('Chi^2 CDF with DF = ', DegreesOfFreedom))
	ChisqCDF <- ChisqCDF + ylab("blue = value of test statistic")
	ChisqCDF <- ChisqCDF + xlab(paste0("Pvalue = ", Pvalue))
	nullminuseins=as.data.frame(cbind(kernels = Data[ValidDataInd],density = OneLogPDF[ValidDataInd]-NullLogPDF[ValidDataInd]))
	# Erstelle plot von LogPDF des OneModels minus LogPDF des NullModels
	LogPDFOneMinusNullplot <- ggplot() + geom_line(data = nullminuseins, aes(x=kernels, y=density), colour="black")
	LogPDFOneMinusNullplot <- LogPDFOneMinusNullplot + ggtitle('LogPDF OneModel - NullModel')
	LogPDFOneMinusNullplot <- LogPDFOneMinusNullplot + ylab("LogPDF OneModel - LogPDF NullModel")
	LogPDFOneMinusNullplot <- LogPDFOneMinusNullplot + xlab("Data")
	# Multiple plot function
	#
	# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
	# - cols:   Number of columns in layout
	# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
	#
	# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
	# then plot 1 will go in the upper left, 2 will go in the upper right, and
	# 3 will go all the way across the bottom.
	#
	#author: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
	multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

	  
	  # Make a list from the ... arguments and plotlist
	  plots <- c(list(...), plotlist)
	  
	  numPlots = length(plots)
	  
	  # If layout is NULL, then use 'cols' to determine layout
	  if (is.null(layout)) {
	    # Make the panel
	    # ncol: Number of columns of plots
	    # nrow: Number of rows needed, calculated from # of cols
	    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
	                     ncol = cols, nrow = ceiling(numPlots/cols))
	  }
	  
	  if (numPlots==1) {
	    print(plots[[1]])
	    
	  } else {
	    # Set up the page
	    grid.newpage()
	    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
	    
	    # Make each plot, in the correct location
	    for (i in 1:numPlots) {
	      # Get the i,j matrix positions of the regions that contain this subplot
	      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
	      
	      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
	                                      layout.pos.col = matchidx$col))
	    }
	  }
	}
	# plotte alle 4 Zeichnungen in einer Uebersicht:	
	multiplot(pdePlot,LogPDFplot, ChisqCDF, LogPDFOneMinusNullplot, layout=matrix(c(1,2,3,4), nrow=2, byrow=TRUE))
}
return(list(Pvalue=Pvalue, NullLogLikelihood=NullLogLikelihood, OneLogLikelihood=OneLogLikelihood))
}#end function



