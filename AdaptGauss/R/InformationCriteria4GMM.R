InformationCriteria4GMM=function(Data,Means,SDs,Weights,IsLogDistribution=Means*0){
# Vres = InformationCriteria4GMM(Data,M,S,W,IsLogDistribution)
# berechnung von AIC (Akaike Information Criterium) und BIC (Bayes Information Criterium)
# fuer GMM in einer Variablen
#
# INPUT
#Data[1:N] 			Vector (1:N) of data points; May Contain NaN
#Means[1:L]			vector[1:L] of Means of Gaussians (of GMM),L ==  Number of #							 Gaussians
#SDs[1:L]			  vector of standard deviations, estimated Gaussian Kernels, #							 has to be the same length as Means
#Weights[1:L]		vector of relative number of points in Gaussians (prior #								probabilities), has to be the same length as Means
# 							sum(Weights)=1
# OPTIONAL
#IsLogDistribution[1:L]		Optional, ==1 if distribution(i) is a LogNormal, #										 			default  vector of zeros of length L}
#
# OUTPUT
# K               Number of gaussian mixtures
# AIC             Akaike Informations criterium
# BIC             Bayes Information criterium
# LogLikelihood   LogLikelihood of GMM, see 
#									\code{\link{LogLikelihood4Mixtures}}
# PDFmixture(1:N) probability density function of GMM, see 
#									\code{\link{Pdf4Mixtures}}
# LogPDFdata(1:N) log(PDFmixture)
# author: MT 01/2016  
# Source
# Aho, K., Derryberry, D., & Peterson, T. (2014). Model selection for ecologists: the worldviews of AIC and BIC. Ecology, 95(3), 631-636.
 
#Note: Variablenbezeichnungen mit package AdaptGauss vereinheitlicht.



#   pdfV=Pdf4Mixtures(Data,Means,SDs,Weights,IsLogDistribution,PlotIt=F)
#   PDFmixture <- PdfForMix$PDF
#   pdf=colSums(PDFmixture)
#   pdf[<=0]=NaN
# LogLikelylihood=sum(log(pdf),na.rm=T)
  kkk=length(Means)
  if(kkk>1) #Normalfall
    K=kkk*3
  else #Bei einem Gauss gibt es kein gewicht!
    K=2
  
  
  LL4GMM=LogLikelihood4Mixtures(Data,Means,SDs,Weights,IsLogDistribution)
  LogLikelihood=LL4GMM$LogLikelihood
  N=length(Data)
  AIC=2*K-2*LogLikelihood
  if(N<40)
    AIC= AIC + (2*K*(K+1))/(N-K-1) 
  #AIC=LogLikelihood-K #S. Bishop 2006, p. 33 (1.73)  
  #BIC=LogLikelihood-1/2*K*log(length(Data))#S. Bishop 2006, p 217 (4.139) 
  BIC= K*log(N)-2*LogLikelihood
	
  return(list(K=K, AIC=AIC, BIC=BIC,LogLikelihood =LogLikelihood, PDFmixture=LL4GMM$PDFmixture,LogPDFdata=LL4GMM$LogPDF))
}