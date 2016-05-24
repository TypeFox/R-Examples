LogLikelihood4Mixtures <- function(Data, Means, SDs, Weights, IsLogDistribution=Means*0){
# LogLikelihood <- LogLikelihood4Mixtures(Data,Means,SDs,Weights,IsLogDistribution)
# berechnung der Loglikelihood fuer ein Mixture model: LogLikelihood = sum(log(PDFmixture))
#
# INPUT
# Data[1:n]									Daten, deren Verteilung verglichen werden soll
# Means[1:L]								Means of Gaussians,  L ==  Number of Gaussians
# SDs[1:L]								estimated Gaussian Kernels = standard deviations
# Weights[1:L]							relative number of points in Gaussians (prior probabilities): sum(Weights) ==1
#
# OPTIONAL
# IsLogDistribution[1:L] 		gibt an, ob die Einzelverteilung einer (generalisierten)Lognormaverteilung ist
#                        		wenn IsLogDistribution[i]==0 dann Mix(i) = W[i] * N(M[i],S[i])
#                        		wenn IsLogDistribution[i]==1 dann Mix(i) = W[i] * LogNormal(M[i],S[i])
#                        		Default: IsLogDistribution = Means*0;
#
# OUTPUT
# LogLikelihood           die Loglikelihood der Verteilung = LogLikelihood = = sum(log(PDFmixture))
# LogPDF(1:n)             = log(PDFmixture);  
# PDFmixture              die Probability density function an jedem Datenpunkt

#	Author: ALU, 2015
#	Uebertrag von Matlab nach R: CL 02/2016
# 1.Editor: MT 02/2016: umbenannt in LogLikelihood4Mixture, da auch LGL Modelle moegliech und analog zu LikelihoodRatio4Mixtures, Chi2testMixtures, KStestMixtures

#Pattern Recogintion and Machine Learning, C.M. Bishop, 2006, isbn: ISBN-13: 978-0387-31073-2, p. 433 (9.14)

PdfForMix = Pdf4Mixtures(Data,Means,SDs,Weights,IsLogDistribution) # PDF ausrechnen
PDFmixture <- PdfForMix$PDFmixture
PDFmixture[PDFmixture<=0] = NaN            # null zu NaN
LogPDF = log(PDFmixture)                   # logarithmieren (natuerlicher Logarithmus)
LogLikelihood = sum(LogPDF, na.rm=TRUE)	# summieren
return(list(LogLikelihood=LogLikelihood, LogPDF = LogPDF, PDFmixture = PDFmixture))
}#end function