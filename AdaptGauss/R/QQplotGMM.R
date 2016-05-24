QQplotGMM=function(Data,Means,SDs,Weights,IsLogDistribution=Means*0,Line=TRUE,PlotSymbol=20,
                   xug=NULL,xog=NULL,LineWidth=2,PointWidth=0.8,
                   ylab='Data',main='QQ-plot Data vs GMM',...)
# QQplotGMM(Data,Means,SDs,Weights,IsLogDistribution,Line,PlotSymbol,xug,xog,LineWidth,PointWidth)
# Quantile/Quantile = QQ-Plot im Vergleich. zu einem Gauss Mixture Model oder LGL Model
# INPUT
# Data(1:n)	             Daten, deren Verteilung verglichen werden soll
# Means(1:L), SDs(1:L), Weights(1:L) die Paramter von Gaussians N(i) = Weights(i) * N(Means(i),SDs(i)
#                        die Gesamtverteilung ergibst sich als Summe der N(i)
# OPTIONAL
# IsLogDistribution(1:L) gibt an ob die Einzelverteilung einer (generalisierten)Lognormaverteilung ist
#                        wenn IsLogDistribution(i)==0 dann Mix(i) = Weights(i) * N(Means(i),SDs(i)
#                        wenn IsLogDistribution(i)==1 dann Mix(i) = Weights(i) * LogNormal(Means(i),SDs(i)
#                        Default: IsLogDistribution = Means*0;
# Line									Line in QQplot: =TRUE (Default), without False
# PlotSymbol             Symbol fur den qqplot, wenn nicht angegeben: PlotSymbol='b.'
# xug,xog                Grenzen der Interpolationsgeraden,  interpoliert wird fuer percentiles(x) in [xug,xog]
#                        Default: xug==min(x),xog==max(x), MT: Noch nicht implementiert!
# LineWidth              Linienbreite der Interpolationsgeraden; Default =3
# PointWidth             Dicke der Punkte im QQplot, existert nicht in Matlab
# LineSymbol             Liniensybol  der Interpolationsgerade;  Default ='r-'   MT: Nicht Implementiert
#
# in \dbt\Plot

# benutzt randomLogMix und qqplotfit
# MT 2014, reimplementiert aus Matlab von ALU 
# Aus historischen Gr?nden QQplotGMM MIT Ausgleichgerade

{
  

#xug = min(Data);  xog = max(Data); zu implementieren
# LineSymbol='r-' nicht implementiert
xlabel='Gaussian Mixture Model'

GMM = RandomLogGMM(Means,SDs,Weights,IsLogDistribution);

 quants<-qqplot(GMM, Data, pch=PlotSymbol, col="blue", cex=PointWidth, xlab=xlabel, ylab=ylab, main=main,...) #MT: na.rm=TRUE argument weglassen
 if(Line){
 fit<-lm(quants$y~quants$x)
 summary(fit)
 abline(fit, col="red", lwd=LineWidth)
 }
 return(invisible(quants))
}