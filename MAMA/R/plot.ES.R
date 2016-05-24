#Ine moznosti grafickeho zobrazenia vysledkov metody kombinovania velkosti efektu
plotES<-function(theScores, ScoresFDR, num.studies, legend.names, colors, which)
{

 par(mfrow=c(1,length(which)))
if (1 %in% which) {
 IDRplot(theScores, CombineExp = 1:num.studies, colPos = colors[2], 
 colNeg = colors[3])
# savePlot("ES2 IDRplot.png", type="png")
}
if (2 %in% which) {
 FDRwholeSettwo <- sort(ScoresFDR$two.sided[, "FDR"])
 experimentstwo <- list()
 for (j in 1:num.studies) 
 {experimentstwo[[j]] <- sort(ScoresFDR$two.sided[, 
 paste("FDR_Ex_", j, sep = "")])}
 theNewC<-colors[2:length(colors)]
 plot(FDRwholeSettwo, pch = ".", col = colors[1], ylab = "FDR", 
 xlab = "Number of genes")
 for (j in 1:num.studies) points(experimentstwo[[j]], pch = ".", 
 col = theNewC[j])
 legend(x="bottomright", legend=legend.names, 
 col=c(colors[1], theNewC[1:num.studies]), pch=19)
# savePlot("ES3 FDR.png", type="png") 
}
if (3 %in% which) {
 CountPlot(ScoresFDR, Score = "FDR", kindof = "two.sided", 
 cols = colors, main = "two sided FDR", 
 xlab = "FDR threshold", ylab = "Number of genes")
 legend(x=0.07,y=600, legend=legend.names, col=colors, pch="*")
# savePlot("ES4 countplot.png",type="png") 
}
}
legend.names=c("Combined set", "Denmark", "Australia", "Japan")
colors=c("red","blue","green","yellow")