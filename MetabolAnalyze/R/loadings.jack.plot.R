loadings.jack.plot <-
function(output)
{
	SignifHighppm<-abs(output$SignifW[,1])>output$cutoff
    SignifHighW<-matrix(output$SignifW[SignifHighppm,], ncol=output$q)
	LowerH<-matrix(output$Lower[SignifHighppm,], ncol=output$q)
	UpperH<-matrix(output$Upper[SignifHighppm,], ncol=output$q)

	## Bar plot of significant high loadings
	barplot2(SignifHighW[,1], ylim=c(min(LowerH[,1],0)-1.5, max(UpperH[,1],0)+1.5), las=2, width=0.5, space=0.5, plot.grid=TRUE, ylab="PC 1 loadings", xlab="Spectral regions", names.arg = rownames(output$SignifW)[SignifHighppm], plot.ci = TRUE,
ci.l = LowerH[,1], ci.u = UpperH[,1], font.main = 2, col="red")
} # close plot.loadings.jack

