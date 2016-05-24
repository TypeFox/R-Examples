plot.ContObservHMM <-
function(x, Series=x$Observations, ylabel="Observation series", xlabel="Time", ...)
{
	plot(Series, type = "l", col="dodgerblue4", ylab=ylabel, xlab=xlabel, lwd=1.5)
	nq<-x$Viterbi
	nq<-ifelse((x$Viterbi)==1, min(Series), max(Series))
	lines(nq, col="firebrick", lwd=1.5)
}
