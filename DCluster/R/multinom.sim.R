multinom.sim<-function(data, mle=NULL)
{
	if(is.null(mle))
		data$Observed<-rmultin(sum(data$Observed), data$Expected/sum(data$Expected))
	else
		data$Observed<-rmultin(mle$n, mle$p)
	return(data)
}
