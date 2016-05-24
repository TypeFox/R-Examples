standardize <-
function(Covars)
{
	for(i in 1:ncol(Covars))
	{
	 rg<-range(Covars[,i])
     Covars[,i]<-(Covars[,i]-min(Covars[,i]))/(rg[2] - rg[1])
    } # end for
    Covars
}

