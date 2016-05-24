windowsize<-function(x,y,winmin,winmax,winstep)
{
	## the value 100 is randomly selected,
    ## usually we will be able to identify the best
    ## window size before the 100nd test. 
	error=rep(0,100)
	win = rep(0,100)
	
	for (i in 1:100)
	{
	dx=winmin+(i-1)*winstep
	
	if (dx > winmax)
		stop
	else
		{win[i]=dx
		gout=proxyJK(x,y,dx)
		error[i]=sum((gout$delw)^2)/length(gout$delw)	
		}
	}
	
	return (list(win=win,error=error))
}

