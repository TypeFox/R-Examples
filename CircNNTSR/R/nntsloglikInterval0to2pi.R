nntsloglikInterval0to2pi <-
function (data, cutpoints, cpars = 1/sqrt(2*pi), M = 0) 
{
    #data <- data/sum(data)
    y <- 0
    for (k in 1:length(data)) 
    {
	if (cutpoints[k]==cutpoints[(k+1)])
		return("Subinterval of length 0")
	a<-(cutpoints[(k+1)]-cutpoints[k])%%(2*pi)
	if (a==0)
		y<-y
	else {
		s1<-cutpoints[(k+1)]%%(2*pi)
		s2<-cutpoints[k]%%(2*pi)
		if (s1 < s2)
        		y <- y + data[k] * log(1 + nntsDistributioninterval0to2pi(cutpoints[(k+1)],cpars, M)
                                        -nntsDistributioninterval0to2pi(cutpoints[k], cpars, M))
		else
        		y <- y + data[k] * log(nntsDistributioninterval0to2pi(cutpoints[(k + 1)], cpars, M)
					-nntsDistributioninterval0to2pi(cutpoints[k], cpars, M))
	}
    }
    res <- Re(y)
    return(res)
}

