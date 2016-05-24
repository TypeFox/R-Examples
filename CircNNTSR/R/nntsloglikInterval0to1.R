nntsloglikInterval0to1 <-
function (data, cutpoints, cpars = 1/sqrt(2*pi), M = 0 ) 
{
#    temp<-cutpoints[(length(data)+1)]%%1
#    if (temp<>0)
#	cutpoints[(length(data)+1)]<-cutpoints[(length(data)+1)]-0.0000000001
    y <- 0
    for (k in 1:length(data)) 
    {
	if (cutpoints[k]==cutpoints[(k+1)])
		return("Subinterval of length 0")
	a<-(cutpoints[(k+1)]-cutpoints[k])%%1
	if (a==0)
		y<-y
	else {
		s1<-cutpoints[(k+1)]%%1
		s2<-cutpoints[k]%%1
		if (s1<s2)
        		y <- y + data[k] * log(1+nntsDistributioninterval0to1(cutpoints[(k + 1)], cpars, M) 
                                        -nntsDistributioninterval0to1(cutpoints[k], cpars, M)) 
		else 
			y <- y + data[k] * log(nntsDistributioninterval0to1(cutpoints[(k + 1)], cpars, M) 
                                        -nntsDistributioninterval0to1(cutpoints[k], cpars, M)) 
	}
    }
    res <- Re(y)
    return(res)
}

