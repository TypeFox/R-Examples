#
#This code returns the clusters found by KUlldorff's test
#

get.knclusters<-function(d, knresults)
{

	ncl<-nrow(knresults)
	clusters<-as.list( rep(NA, ncl) )

	for(i in 1:ncl)
	{
		x<-knresults$x[i]
		y<-knresults$y[i]
		s<-knresults$size[i]

		dst<-(d$x-x)^2+(d$y-y)^2
		clusters[[i]]<-(order(dst))[1:s]
	}

	return(clusters)
}

