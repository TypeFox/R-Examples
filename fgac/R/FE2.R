"FE2" <-
function(u,v,x,y)
{n<-length(u);
	secdirac2<-matrix(nrow=n, ncol=1);
	for(i in 1:n)
		secdirac2[i]<-dirac2(u[i], v[i], x, y);
		FE2<-1/n*sum(secdirac2[1:n])
		resultFE2<-FE2
}

