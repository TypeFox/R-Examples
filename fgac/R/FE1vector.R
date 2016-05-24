"FE1vector" <-
function(u,x)
{n<-length(u);
m<-length(x);
	secdirac1vector<-matrix(nrow=n,ncol=m);
	FE1vec<-matrix(nrow=1,ncol=m)
	for(j in 1:m)
	for(i in 1:n)
	{secdirac1vector[i,j]<-dirac1(u[i],x[j]);
		FE1vec[j]<-1/n*sum(secdirac1vector[1:n,j])};
		resultFE1vec<-FE1vec
}

