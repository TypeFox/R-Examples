"SOB2" <-
function(u,v,x,y)
{n<-length(u);
	secdiracS2<-matrix(nrow=n, ncol=1);
	for(i in 1:n)
		secdiracS2[i]<-diracS2(u[i], v[i], x, y);
		SOB2<-1/n*sum(secdiracS2[1:n])
		resultSOB2<-SOB2
}

