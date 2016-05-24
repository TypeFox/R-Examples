dist.BC<-function(x)
{
	if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
	}
	nr<-nrow(x)
	wynik<-matrix(nrow=nr,ncol=nr,dimnames=names(x))
	for (i in 1:nr)
	for (j in 1:nr)
	{
		suma1=0
		suma2=0
		if (i!=j)
		for (k in 1:ncol(x))
		{
			suma1=suma1+abs(x[i,k]-x[j,k])
			suma2=suma2+(x[i,k]+x[j,k])
			if(suma2!=0)
				wynik[i,j]=suma1/suma2
		}
			
	}
	as.dist(wynik)
}

