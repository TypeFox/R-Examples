#Hinda Haned
# December 2009
PV <-function(mat,prior)
{
	#mat is the matrix of the Pr(\hat x=i|x=i)
	nc<-ncol(mat)#maximum number of contributors
	nr<-nrow(mat)#estimates range
	#standards verif
	if(!is.vector(prior) | length(prior)<nr|!is.numeric(prior) | sum(prior)!=1)
	{
		stop("argument prior must be a vector of numerics summing to 1")
	}
	#in case of over estimates of the ML, or underestimates

	deno<-NULL
	num<-NULL
	if(nr<=nc)
	{
		for(i in 1:nr)
		{
		num<-c(num,mat[i,i]*prior[i])
		deno<-c(deno,sum(mat[,i]*prior))

		}
		return(num/deno)
	}

	else
	{
		for(i in 1:nc)
		{
		num<-c(num,mat[i,i]*prior[i])
		deno<-c(deno,sum(mat[,i]*prior))

		}
		return(num/deno)

	}
}
