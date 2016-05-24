pcscorebootstrapdata <- function(dat, bootrep, statistic, bootmethod = c("st", "sm", "mvn", "stiefel", "meboot"), smo)
{
	dat = t(dat)
	n = dim(dat)[1]
	p = dim(dat)[2]
	boot = matrix(, p, bootrep)
	bootdataarray = array(, dim = c(p, n, bootrep))

	data = scale(dat, center = TRUE, scale = FALSE)
	meandata = matrix(rep(as.matrix(colMeans(dat)), n), p, n)

	dummy = svd(data)
	U = dummy$u
	D = diag(dummy$d)
	V = dummy$v
	if(bootmethod == "st")
	{
		for(b in 1:bootrep)
		{
			drawU = U[sample(1:n, n, replace = TRUE), ]
			bootdata = t(drawU%*%D%*%t(V)) + meandata
			bootdataarray[,,b] = bootdata
			boot[,b] = centre(bootdata, type = statistic)
		}
	}
	if(bootmethod == "sm")
	{
		for(b in 1:bootrep)
		{
			drawU = U[sample(1:n, n, replace = TRUE), ] + mvrnorm(n, rep(0, min(p, n)), var(U) * smo)
			bootdata = t(drawU%*%D%*%t(V)) + meandata
			bootdataarray[,,b] = bootdata
			boot[,b] = centre(bootdata, type = statistic)
		}
	}
	if(bootmethod == "mvn")
	{
		for(b in 1:bootrep)
		{
			drawU <- mvrnorm(n, colMeans(U), var(U))
			bootdata = t(drawU%*%D%*%t(V)) + meandata
			bootdataarray[,,b] = bootdata
			boot[,b] = centre(bootdata, type = statistic)
		}
	}
	if(bootmethod == "stiefel")
	{
		for(b in 1:bootrep)
		{
			X <- mvrnorm(n, colMeans(U), var(U))
			tmp <- eigen(t(X)%*%X)
			drawU = (X%*%(tmp$vec%*%sqrt(diag(1/tmp$val))%*%t(tmp$vec)))[sample(1:n,n,replace=TRUE),]
			bootdata = t(drawU%*%D%*%t(V)) + meandata
			bootdataarray[,,b] = bootdata
			boot[,b] = centre(bootdata, type = statistic)
		}
	}
	if(bootmethod == "meboot")
	{
		out <- meboot::meboot(ts(as.numeric(U),start=1,end=n*(min(n,p))), reps = bootrep)
  		for(b in 1:bootrep)
  		{
  			drawU = matrix(out$ensemble[,b], n, min(n,p))
        	bootdata = t(drawU%*%D%*%t(V)) + meandata
        	bootdataarray[,,b] = bootdata
			boot[,b] = centre(bootdata, type = statistic)
    	}
	}	
	return(list(bootdata = bootdataarray, meanfunction = boot))
}
