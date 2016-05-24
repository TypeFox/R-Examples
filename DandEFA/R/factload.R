factload <- function(data,cormeth="spearman",data.cor,method="pc",nfac=1,rotation="varimax")
{
	if(!missing(data))
	{

	if(cormeth=="polycor")
	{
        poly<-matrix(rep(1,ncol(data)),nrow=ncol(data),ncol=ncol(data))
	  
	  for (i in 1:(ncol(data)-1))
	  {
	    for(j in (i+1):ncol(data))
	    {
	      
	      poly[i,j]<-polychor(data[,i],data[,j])
	      poly[j,i]<-poly[i,j]
	    }
	  }
      data.cor <- poly
 	}
	else data.cor <- cor(data,method=cormeth)
	temp <- colnames(data)
	}

	if(missing(data.cor))
	{
	stop("please identify either a data or a correlation matrix")
	}

	
	if(method=="pc")
	{
		load <- eigen(data.cor)$vectors[,1:nfac]
		val <- eigen(data.cor)$values[1:nfac]
		load <- load %*% diag(sqrt(val))
	}
	else if(method=="mle")
	{
		load <- factanal(covmat=data.cor,factors=nfac, rotation="none")$loadings
      }
	else if(method=="prax")
	{
	  p <- nrow(data.cor)
	  cm       <- data.cor^2
	  diag(cm) <- 0
	  h       <- apply(cm,2,max)
        e  <- 1-h
        for( i in 1:1000)
        {
            eold <- e 
            r <- data.cor-diag(e)
  	      load  <- eigen(r)$vectors
            val <- eigen(r)$values
   	      load <- load%*%diag(sqrt(val))[,1:nfac]
            e  <- diag(data.cor-load%*%t(load))   
            if(max(abs(e-eold)) < 0.001) break    
	  }     
      }
      if(rotation=="varimax")
	{
		load <- varimax(load)$loadings
	}	
      else if(rotation=="promax")
	{
		load <- promax(load)$loadings
	}
	else if(rotation=="none") class(load) <- "loadings"
	else stop("please enter none, varimax or promax")
	if(exists("temp"))
		rownames(load) <- temp
	load
}


	