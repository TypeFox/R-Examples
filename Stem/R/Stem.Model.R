`Stem.Model` <-
function(...) {
	
     if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

    #Stem.Model components : skeleton and data
    skeleton <- Stem.Skeleton(x) #(phi=phi, p=p, K=K)
    
    data <- Stem.Data(x) #(z=z, coordinates=coordinates, covariates=covariates)
    
    if(length(skeleton$phi$beta) != ncol(data$covariates)) stop("The length of Beta must be equal to the number of columns of covariates")
    if(!(nrow(skeleton$K) == data$d && ncol(skeleton$K) == skeleton$p)) stop("The dimension of matrix K must be d*p")
    
	  
    x=list(skeleton=skeleton,data=data)
    class(x) <- "Stem.Model"
    return(x)	
}

