`Stem.Data` <-
function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

# Stem.Data components
  comp <- c("z", "coordinates", "covariates")

#Verify if all the required components are given by user	
  compInd <- match(comp, names(x))
    if (any(is.na(compInd)))
        stop(paste("Component(s)", paste(comp[is.na(compInd)], collapse=", "),
                   "is (are) missing"))
                   
#Controls over the dimensions of the components                    
   if (sum(!(apply(x$z,1,is.numeric)))>0)  stop("Component z must be numeric")
   if (sum(!(apply(x$coordinates,1,is.numeric)))>0)  stop("Component coordinates must be numeric")
   if (sum(!(apply(x$covariates,1,is.numeric)))>0)  stop("Component covariates must be numeric")
   r <- ncol(x$covariates)
   n <- nrow(x$z)
   d <- ncol(x$z)
      
   if(nrow(x$covariates) != n*d) stop("The number or row of covariates must be n*d")   
   if (any( c(is.na(x$z), is.na(x$coordinates), is.na(x$covariates))))
        stop("Missing values are not allowed in components z, coordinates, covariates")
   if(!(nrow(x$coordinates) == d && ncol(x$coordinates) == 2)) stop("The dimension of matrix coordinates must be d*2")
  
#Definition of the class Stem.Data  
  x$r=r
  x$n=n
  x$d=d
  class(x) <- "Stem.Data"
   return(list(z=x$z,
		coordinates=x$coordinates,
		covariates=x$covariates,
		r=x$r,
		n=x$n,
		d=x$d))
}

