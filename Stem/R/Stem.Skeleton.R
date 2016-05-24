`Stem.Skeleton` <-
function(...) {
    if (nargs() == 1)
        x <- as.list(...)
    else
        x <- list(...)

# Stem.Skeleton components, phi and K must be provided by the user, p is set to 1 by default
   comp <- c("phi", "p", "K")
   if(!any(names(x)== "p")) x$p = 1
	
#Verify if all the required components are given by user 
	compInd <- match(comp, names(x))
   if (any(is.na(compInd)))
        stop(paste("Component(s)", paste(comp[is.na(compInd)], collapse=", "),
                   "is (are) missing"))
                   
#Controls over the dimensions of the components 
  if (sum(!unlist((lapply(x$phi,is.numeric))))>0)  stop("Component of phi must be numeric and positive")
  if(!(nrow(x$phi$G) == x$p && ncol(x$phi$G) == x$p)) stop("The dimension of matrix G must be p*p")
  if(!(nrow(x$phi$Sigmaeta) == x$p && ncol(x$phi$Sigmaeta) == x$p)) stop("The dimension of matrix Sigmaeta must be p*p")
  if(!(nrow(x$phi$m0)== x$p && ncol(x$phi$m0)== 1)) stop("The dimension of matrix mu0 must be p*1")
  if(!(nrow(x$phi$C0) == x$p && ncol(x$phi$C0) == x$p)) stop("The dimension of matrix C0 must be p*p")
  if(!(length(x$phi$theta) == 1 && (x$phi$theta) > 0)) stop("Theta must be scalar and positive")
  if(!(length(x$phi$sigma2eps) == 1 && (x$phi$sigma2eps) > 0)) stop("Sigma2eps must be scalar and positive")
  if(!(length(x$phi$sigma2omega) == 1 && (x$phi$sigma2omega) > 0)) stop("Sigma2omega must be scalar and positive")
  
#Definition of the class Stem.Skeleton
class(x) <- "Stem.Skeleton"
return(list(phi=x$phi,
		p=x$p,
		K=x$K))
}

