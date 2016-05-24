# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  March 2010
# Version beta
# Licence GPL v3

setGeneric("normalize", function(x, ...) standardGeneric("normalize"))

setMethod("normalize", signature(x = "TransitionLayer"), def = function(x, method="row")
	{
		tr <- transitionMatrix(x)
		tr <- .normalize(tr, method)
		transitionMatrix(x) <- tr
		return(x)
	}
)

.normalize <- function(x, method)
	{
		
		if(!(method %in% c("row","col","symm"))){stop("invalid method argument")}
		if(method=="symm")
		{
			rs <- rowSums(x)^-.5
			cs <- colSums(x)^-.5

      tr <- x * rs
      tr <- t(tr)
			tr <- tr * cs
    
      tr <- t(tr)

      if(isSymmetric(x)) 
      {
        tr <- forceSymmetric(tr)
        tr <- as(tr, "CsparseMatrix")
      }
		}

    if(method=="row")
		{
			rs <- 1 / rowSums(x)
			rs[rs == Inf] <- 0
			tr <- x * rs
		}

    if(method=="col")
		{
			rs <- 1 / colSums(x)
			rs[rs == Inf] <- 0
			tr <- t(t(x) * rs)
		}

		return(tr)
	}