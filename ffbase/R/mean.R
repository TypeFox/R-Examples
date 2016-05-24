#' Mean of ff vector
#' 
#' @method mean ff
#' @export
#' @export mean.ff
#' @example ../examples/mean.R
#' @param x a ff vector
#' @param trim percentage of robustness, between 0 and 1
#' @param ... other arguments passed to \code{mean}
#' @param range a \code{ri} or an \code{integer} vector of \code{length==2} giving a range restriction for chunked processing
#' @importFrom stats weighted.mean
#' @return mean value
mean.ff <- function(x, trim=0, ..., range=NULL){
    r <- checkRange(range, x)
    
    if (trim > 0){
	   trim <- min(trim, 0.5)
      #TODO if range then make a selection
	   x <- ffsort(x)  #sort x
	   n <- length(x)     #calculate length
	   r[1] <- floor(n * trim)
	   r[2] <- n - r[1]
	  }
	
   
   res <- sapply( chunk(x, from=min(r), to=max(r))
                , function(i){
                     Log$chunk(i)
                     c( mean=mean(x[i], ...)
                      , w = sum(i)/max(r)
                      )
                  }
				    )
    weighted.mean(res['mean',], res['w',])
}
