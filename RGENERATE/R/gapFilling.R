
# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL

#' gapFilling
#' 
#' 
#' 
#' @rdname gapFilling
#' @export

gapFilling <- function (x=NULL,...)  {
	
	
	return(UseMethod("gapFilling",x))
	
}




NULL
#' gabFilling method 
#' 
#' It fills in a gab of a data frame by using  \code{\link{generate}} method
#' 
#' @param x object with gaps to fill 
#' @param objectForGeneration object used for \code{\link{generate}} method
#' @param max.filling integer values: gap are filled if the previous \code{max.filling} values are not \code{NA} or  \code{nofill.code} 
#' @param nofill.code Alternative value to \code{NA} which indicates the gaps which are not filled
#' @param ... further argument for \code{\link{generate}} method
#' 
#' 
#' @examples 
#' 
#' set.seed(122)
#' NSTEP <- 1000
#' x <- rnorm(NSTEP)
#' y <- x+rnorm(NSTEP)
#' z <- c(rnorm(1),y[-1]+rnorm(NSTEP-1))
#' df <- data.frame(x=x,y=y,z=z)
#' var <- VAR(df,type="none")
#' 
#' dfobs <- df
#' dfobs[20:30,2] <- NA
#' n <- nrow(df)
#' gp <- gapFilling(x=dfobs,objectForGeneration=var,max.filling=2)
#' 
#' 
#' 
#' 
#' 
#' @rdname gapFilling
#' @method gapFilling default
#' @S3method gapFilling default
#' @aliases gapFilling 
#' @export



gapFilling.default <- function (x,objectForGeneration=NULL,...)  {
	
	
	out <- generate(x=objectForGeneration,gap.filling=x,...)
	
	
	
	return(out)
	
}


NULL
#' @rdname gapFilling
#' @method gapFilling data.frame
#' @S3method gapFilling data.frame
#' @aliases gapFilling 
#' @export


gapFilling.data.frame <- function (x,objectForGeneration=NULL,max.filling=2,nofill.code=-9999,...)  {
	
	isna <- is.na(x)
	
	for (c in 1:ncol(x)) { 
		
#		
		if (nrow(x)>max.filling) {
			

			
			for (r in (max.filling+1):nrow(x)) {
			
				cond <- isna[r,c] & isna[r-max.filling,c]
			
				if (cond) x[r,c] <- nofill.code
				
			}
		}
		
	
	
	}


	out <- generate(x=objectForGeneration,gap.filling=x,n=nrow(x),...)
	
	out[out==nofill.code] <- NA 
	
	return(out)
	
}








