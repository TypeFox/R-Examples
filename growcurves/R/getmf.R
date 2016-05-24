######################################################
## Take formula input for X and Z and process
## into design matrices
######################################################
#' Produce fixed and random effects design matrices from single formula input
#'
#' An internal function to \code{\link{dpgrow}} and  \code{\link{dpgrowmm}}
#'
#' @export
#' @aliases getmf get.mf
#' @param formula A formula of format \code{y ~ x_1 + x_2*x_3 | z_1*z_2 }
#'	where \code{|} separates fixed (to the left of \code{|}) and random effects.
#' @param random.only A boolean scalar used in the case that either fixed or random effects
#'	are entered in \code{formula}, but not both, which case the \code{|} is not entered
#'	(e.g. \code{y ~ x_1 + x_2*x_3 }.  Then, if \code{random.only == TRUE} the variables
#'	on the right-hand side are interpreted to be random effects; otherwise fixed
#'	for use in \code{\link{dpgrow}} and  \code{\link{dpgrowmm}}.
#' @param data Associated data.frame containing names variables in \code{formula}
#' @param na.action Should be left blank for use in \code{\link{dpgrow}} and \code{\link{dpgrowmm}},
#'	where is automatically set to \code{na.fail}.
#' @return res A list object containing \code{list(y = y, x = x, z = z, m = m, mf = mf)}.
#' @seealso \code{\link{dpgrow}} \code{\link{dpgrowmm}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrow}} and \code{\link{dpgrowmm}}
getmf <- function(formula, random.only, data, na.action)
{
	mf <- match.call(expand.dots = FALSE)
      	m <- match(c("formula", "data", "na.action"), names(mf), 0)
      	mf <- mf[c(1, m)]

	f	<- Formula(formula)
      	mf[[1]] <- as.name("model.frame")
	      mf$formula <- f
	      mf <- eval(mf, parent.frame())

     	y <- model.response(mf)
	if(length(f)[2] == 2) ## fixed + random
	{
     		X <- model.matrix(f,data = mf, rhs = 1)
     		Z <- model.matrix(f, data = mf, rhs = 2) ## includes an intercept, unless formula 
     		                                                  ## includes -1
      	}else{
		if(random.only == FALSE)
		{
			X <- model.matrix(f,data = mf, rhs = 1)
			Z <- NULL
		}else{ ## random.only == TRUE
			X <- NULL
			Z <- model.matrix(f,data = mf, rhs = 1)
		}
	} 
	
	if( !is.null(X) )
	{
		names.x		<- colnames(X)
		X 			<- as.matrix(X[,-1]) ## remove intercept since model has global intercept
		colnames(X)		<- names.x[-1]
	}

     	res 		<- list(y = y, X = X, Z = Z, m = m, mf = mf)

	return(res)
} ## end function get.mf