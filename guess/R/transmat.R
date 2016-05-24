#' transmat: Cross-wave transition matrix
#'
#' @description Prints Cross-wave transition matrix and returns the vector behind the matrix.  
#' Missing values are treated as ignorance. Don't know responses need to be coded as 'd'.
#' @param pre_test_var Required. A vector carrying pre-test scores of a particular item. Only 
#' @param pst_test_var Required. A vector carrying post-test scores of a particular item
#' @param subgroup     Optional. A Boolean vector indicating rows of the relevant subset.  
#' @param force9       Optional. There are cases where DK data doesn't have DK. But we need the entire matrix. By default it is FALSE.
#' @return a numeric vector. 
#' Assume 1 denotes correct answer, 0 and NA incorrect, and d 'don't know.'
#' When there is no don't know option and no missing, the entries are: x00, x10, x01, x11
#' When there is a don't know option, the entries of the vector are: x00, x10, xd0, x01, x11, xd1, xd0, x1d, xdd
#' @export
#' @examples
#' pre_test_var <- c(1,0,0,1,0,1,0)
#' pst_test_var <- c(1,0,1,1,0,1,1)
#' transmat(pre_test_var, pst_test_var)
#' 
#' # With NAs
#' pre_test_var <- c(1,0,0,1,"d","d",0,1,NA)
#' pst_test_var <- c(1,NA,1,"d",1,0,1,1,"d") 
#' transmat(pre_test_var, pst_test_var)

transmat <- function(pre_test_var, pst_test_var, subgroup=NULL, force9=FALSE) 
{	

	if (!is.null(subgroup))
	{
		pre_test_var <-subset(pre_test_var, subgroup)
		pst_test_var <-subset(pst_test_var, subgroup)
	}

	# No NAs
	pre_test_nona <- nona(pre_test_var)
	pst_test_nona <- nona(pst_test_var)
	
	# Check if the vector has only one of the 4 values
	if (!all(unique(c(pre_test_nona, pst_test_nona)) %in% c(NA, "1", "0", "d"))) stop("The input vectors can only contain: 0, 1, NA, d")

	pre_pst <- paste0(pre_test_nona, pst_test_nona)

	# Get the numbers
	x00     <- sum(pre_pst=="00")
	x10     <- sum(pre_pst=="10")  
	x01     <- sum(pre_pst=="01")  
	x11     <- sum(pre_pst=="11")  
	xd0     <- sum(pre_pst=="d0") 
	xd1     <- sum(pre_pst=="d1") 
	x1d     <- sum(pre_pst=="1d") 
	xdd     <- sum(pre_pst=="dd")
	x0d     <- sum(pre_pst=="0d")
	
	res <- c(x00, x01, x10, x11)
	names(res) <- c("x00", "x01", "x10", "x11")

	if ((xd0 + xd1 + x1d + xdd != 0) | force9) {

		res <- c(x00, x01, x0d, x10, x11, x1d, xd0, xd1, xdd)
		names(res) <- c("x00", "x01", "x0d", "x10", "x11", "x1d", "xd0", "xd1", "xdd")
	} 

    return(invisible(res))

}