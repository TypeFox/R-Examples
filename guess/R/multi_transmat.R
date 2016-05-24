#' multi_transmat: transition matrix of all the items
#'
#' @title Creates a transition matrix for each item.
#' @description Needs an 'interleaved' dataframe (see interleave function). Pre-test item should be followed by corresponding post-item item etc. 
#' Don't knows must be coded as NA. Function handles items without don't know responses.
#' The function is used internally. It calls transmat.
#' @param pre_test Required. data.frame carrying responses to pre-test questions.
#' @param pst_test Required. data.frame carrying responses to post-test questions.
#' @param subgroup a Boolean vector identifying the subset. Default is NULL.
#' @param force9       Optional. There are cases where DK data doesn't have DK. But we need the entire matrix. By default it is FALSE.
#' @return matrix with rows = total number of items + 1 (last row contains aggregate distribution across items)
#' number of columns = 4 when no don't know, and 9 when there is a don't know option
#' @export
#' @examples
#' pre_test <- data.frame(pre_item1=c(1,0,0,1,0), pre_item2=c(1,NA,0,1,0)) 
#' pst_test <- data.frame(pst_item1=pre_test[,1] + c(0,1,1,0,0), 
#'						  pst_item2 = pre_test[,2] + c(0,1,0,0,1))
#' multi_transmat(pre_test, pst_test)

multi_transmat <- function (pre_test = NULL, pst_test=NULL, subgroup=NULL, force9=FALSE) 
{

	# Checks
	if (!is.data.frame(pre_test)) stop("Specify pre_test data.frame.") # pre_test data frame is missing
	if (!is.data.frame(pst_test)) stop("Specify pst_test data.frame.") # post_test data frame is missing
	if(length(pre_test)!=length(pst_test)) stop("Lengths of pre_test and pst_test must be the same.") # If different no. of items
	
	# Subset
	if (!is.null(subgroup))
	{
		pre_test <- subset(pre_test, subgroup)
		pst_test <- subset(pst_test, subgroup)
	}
	
	# No. of items
	n_items <- length(pre_test)

	# Initialize results
	res <- list()
	
	# Get transition matrix for each item pair
	for (i in 1:n_items)
	{
		# cat("\n Item", i, "\n")
		res[[i]] <- transmat(pre_test[,i], pst_test[,i], force9=force9)
	}

	# Prepping results
	row_names <- paste0("item", 1:n_items) 
	col_names <- names(res[[1]])

	res       <- matrix(unlist(res), nrow=n_items, byrow=T, dimnames=list(row_names, col_names))
	res       <- rbind(res, colSums(res, na.rm=T))
	rownames(res)[nrow(res)] <- "agg"

	#cat("\n Aggregate \n")
	#prmatrix(res)
	#cat("\n")

	return(invisible(res))
}