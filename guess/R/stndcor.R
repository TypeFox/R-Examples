#' Standard Guessing Correction for Learning
#'
#' Estimate of learning adjusted with standard correction for guessing. Correction is based on number of options per question.
#' The function takes separate pre-test and post-test dataframes. Why do we need dataframes? To accomodate multiple items.
#' The items can carry NA (missing). Items must be in the same order in each dataframe. Assumes that respondents are posed same questions twice. 
#' The function also takes a lucky vector -- the chance of getting a correct answer if guessing randomly. Each entry is 1/(no. of options).
#' The function also optionally takes a vector carrying names of the items. By default, the vector carrying adjusted learning estimates takes same item names as the pre_test items.  
#' However you can assign a vector of names separately via item_names.
#' @param pre_test Required. data.frame carrying responses to pre-test questions.
#' @param pst_test Required. data.frame carrying responses to post-test questions.
#' @param lucky    Required. A vector. Each entry is 1/(no. of options)
#' @param item_names Optional. A vector carrying item names.
#' @return   a list of three vectors, carrying pre-treatment corrected scores, post-treatment scores, and adjusted estimates of learning
#' @export
#' @examples
#' pre_test <- data.frame(item1=c(1,0,0,1,0), item2=c(1,NA,0,1,0)); 
#' pst_test <- pre_test + cbind(c(0,1,1,0,0), c(0,1,0,0,1))
#' lucky <- rep(.25, 2); stndcor(pre_test, pst_test, lucky)
 
stndcor <- function(pre_test=NULL, pst_test=NULL, lucky=NULL, item_names=NULL)
{	
	if (!is.data.frame(pre_test)) stop("Specify pre_test data.frame.") # pre_test data frame is missing
	if (!is.data.frame(pst_test)) stop("Specify pst_test data.frame.") # post_test data frame is missing
	if (is.null(lucky))    stop("Specify lucky vector.")           # lucky vector is missing
	
	if (length(unique(length(pre_test), length(pst_test), length(lucky)))!=1) stop("Length of input varies. Length of pre_test, pst_test, and lucky must be the same.")

	n_items <- length(pre_test) # total number of items
	pre_test_cor  <- pst_test_cor <- stnd_cor <- rep(NA, length=n_items) # initialize some results vectors
	
	# Names of the return vector
	if (is.null(item_names)) {
		names(pre_test_cor) <- names(pst_test_cor) <- names(stnd_cor) <- names(pre_test)
	} else {
		names(pre_test_cor) <- names(pst_test_cor) <- names(stnd_cor) <- item_names
	}

	for (i in 1:n_items)
	{
		pre_test_cor[i] <- sum(pre_test[,i], na.rm = TRUE) - length(which(pre_test[,i] == 0))/(1/lucky[i] - 1)
		pst_test_cor[i] <- sum(pst_test[,i], na.rm = TRUE) - length(which(pst_test[,i] == 0))/(1/lucky[i] - 1)
		stnd_cor[i]     <- pst_test_cor[i] - pre_test_cor[i]
	}

	pre   <- pre_test_cor/nrow(pre_test)
	pst   <- pst_test_cor/nrow(pst_test)
	learn <- stnd_cor/nrow(pre_test)

	return(list(pre=pre, pst=pst, learn=learn))
}