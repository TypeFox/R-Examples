#' fit_nodk
#'
#' @title Goodness of fit statistics for data with don't know
#' 
#' @description For data with Don't Know, chi-square goodness of fit between true and model based multivariate distribution
#'
#' @param pre_test data.frame carrying pre_test items
#' @param pst_test data.frame carrying pst_test items 
#' @param g 		estimates of \eqn{\gamma} produced from \code{\link{guesstimate}}
#' @param est.param estimated parameters produced from \code{\link{guesstimate}}
#' @param force9       Optional. There are cases where DK data doesn't have DK. But we need the entire matrix. By default it is FALSE.
#' @return matrix with two rows: top row carrying chi-square value, and bottom row probability of observing that value
#' @export
#' @examples
#' \dontrun{fit_dk(pre_test, pst_test, g, est.param)}

fit_dk <- function(pre_test, pst_test, g, est.param, force9=FALSE) 
{

	data    <- multi_transmat(pre_test, pst_test, force9=force9)
	data    <- data[(1:nrow(data)-1),] # remove the agg.
	expec	<- matrix(ncol=nrow(data), nrow=9)
	fit		<- matrix(ncol=nrow(data), nrow=2)
	colnames(fit) <- rownames(data)
	rownames(fit) <- c("chi-square", "p-value")

	for(i in 1:nrow(data)){
		gi			<- g[[i]]
		expec[1, i]	<- (1 - gi)*(1-gi)*est.param[1,i]*sum(data[i,])
		expec[2, i]	<- ((1 - gi)*gi*est.param[1,i] + (1 - gi)*est.param[2,i])*sum(data[i,])
		expec[3, i]	<- ((1 - gi)*est.param[3,i])*sum(data[i,])
		expec[4, i]	<- ((1 - gi)*gi*est.param[1,i])*sum(data[i,])
		expec[5, i]	<- (gi*gi*est.param[1,i]+gi*est.param[2,i]+est.param[4,i])*sum(data[i,])
		expec[6, i]	<- (gi*est.param[3,i])*sum(data[i,])
		expec[7, i]	<- ((1 - gi)*est.param[5,i])*sum(data[i,])
		expec[8, i]	<- (gi*est.param[5,i] + est.param[6,i])*sum(data[i,])
		expec[9, i]	<- est.param[7,i]*sum(data[i,])
		test 		<- suppressWarnings(chisq.test(expec[,i], p=data[i,]/sum(data[i,])))
		fit[1:2,i]	<- round(unlist(test[c(1,3)]),3)
	}


	fit
}
