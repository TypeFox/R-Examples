#' Confusion Matrix
#' 
#' \code{confusion.matrix} calculates a confusion matrix. \cr \cr \bold{Note:}
#' this method will exclude any missing data
#' 
#' 
#' @param obs a vector of observed values which must be 0 for absences and 1
#' for occurrences
#' @param pred a vector of the same length as \code{obs} representing the
#' predicted values. Values must be between 0 & 1 prepresenting a likelihood.
#' @param threshold a single threshold value between 0 & 1
#' @return Returns a confusion matrix (table) of class 'confusion.matrix'
#' representing counts of true & false presences and absences.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{auc}}, \code{\link{Kappa}}, \code{\link{omission}},
#' \code{\link{sensitivity}}, \code{\link{specificity}},
#' \code{\link{prop.correct}}, \code{\link{accuracy}}
#' @examples
#' 
#' 
#' #create some data
#' obs = c(sample(c(0,1),20,replace=TRUE),NA); obs = obs[order(obs)]
#' pred = runif(length(obs),0,1); pred = pred[order(pred)]
#' 
#' #calculate the confusion matrix
#' confusion.matrix(obs,pred,threshold=0.5)
#' 
#' 
#' @export 
confusion.matrix <- function(obs,pred,threshold=0.5){
	#input checks
	if (length(obs)!=length(pred)) stop('this requires the same number of observed & predicted values')
	if (!(length(threshold)==1 & threshold[1]<=1 & threshold[1]>=0)) stop('inappropriate threshold value... must be a single value between 0 & 1.')
	n = length(obs); if (length(which(obs %in% c(0,1,NA)))!=n) stop('observed values must be 0 or 1') #ensure observed are values 0 or 1
	
	#deal with NAs
	if (length(which(is.na(c(obs,pred))))>0) {
		na = union(which(is.na(obs)),which(is.na(pred)))
		warning(length(na),' data points removed due to missing data')
		obs = obs[-na]; pred = pred[-na]
	}
	#apply the threshold to the prediction
	if (threshold==0) {
		pred[which(pred>threshold)] = 1; pred[which(pred<=threshold)] = 0
	} else {
		pred[which(pred>=threshold)] = 1; pred[which(pred<threshold)] = 0
	}
	#return the confusion matrix
	mat = table(pred=factor(pred,levels=c(0,1)),obs=factor(obs,levels=c(0,1)))
	attr(mat,'class') = 'confusion.matrix'
	return(mat)
}

