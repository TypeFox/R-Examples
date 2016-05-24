#' Area Under the Curve of the Reciever Operating Curve
#' 
#' \code{auc} estimates the AUC of the ROC using a Mann-Whitney U statistic.
#' \cr \cr \bold{Note:} this method will exclude any missing data.
#' 
#' 
#' @param obs a vector of observed values which must be 0 for absences and 1
#' for occurrences
#' @param pred a vector of the same length as \code{obs} representing the
#' predicted values. Values must be between 0 & 1 representing a likelihood.
#' @return Returns a single value represting the AUC value.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{Kappa}}, \code{\link{omission}},
#' \code{\link{sensitivity}}, \code{\link{specificity}},
#' \code{\link{prop.correct}}, \code{\link{confusion.matrix}},
#' \code{\link{accuracy}}
#' @examples
#' 
#' 
#' #create some data
#' obs = c(sample(c(0,1),20,replace=TRUE),NA)
#' pred = runif(length(obs),0,1)
#' 
#' #calculate AUC from the random data
#' auc(obs,pred)
#' 
#' #calculate an example 'perfect' AUC
#' obs = obs[order(obs)]
#' pred = pred[order(pred)]
#' auc(obs,pred)
#' 
#' 
#' @export
auc <- function(obs,pred){
	#input checks
	if (length(obs)!=length(pred)) stop('this requires the same number of observed & predicted values')
	
	
	#deal with NAs
	if (length(which(is.na(c(obs,pred))))>0) {
		na = union(which(is.na(obs)),which(is.na(pred)))
		warning(length(na),' data points removed due to missing data')
		obs = obs[-na]; pred = pred[-na]
	}

	#define the n's and do checks
	n = length(obs); if (length(which(obs %in% c(0,1)))!=n) stop('observed values must be 0 or 1') #ensure observed are values 0 or 1
	n1 = as.double(length(which(obs==1))); n0 = as.double(length(which(obs==0)))
	if (n1==0 || n1==n) return( NaN ) #if all observed 1's or 0's return NaN

	###calc AUC
	pred0 = pred[which(obs==0)]
	pred1 = pred[which(obs==1)]
	ranks = rank(pred,ties.method='average')#define ranks
	ranks0 = ranks[which(obs==0)]
	ranks1 = ranks[which(obs==1)]
	U = n0*n1 + (n0*(n0+1))/2 - sum(ranks0) #calc U stat
	AUC = U/(n0*n1) #estimate AUC
	if (AUC<.5) AUC = 1-AUC
	
	#return the auc value
	return(AUC)
}

