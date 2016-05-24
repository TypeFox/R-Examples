#' Measures of Model Accuracy
#' 
#' \code{accuracy} estimates six measures of accuracy for presence-absence or
#' presence-psuedoabsence data. These include AUC, ommission rates,
#' sensitivity, specificity, proportion correctly identified and Kappa. \cr \cr
#' \bold{Note:} this method will exclude any missing data.
#' 
#' 
#' @param obs a vector of observed values which must be 0 for absences and 1
#' for occurrences
#' @param pred a vector of the same length as \code{obs} representing the
#' predicted values.  Values must be between 0 & 1 prepresenting a likelihood.
#' @param threshold this can be: \cr a) a single value representing a single
#' threshold between 0 & 1; \cr b) a vector of threshold values between 0 & 1;
#' OR \cr c) an integer value representing the number of equal interval
#' threshold values between 0 & 1
#' @return a data.frame with seven columns: \item{threshold}{the threshold
#' values representing each row of data} \item{AUC}{the AUC given the defined
#' threshold value} \item{ommission.rate}{the ommission rate as a proportion of
#' true occurrences misidentified given the defined threshold value}
#' \item{sensitivity}{the sensitivity given the defined threshold value}
#' \item{specificity}{the specificity given the defined threshold value}
#' \item{prop.correct}{the proportion of the presence and absence records
#' correctly identified given the defined threshold value} \item{Kappa}{the
#' Kappa statistic of the model given the defined threshold value}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{auc}}, \code{\link{Kappa}}, \code{\link{omission}},
#' \code{\link{sensitivity}}, \code{\link{specificity}},
#' \code{\link{prop.correct}}, \code{\link{confusion.matrix}}
#' @examples
#' 
#' 
#' #create some data
#' obs = c(sample(c(0,1),20,replace=TRUE),NA); obs = obs[order(obs)]
#' pred = runif(length(obs),0,1); pred = pred[order(pred)]
#' 
#' #calculate accuracy of the model with a single threshold value
#' accuracy(obs,pred,threshold=0.5)
#' 
#' #calculate accuracy given several defined thresholds
#' accuracy(obs,pred,threshold=c(0.33,0.5,0.66))
#' 
#' #calculate accuracy given a number of equal interval thresholds
#' accuracy(obs,pred,threshold=20)
#' 
#' 
#' @export 
accuracy <- function(obs,pred,threshold=0.5){
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

	# check / setup the threshold values
	if (length(threshold)==1 & threshold[1]<=1 & threshold[1]>=0) { 
		thresholds = threshold
	} else if (length(threshold)==1 & threshold[1]>1) {
		thresholds = seq(0,1,length=threshold)
	} else if (length(threshold)>1 & max(threshold)<=1 & min(threshold)>=0) {
		thresholds = threshold
	} else { stop('inappropriate threshold values used as input. See help file.') }
	
	#cycle through each of the helpfiles
	out = data.frame(threshold=as.double(thresholds),AUC=NA,omission.rate=NA,sensitivity=NA,
		specificity=NA,prop.correct=NA,Kappa=NA)
	for (ii in 1:length(thresholds)) {
		threshold = thresholds[ii]
		#convert preditions to binary based on threshold
		bin.pred = pred;
		if (threshold==0) {
			bin.pred[which(bin.pred>threshold)] = 1; bin.pred[which(bin.pred<=threshold)] = 0
		} else {
			bin.pred[which(bin.pred>=threshold)] = 1; bin.pred[which(bin.pred<threshold)] = 0
		}
		#create the confusion matrix ...just like using confusion.matrix command
		mat = table(bin.pred=factor(bin.pred,levels=c(0,1)),obs=factor(obs,levels=c(0,1)))
		attr(mat,'class') = 'confusion.matrix'
		#sppend data to output dataframe
		out$AUC[ii] = auc(obs,bin.pred)
		out$omission.rate[ii] = omission(mat)
		out$sensitivity[ii] = sensitivity(mat)
		out$specificity[ii] = specificity(mat)
		out$prop.correct[ii] = prop.correct(mat)
		out$Kappa[ii] = Kappa(mat)
	}
	#return the output data
	return(out)
}

