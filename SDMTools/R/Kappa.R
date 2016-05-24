#' Kappa Statistic
#' 
#' \code{Kappa} estimates the Kappa statistic for model accuracy.
#' 
#' 
#' @param mat a confusion matrix of class 'confusion.matrix' from
#' \code{confusion.matrix}
#' @return Returns a single value represting the Kappa statistic.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{auc}}, \code{\link{omission}},
#' \code{\link{sensitivity}}, \code{\link{specificity}},
#' \code{\link{prop.correct}}, \code{\link{confusion.matrix}},
#' \code{\link{accuracy}}
#' @examples
#' 
#' 
#' #create some data
#' obs = c(sample(c(0,1),20,replace=TRUE),NA); obs = obs[order(obs)]
#' pred = runif(length(obs),0,1); pred = pred[order(pred)]
#' 
#' #calculate the confusion matrix
#' mat = confusion.matrix(obs,pred,threshold=0.5)
#' 
#' #calculate the Kappa statistic
#' Kappa(mat)
#' 
#' 
#' @export 
Kappa <- function(mat){
	#input checks
	if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
	#calculate Kappa
	n<-sum(mat)
	colsums<-as.double(apply(mat,2,sum))
	rowsums<-as.double(apply(mat,1,sum))
	t1 = sum(diag(mat))/n; 	t2 = sum(rowsums*colsums)/(n^2)
	#return the value
	return((t1-t2)/(1-t2))
}
