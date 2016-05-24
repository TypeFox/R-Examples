#' Measures of Accuracy
#' 
#' Estimates different measures of accurracy given a confusion matrix.
#' 
#' 
#' @param mat a confusion matrix of class 'confusion.matrix' from
#' \code{confusion.matrix}
#' @return returns single values representing the: \item{ommission}{the
#' ommission rate as a proportion of true occurrences misidentified given the
#' defined threshold value} \item{sensitivity}{the sensitivity given the
#' defined threshold value} \item{specificity}{the specificity given the
#' defined threshold value} \item{prop.correct}{the proportion of the presence
#' and absence records correctly identified given the defined threshold value}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{auc}}, \code{\link{Kappa}},
#' \code{\link{confusion.matrix}}, \code{\link{accuracy}}
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
#' #calculate the accuracy measures
#' omission(mat)
#' sensitivity(mat)
#' specificity(mat)
#' prop.correct(mat)
#' 
#' 
#' @export 
omission <- function(mat){
	#input checks
	if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
	#return the value
	return(mat[1,2]/sum(mat[,2]))
}

#' @rdname omission
#' @export
sensitivity = function(mat) { 
	#input checks
	if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
	#return the value
	return(mat[2,2]/sum(mat[,2])) 
}

#' @rdname omission
#' @export
specificity = function(mat) { 
	#input checks
	if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
	#return the value
	return(mat[1,1]/sum(mat[,1])) 
}

#' @rdname omission
#' @export
prop.correct = function(mat) { 
	#input checks
	if (attr(mat,'class')!='confusion.matrix') stop('mat must be of class confusion.matrix')
	#return the value
	return(sum(diag(mat))/sum(mat)) 
}
