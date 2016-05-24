#' Estimates propensity scores using the recursive partitioning in a conditional inference framework.
#' 
#' This function will estimate propensity scores using the conditional inference
#' framework as outlined in the \code{party} package. Specifically, a separate
#' tree will be estimated for each level 2 (or cluster). A key advantage of this
#' framework over other methods for estimating propensity scores is that this
#' method will work on data sets containing missing values.
#' 
#' @param vars a data frame containing the covariates to use for estimating the propensity scores.
#' @param formula the model for estimating the propensity scores. For example,
#'        treat ~ .
#' @param level2 the name of the column in \code{vars} specifying the level 2 (or cluster).
#' @param ... currently unused.
#' @references Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased 
#'             Recursive Partitioning: A Conditional Inference Framework. Journal
#'             of Computational and Graphical Statistics, 15(3), 651--674.
#' @return a list of BinaryTree-class classes for each level 2
#' @seealso \code{\link{getStrata}}
#' @seealso \code{\link{tree.plot}}
#' @export
mlpsa.ctree <- function(vars, formula, level2, ...) {
	partyPlyr <- function(x) {
		excludeVars = names(x) %in% c(level2)
		x = x[,!excludeVars]
		tmp.party = ctree(formula, data=x)
		return(tmp.party)
	}
	party.results = dlply(vars, level2, partyPlyr, .progress='text')
	attr(party.results, 'formula') <- formula
	attr(party.results, 'level2') <- vars[,level2]
	attr(party.results, 'treatment') <- vars[,as.character(formula)[[2]]]
	attr(party.results, 'covars') <- vars[,!names(vars) %in% c(level2, as.character(formula)[[2]])]
	class(party.results) <- c('mlpsa.tree', 'list')
	return(party.results)
}

#' Returns a data frame with two columns corresponding to the level 2 variable
#' and the leaves from the conditional inference trees.
#' 
#' @seealso mlpsa.ctree
#' @param party.results the results of \code{\link{mlpsa.ctree}}
#' @param data the data frame to merge results to
#' @param level2 the name of the level 2 variable.
#' @return a data frame
#' @export
getStrata <- function(party.results, data, level2) {
	data$strata = as.numeric(NA)
	for(i in names(party.results)) {
		data[which(data[,level2] == i),]$strata = where(party.results[i][[1]], 
									newdata=data[which(data[,level2] == i),])
	}
	return(data)
}
