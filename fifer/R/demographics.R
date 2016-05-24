##' Summarize a continuous variable with mean plus/minus standard deviation.
##'
##' Default \code{continuous.summary.function} for use in \code{\link{SummarizeVar}}. Returns formatted text of mean plus/minus standard deviation, possibly by group. For use in construction of demographics tables.
##' @title Summarize a continuous vector with mean plus/minus standard deviation
##' @param x Vector of values.
##' @param group Group identifier to return summaries by group.
##' @param decimal The number of decimal values to format the results; defaults to 2.
##' @param latex Return LaTeX characters if \code{TRUE}; for example, the LaTeX code for the plus-minus symbol.
##' @param na.rm Remove missing values if \code{TRUE}.
##' @param ... Nothing.
##' @return Formatted text of mean plus/minus standard deviation in a vector or matrix.
##' @author Vinh Nguyen
##' @examples
##' SummarizeContinuousDefault(x=c(rnorm(100, 5), rnorm(100, 0)), group=rep(0:1, each=100))
##' @export
SummarizeContinuousDefault <- function(x, group=rep(1, length(x)), decimal=2, latex=TRUE, na.rm=TRUE, ...){
  mu <- formatC(tapply(x, group, mean, na.rm=na.rm), format="f", digits=decimal)
  sd <- formatC(tapply(x, group, sd, na.rm=na.rm), format="f", digits=decimal)
  return(paste(mu, ifelse(latex, "$\\pm$", "\u00B1"), sd, sep=" ")) ## ±
}

##' Summarize a factor variable with count and percentages.
##'
##' Default \code{factor.summary.function} for use in \code{\link{SummarizeVar}}. Returns formatted text of count and percentages. For use in construction of demographics tables.
##' @title Summarize a factor vector with count and percentages
##' @param x Vector of values.
##' @param group Group identifier to return summaries by group.
##' @param decimal The number of decimal values to format the results; defaults to 0.
##' @param latex Return LaTeX characters if \code{TRUE} (default). For example, the LaTeX code for the percentage symbol should be preceeded by the escape character \code{\}}.
##' @param useNA Defaults to \code{ifany} and passed to \code{\link{table}}.
##' @param ... Nothing.
##' @return Formatted text of counts with percentages in parentheses, in a vector or matrix.
##' @author Vinh Nguyen
##' @examples
##' SummarizeFactorDefault(x=c(sample(1:5, 100, replace=TRUE), sample(1:5, 100, replace=TRUE)), 
##'   group=rep(0:1, each=100))
##' @export
SummarizeFactorDefault <- function(x, group=rep(1, length(x)), decimal=0, latex=TRUE, useNA="ifany", ...){
  counts <- table(x, group, useNA=useNA)
  pct <- formatC(counts / matrix(colSums(counts), nrow=nrow(counts), ncol=ncol(counts), byrow=TRUE)*100, format="f", digits=decimal)
  rslt <- matrix(paste(counts, " (", pct, ifelse(latex, "\\%", "%"), ")", sep=""), ncol=ncol(counts))
  rownames(rslt) <- rownames(counts)
  colnames(rslt) <- colnames(counts)
  return(rslt)
}

##' Summarize a continuous or discrete (factor) vector.
##' @title Summarize a vector (continuous or factor)
##' @param x Vector of values.
##' @param group Group identifiers to return summaries by group.
##' @param latex Return LaTeX characters if \code{TRUE} (default). For example, the LaTeX code for the percentage symbol should be preceeded by the escape character \code{\\}.
##' @param decimalFactor The number of decimals to display in percentages for factor variables. This is passed to the \code{decimal} in \code{factor.summary.function}.
##' @param decimalContinuous The number of decimals to display in percentages for numeric variables. This is passed to the \code{decimal} in \code{ContinuousSummaryFunction}.
##' @param ContinuousSummaryFunction Function to use to summarize a continuous variable; defaults to \code{\link{SummarizeContinuousDefault}}. Function must take in the following arguments:
##' \code{x}: a vector of values.
##' \code{group}: a vector that identifies group.
##' \code{decimal}: a numeric value to indicate the decimal places in the formatted output.
##' \code{latex}: a logical value that indicates whether the resulting output contains LaTeX code; should default to \code{TRUE}.
##' \code{...}: additional arguments.
##' @param FactorSummaryFunction Function to use to summarize a factor variable; defaults to \code{\link{SummarizeFactorDefault}}. See \code{ContinuousSummaryFunction}.
##' @param ... Arguments to be passed to \code{ContinuousSummaryFunction} and \code{FactorSummaryFunction}.
##' @return Formatted text in a vector or matrix.
##' @author Vinh Nguyen
##' @references This function was borrowed (and modified) from Vinh Nguyen's \code{day2day} package. 
##' @export
SummarizeVar <- function(x, group=rep(1, length(x)), latex=TRUE, decimalFactor=0, decimalContinuous=2, ContinuousSummaryFunction=SummarizeContinuousDefault, FactorSummaryFunction=SummarizeFactorDefault, ...){
  if(any(is.na(group))) stop("group must not contain NA.")
  if(is.numeric(x)) rslt <- ContinuousSummaryFunction(x=x, group=group, latex=latex, decimal=decimalContinuous, ...)
  else if(is.factor(x) | is.character(x)) rslt <- FactorSummaryFunction(x=x, group=group, latex=latex, decimal=decimalFactor, ...)
  else stop("x needs to be numeric or factor/character.")
  return(rslt)
}

##' Creates a summary table of a data set in a matrix object for pretty printing via \code{\link[xtable]{xtable}}.
##'
##' This is generally used to create demographics table and used with the package \link[xtable]{xtable} to print. To get proper names to display,
##' a \code{data.frame} should be constructed such that the variable names are what the users want to be displayed. For \code{factor} variables,
##' the user should make use of the \code{levels} and \code{labels} arguments in \code{\link{factor}}.
##' @title Summarize a Data Set (Demographics)
##' @param formula A \code{\link[stats]{formula}}, with the left-hand side being empty or a group variable to summarize by. The right-hand side should include variables to summarize by; they should be either continuous variables or factors/characters.
##' @param data A \code{\link{data.frame}} where the variables in \code{formula} come from; if not specified, variables are looked for in the parent environment.
##' @param latex A \link{logical} variable that determines whether the resulting output will be part of a LaTeX document. Defaults to \code{TRUE}.
##' @param na.action A function to handle missing data. See \code{\link[stats]{na.pass}}.
##' @param ... Additional arguments passed to \link{SummarizeVar}, such as \code{decimalFactor}, \code{decimalContinuous}, \code{ContinuousSummaryFunction}, and \code{FactorSummaryFunction}.
##' @return A matrix to be used with \code{\link[xtable]{xtable}}, which in turn should be used in \code{\link[xtable]{print.xtable}}.
##' @author Vinh Nguyen
##' @examples
##' set.seed(1)
##' n <- 50
##' df <- data.frame(trt=sample(0:1, 2*n, replace=TRUE), x1=runif(2*n), x2=rnorm(2*n), 
#'    x3=sample(c("a", "b", "c"), 2*n, replace=TRUE))
##' demographics(~x1+x2+x3, data=df)
##' demographics(trt~x1+x2+x3, data=df)
##' demographics(~., data=df)
##' demographics(trt~., data=df, decimalFactor=2)
##' \dontrun{print(xtable(Summarize(trt~., data=df)), sanitize.text.function=identity)}
##' @export
demographics <- function(formula, data, latex=TRUE, na.action=na.pass, ...) {
  if(missing(data)) {
    data <- environment(formula)
  } else {
    stopifnot(is.data.frame(data))
  }
  mf <- model.frame(formula, data=data, na.action=na.action)
  if(attr(terms(formula, data=mf), "response")) { ## LHS ~ RHS. tells me LHS specified
    group <- factor(model.response(data=mf))##group <- factor(mf[, 1]) ## first column is LHS, assuming something like "x~" and not "x+y~"
    mf <- mf[, -1]
  } else {
    group <- rep(1, nrow(mf))
  }
  if(length(attr(terms(formula, data=mf), "term.labels")) != ncol(mf)) {
    stop("Something wrong with the formula. LHS of '~' should be empty or a single variable.")
  }
  nGroups <- length(unique(group))
  rslt <- NULL
  for(var in colnames(mf)) {
    curr.var.summary <- SummarizeVar(x=mf[, var], group=group, ...)
    if(is.matrix(curr.var.summary)){ ## factor/character variable
      rslt <- rbind(rslt, rep("", nGroups)) ## variable header
      rownames(rslt)[nrow(rslt)] <- var
      rownames(curr.var.summary) <- paste(ifelse(latex, "~~~~", " "), rownames(curr.var.summary), sep="")
      rslt <- rbind(rslt, curr.var.summary)
    } else {
      rslt <- rbind(rslt, curr.var.summary)
      rownames(rslt)[nrow(rslt)] <- var
    }
  }
  colnames(rslt) <- paste(levels(group), " (n=", table(group), ")", sep="")
  if (!latex){
  	rslt = data.frame(rslt)
	names(rslt) = subsetString(names(rslt), "..",position=1)
	for (i in 1:ncol(rslt)){
		rslt[,i] = gsub("$\\pm$", " sd =", x= rslt[,i], fixed=T)	
		rslt[,i] = gsub("\\%", " percent", x= rslt[,i], fixed=T)
	}
  }
  return(rslt)
}
