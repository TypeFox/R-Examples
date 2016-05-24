
#' converts to character with minimal loss of precision for numeric variables
#' @param x the value
#' @param ... passed to either format or as.character.
#' @export
as.char <- function(x,...){
	if(is.numeric(x))
		if(is.integer(x))
			ifelse(is.na(x),NA,format(x,trim=TRUE,scientific=FALSE,...))
		else
			ifelse(is.na(x),NA,format(x,trim=TRUE,scientific=FALSE,digits=15,...))
	else
		as.character(x,...)
}


orWaldCI <- function(tables, alpha = 0.05){
# a = extract.counts(tables[[1]])
  a = tables[[1]][[1]]$table
  n00 = a[1,1]
  n01 = a[1,2]
  n10 = a[2,1]
  n11 = a[2,2]
  #
  #  Compute the odds ratio between two binary variables, x and y,
  #  as defined by the four numbers nij:
  #
  #    n00 = number of cases where x = 0 and y = 0
  #    n01 = number of cases where x = 0 and y = 1
  #    n10 = number of cases where x = 1 and y = 0
  #    n11 = number of cases where x = 1 and y = 1
  #
  OR <- (n00 * n11)/(n01 * n10)
  #
  #  Compute the Wald confidence intervals:
  #
  siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
  zalph <- qnorm(1 - alpha/2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  #
  ORlo <- exp(loglo)
  ORhi <- exp(loghi)
  #
# oframe <- data.frame(LowerCI = ORlo, OR = OR, UpperCI = ORhi, alpha = alpha)
  cat(sprintf("The odds ratio is %f. A 95%% confidence interval is (%f, %f).\n",
     OR,ORlo,ORhi))
  invisible()
}

.catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 
