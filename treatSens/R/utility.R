## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

## use this to produce calls of the form
##  treatSens:::functionName
## so that we can evaluate non-exported functions in
## the user's environment
quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1]] <- as.symbol(":::")
  result[[2]] <- as.symbol("treatSens")
  
  result[[3]] <- if (character.only) name else match.call()[[2]]
  result
}
