## extend mean, sd, ... to have a formula interface
## Idea taken from the mosaic package, but put here to
## keep this more self contained

##' Extend function to have formula interface
##'
##' This function wraps a formula interface around a function. It
##' basically calls aggregate and then simplifies the output.
##' @param FUN function to add formula interface for
##' @return a function with a formula interface for the original
## Formulaize <- function(FUN) {
##   flatten <- function(x) if(length(dim(x)) > 2) ftable(x) else x
  
##   function(formula, data, subset, weights, na.action, ...) {
    
##     mf <- match.call(expand.dots = FALSE)
##     m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
##     mf <- mf[c(1L, m)]
##     mf[[1L]] <- as.name("model.frame")
##     d <- eval(mf, parent.frame())
    
##     out <- aggregate(formula, data=d, FUN, ...)
##     flatten(xtabs(formula, out))
##   }
## }

## export these
##' export
##mean.formula   <- Formulaize(mean)
##' export
##median.formula <- Formulaize(median)


