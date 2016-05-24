#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       simult.pvalue.R                     ####
####                                                 ####
#### FUNCTIONS:  simult.pvalue                       ####
#########################################################

### ======================================
### simult.pvalue
### ======================================
simult.pvalue <- function(sample, precision=0.001, prob=0.95)
{
  PROBS0 <- seq(0.5, 1-precision, by=precision)
  nPROBS0 <- length(PROBS0)
  PROBS <- c(prob, PROBS0)
  CRS <- credible.region(sample, probs=PROBS)

  RES <- list(CR=CRS[[1]], prob=prob)

  contain.zero <- function(MM)
  {
    if (!is.matrix(MM)) stop("MM must be a matrix")
    if (nrow(MM) != 2) stop("MM must have exactly 2 rows")
    DIM <- ncol(MM)
    low.neg <- MM[1,] < 0
    up.pos <- MM[2,] > 0
    contain.ZERO <- low.neg + up.pos
    contain.ZERO <- (sum(contain.ZERO==2) == DIM)

    return(contain.ZERO)
  }

  ## p-value > 0.5?
  is.zero <- contain.zero(CRS[[2]])
  if (is.zero) pval <- ">0.5"
  else{
    ## p.val < precision?
    is.zero <- contain.zero(CRS[[length(CRS)]])
    if (!is.zero) pval <- paste("<", precision, sep="")
    else{
      ## precision <= p-value <= 0.5
      is.zero <- sapply(CRS, contain.zero)
      is.zero <- is.zero[-1]
      first.zero <- nPROBS0 - sum(is.zero) + 1
      pval <- paste(1 - PROBS0[first.zero])
    }
  }

  RES$p.value <- pval
  class(RES) <- "simult.pvalue"
  return(RES)
}  


print.simult.pvalue <- function(x, ...)
{
  level <- paste(x$prob*100, "%", sep="")
  cat("Simultaneous ", level, " rectangular credible region:\n")
  print(x$CR)
  cat("\n")
  cat("Simultaneous Bayesian p-value (based on the smallest rectangle):\n    ", x$p.value, "\n")

  return(invisible(x))
}
