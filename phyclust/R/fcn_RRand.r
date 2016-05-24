### This file contains a wrap to call C function in "src/R_RRand.c".
### Written: Wei-Chen Chen on 2008/10/27.


# Call:
#   SEXP R_RRand(SEXP N, SEXP TRUK, SEXP PREDK, SEXP trcl, SEXP prcl)
# Input:
#   N: SEXP[1], number of observations.
#   TRUK: SEXP[1], number of true clusters.
#   PREDK: SEXP[1], number of predicted clusters.
#   trcl: SEXP[N], true cluster ids.
#   prcl: SEXP[N], predicted cluster ids.
# Output in C:
#   ret: a list contains
#      Rand: SEXP[1], Rand index.
#      adjRand: SEXP[1], adjust Rand index.
#      Eindex: SEXP[1], Eindex.
RRand <- function(trcl, prcl, lab = NULL){
  if(! is.null(lab)){
    trcl <- trcl[lab == 0]
    prcl <- prcl[lab == 0]
  }

  N <- length(trcl)

  if(length(trcl) != N || length(prcl) != N){
    stop("The lengths of trcl and prcl do not agree!")
  } 

  tmp.TRUK <- unique(trcl)
  tmp.PREDK <- unique(prcl)
  TRUK <- max(tmp.TRUK)
  PREDK <- max(tmp.PREDK)

  if(min(tmp.TRUK) < 1 || min(tmp.PREDK) < 1){
    stop("The minimum ID is smaller than 1!")
  }
  if(min(tmp.TRUK) != 1 || min(tmp.PREDK) != 1){
    warnings("The minimum ID is not 1!")
  }

  ret <- .Call("R_RRand",
               as.integer(N),
               as.integer(TRUK),
               as.integer(PREDK),
               as.integer(trcl - 1),
               as.integer(prcl - 1),
               PACKAGE = "phyclust")

  class(ret) <- "RRand"
  ret
} # End of RRand().

print.RRand <- function(x, digits = max(4, getOption("digits") - 3), ...){
  RRand <- x
  my.print(unlist(RRand), digits = digits)
} # End of print.RRand().
