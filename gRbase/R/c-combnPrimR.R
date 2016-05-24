combnPrim <- function(x, m, simplify=TRUE){
  ## FIXME: combnPrim: Could take a FUN argument.
  if (length(x)==1 && is.numeric(x))
    x <- seq(x)
  if (length(x) < m)
    stop("Error in combnPrim: n < m\n")

##   nofun <- is.null(FUN)
##   if (!nofun && !is.function(FUN))
##     stop("'FUN' must be a function or NULL")

  NCAND <- length(x)
  NSEL  <- as.integer(m)
  NSET <- as.integer(choose(NCAND,NSEL))
  ANS  <- rep.int(0L, NSET*NSEL)
  res <- .C("combnC", NSEL, NCAND, NSET, ANS
            ,PACKAGE="gRbase"
  )[[4]]



  if (simplify){
    matrix(x[res], nrow=NSEL, ncol=NSET)
  } else {
    res <- matrix(x[res], nrow=NSEL, ncol=NSET)
    ##res <- split(res, col(res))
    res <- colmat2list( res )
    names(res) <- NULL
    res
  }
}
