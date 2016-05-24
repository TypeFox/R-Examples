### Generalized Schur form in double precision.
#
# Call "src/R_dgees.c" wrap "src/dgees.f".
# See "src/dgees.f" for details notations.
#
#   A = Q * T * (Q)**T
#
# where (Q)**T is the transpose of Q.
#

qz.dgees <- function(A, vs = TRUE, LWORK = NULL){
  # Check
  if(!is.double(A)){
    stop("A should be in double")
  }
  if(!is.matrix(A)){
    stop("A should be in matrix.")
  }
  if(dim(A)[1] != dim(A)[2]){
    stop("Squared matrices are required.")
  }

  # Prepare
  JOBVS <- ifelse(vs, "V", "N")
  SORT <- "N"                        # WCC: TODO
  SELECT <- 0L                       # WCC: TODO, no effect if SORT = "N".

  N <- as.integer(ncol(A))
  # S <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 

  SDIM <- as.integer(0L)             # WCC: TODO, no effect if SORT = "N".
  WR <- double(N)
  WI <- double(N)

  if(vs){
    LDVS <- as.integer(N)
    VS <- double(LDVS * N)
    dim(VS) <- c(LDVS, N)
  } else{
    LDVS <- as.integer(1)
    VS <- double(1)
  }

  if(is.null(LWORK) || LWORK < 3 * N){
    LWORK <- as.integer(3 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  BWORK <- integer(1)                # WCC: TODO, no effect if SORT = "N".

  INFO <- integer(1)

  # Run
  ret <- .Call("R_dgees",
               JOBVS, SORT, SELECT, N,
               A, LDA, SDIM,
               WR, WI, VS, LDVS,
               WORK, LWORK, BWORK,
               INFO,
               PACKAGE = "QZ")

  # Return
  ret$WR <- WR
  ret$WI <- WI
  if(vs){
    ret$VS <- VS
  } else{
    ret$VS <- NULL
  }
  ret$WORK <- WORK[1]
  ret$INFO <- INFO

  # Extra returns.
  if(all(WI == 0)){
    ret$W <- WR
  } else{
    ret$W <- complex(real = WR, imaginary = WI)
  }
  if(vs){
    ret$Q <- VS
  } else{
    ret$Q <- NULL
  }

  class(ret) <- "dgees"
  ret
} # End of qz.dgees().


### S3 method.
print.dgees <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("W:\n")
  print(x$W)
  cat("\nT:\n")
  print(x$T)
  if(! is.null(x$Q)){
    cat("\nQ:\n")
    print(x$Q)
  }

  options(op.org)
  invisible()
} # end of print.dgees().
