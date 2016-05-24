### Generalized Schur form in complex*16 precision.
#
# Call "src/R_zgees.c" wrap "src/zgees.f".
# See "src/zgees.f" for details notations.
#
#   A = Q *T * (Q)**H.
#
# where (Q)**H is the conjugate-transpose of Q.
#

qz.zgees <- function(A, vs = TRUE, LWORK = NULL){
  # Check
  if(!is.complex(A)){
    stop("A should be in complex.")
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
  # T <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 

  SDIM <- as.integer(0L)             # WCC: TODO, no effect if SORT = "N".
  W <- complex(N)

  if(vs){
    LDVS <- as.integer(N)
    VS <- complex(LDVS * N)
    dim(VS) <- c(LDVS, N)
  } else{
    LDVS <- as.integer(1)
    VS <- complex(LDVS)
  }

  if(is.null(LWORK) || LWORK < 2 * N){
    LWORK <- as.integer(2 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- complex(LWORK)

  RWORK <- double(N)

  BWORK <- integer(1)                # WCC: TODO, no effect if SORT = "N".

  INFO <- integer(1)

  # Run
  ret <- .Call("R_zgees",
               JOBVS, SORT, SELECT, N,
               A, LDA, SDIM,
               W, VS, LDVS,
               WORK, LWORK, RWORK, BWORK,
               INFO,
               PACKAGE = "QZ")

  # Return
  ret$W <- W 
  if(vs){
    ret$VS <- VS
  } else{
    ret$VS <- NULL
  }
  ret$WORK <- WORK[1]
  ret$INFO <- INFO

  # Extra returns
  if(vs){
    ret$Q <- VS
  } else{
    ret$Q <- NULL
  }

  class(ret) <- "zgees"
  ret
} # End of qz.zgees().


### S3 method.
print.zgees <- function(x, digits = max(4, getOption("digits") - 3), ...){
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
} # end of print.zgees().

