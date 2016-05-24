### Generalized Schur form in complex*16 precision.
#
# Call "src/R_zgges.c" wrap "src/zgges.f".
# See "src/zgges.f" for details notations.
#
# (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
#
# where (VSR)**H is the conjugate-transpose of VSR.
#

qz.zgges <- function(A, B, vsl = TRUE, vsr = TRUE, LWORK = NULL){
  # Check
  if(!(is.complex(A) && is.complex(B))){
    stop("A and B should be in complex.")
  }
  if(!(is.matrix(A) && is.matrix(B))){
    stop("A and B should be in matrix.")
  }
  if(any(dim(A) != dim(B))){
    stop("Dimensions of A and B should be equal.")
  }
  if(dim(A)[1] != dim(A)[2]){
    stop("Squared matrices are required.")
  }

  # Prepare
  JOBVSL <- ifelse(vsl, "V", "N")
  JOBVSR <- ifelse(vsr, "V", "N")
  SORT <- "N"                        # WCC: TODO
  SELCTG <- 0L                       # WCC: TODO, no effect if SORT = "N".

  N <- as.integer(ncol(A))
  # S <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 
  # T <- B                           # WCC: memory copy, done in C.
  LDB <- as.integer(N)

  SDIM <- as.integer(0L)             # WCC: TODO, no effect if SORT = "N".
  ALPHA <- complex(N)
  BETA <- complex(N)

  if(vsl){
    LDVSL <- as.integer(N)
    VSL <- complex(LDVSL * N)
    dim(VSL) <- c(LDVSL, N)
  } else{
    LDVSL <- as.integer(1)
    VSL <- complex(LDVSL)
  }

  if(vsr){
    LDVSR <- as.integer(N)
    VSR <- complex(LDVSR * N)
    dim(VSR) <- c(LDVSR, N)
  } else{
    LDVSR <- as.integer(1)
    VSR <- complex(LDVSR)
  }

  if(is.null(LWORK) || LWORK < 2 * N){
    LWORK <- as.integer(2 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- complex(LWORK)

  RWORK <- double(8 * N)

  BWORK <- integer(1)                # WCC: TODO, no effect if SORT = "N".

  INFO <- integer(1)

  # Run
  ret <- .Call("R_zgges",
               JOBVSL, JOBVSR, SORT, SELCTG, N,
               A, LDA, B, LDB, SDIM,
               ALPHA, BETA, VSL, LDVSL, VSR, LDVSR,
               WORK, LWORK, RWORK, BWORK,
               INFO,
               PACKAGE = "QZ")

  # Return
  ret$ALPHA <- ALPHA
  ret$BETA <- BETA
  if(vsl){
    ret$VSL <- VSL
  } else{
    ret$VSL <- NULL
  }
  if(vsr){
    ret$VSR <- VSR
  } else{
    ret$VSR <- NULL
  }
  ret$WORK <- WORK[1]
  ret$INFO <- INFO

  # Extra returns
  if(vsl){
    ret$Q <- VSL
  } else{
    ret$Q <- NULL
  }
  if(vsr){
    ret$Z <- VSR
  } else{
    ret$Z <- NULL
  }

  class(ret) <- "zgges"
  ret
} # End of qz.zgges().


### S3 method.
print.zgges <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("ALPHA:\n")
  print(x$ALPHA)
  cat("\nBETA:\n")
  print(x$BETA)
  cat("S:\n")
  print(x$S)
  cat("\nT:\n")
  print(x$T)
  if(! is.null(x$Q)){
    cat("\nQ:\n")
    print(x$Q)
  }
  if(! is.null(x$Z)){
    cat("\nZ:\n")
    print(x$Z)
  }

  options(op.org)
  invisible()
} # end of print.zggev().

