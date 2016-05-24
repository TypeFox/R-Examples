### Generalized Schur form in double precision.
#
# Call "src/R_dgges.c" wrap "src/dgges.f".
# See "src/dgges.f" for details notations.
#
#   (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
#
# where (VSR)**T is the transpose of VSR.
#

qz.dgges <- function(A, B, vsl = TRUE, vsr = TRUE, LWORK = NULL){
  if(!(is.double(A) && is.double(B))){
    stop("A and B should be in double")
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
  ALPHAR <- double(N)
  ALPHAI <- double(N)
  BETA <- double(N)

  if(vsl){
    LDVSL <- as.integer(N)
    VSL <- double(LDVSL * N)
    dim(VSL) <- c(LDVSL, N)
  } else{
    LDVSL <- as.integer(1)
    VSL <- double(1)
  }

  if(vsr){
    LDVSR <- as.integer(N)
    VSR <- double(LDVSR * N)
    dim(VSR) <- c(LDVSR, N)
  } else{
    LDVSR <- as.integer(1)
    VSR <- double(1)
  }

  if(is.null(LWORK) || LWORK < 8 * N + 16){
    LWORK <- as.integer(8 * N + 16)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  BWORK <- integer(1)                # WCC: TODO, no effect if SORT = "N".

  INFO <- integer(1)

  ret <- .Call("R_dgges",
               JOBVSL, JOBVSR, SORT, SELCTG, N,
               A, LDA, B, LDB, SDIM,
               ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR,
               WORK, LWORK, BWORK,
               INFO,
               PACKAGE = "QZ")

  ret$ALPHAR <- ALPHAR
  ret$ALPHAI <- ALPHAI
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

  # For returns.
  if(all(ALPHAI == 0)){
    ret$ALPHA <- ALPHAR
  } else{
    ret$ALPHA <- complex(real = ALPHAR, imaginary = ALPHAI)
  }
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

  class(ret) <- "dgges"
  ret
} # End of qz.dgges().


### S3 method.
print.dgges <- function(x, digits = max(4, getOption("digits") - 3), ...){
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
} # end of print.dgges().
