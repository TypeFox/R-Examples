### Generalized eigen decomposition in complex*16 precision.
#
# Call "src/R_zggev.c" wrap "src/zggev.f".
# See "src/zggev.f" for details notations.
#
# The right generalized eigenvector v(j) corresponding to the
# generalized eigenvalue lambda(j) of (A,B) satisfies
#
#   A * v(j) = lambda(j) * B * v(j).
#
# The left generalized eigenvector u(j) corresponding to the
# generalized eigenvalues lambda(j) of (A,B) satisfies
#
#   u(j)**H * A = lambda(j) * u(j)**H * B
#
# where u(j)**H is the conjugate-transpose of u(j).
#

qz.zggev <- function(A, B, vl = TRUE, vr = TRUE, LWORK = NULL){
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
  JOBVL <- ifelse(vl, "V", "N")
  JOBVR <- ifelse(vr, "V", "N")

  N <- as.integer(ncol(A))
  # S <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 
  # T <- B                           # WCC: memory copy, done in C.
  LDB <- as.integer(N)

  ALPHA <- complex(N)
  BETA <- complex(N)

  if(vl){
    LDVL <- as.integer(N)
    VL <- complex(LDVL * N)
    dim(VL) <- c(LDVL, N)
  } else{
    LDVL <- as.integer(1)
    VL <- complex(1)
  }

  if(vr){
    LDVR <- as.integer(N)
    VR <- complex(LDVR * N)
    dim(VR) <- c(LDVR, N)
  } else{
    LDVR <- as.integer(1)
    VR <- complex(1)
  }

  if(is.null(LWORK) || LWORK < 2 * N){
    LWORK <- as.integer(2 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- complex(LWORK)

  RWORK <- double(8 * N)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_zggev",
               JOBVL, JOBVR, N,
               A, LDA, B, LDB,
               ALPHA, BETA, VL, LDVL, VR, LDVR,
               WORK, LWORK, RWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$ALPHA <- ALPHA
  ret$BETA <- BETA
  if(vl){
    ret$VL <- VL
  } else{
    ret$VL <- NULL
  }
  if(vr){
    ret$VR <- VR
  } else{
    ret$VR <- NULL
  }
  ret$WORK <- WORK[1]
  ret$INFO <- INFO

  # Extra returns
  if(vl){
    ret$U <- VL
  } else{
    ret$U <- NULL
  }
  if(vr){
    ret$V <- VR
  } else{
    ret$V <- NULL
  }

  class(ret) <- "zggev"
  ret
} # End of qz.zggev().


### S3 method.
print.zggev <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("ALPHA:\n")
  print(x$ALPHA)
  cat("\nBETA:\n")
  print(x$BETA)
  if(! is.null(x$U)){
    cat("\nU:\n")
    print(x$U)
  }
  if(! is.null(x$V)){
    cat("\nV:\n")
    print(x$V)
  }

  options(op.org)
  invisible()
} # end of print.zggev().

