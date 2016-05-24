### Generalized eigen decomposition in complex*16 precision.
#
# Call "src/R_zgeev.c" wrap "src/zgeev.f".
# See "src/zgeev.f" for details notations.
#
# The right eigenvector v(j) of A satisfies
#
#                  A * v(j) = lambda(j) * v(j)
#
# where lambda(j) is its eigenvalue.
# The left eigenvector u(j) of A satisfies
#
#               u(j)**H * A = lambda(j) * u(j)**H
#
# where u(j)**H denotes the conjugate transpose of u(j).
#

qz.zgeev <- function(A, vl = TRUE, vr = TRUE, LWORK = NULL){
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

  # Prapare
  JOBVL <- ifelse(vl, "V", "N")
  JOBVR <- ifelse(vr, "V", "N")

  N <- as.integer(ncol(A))
  # T <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 

  W <- complex(N)

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

  RWORK <- double(2 * N)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_zgeev",
               JOBVL, JOBVR, N,
               A, LDA,
               W, VL, LDVL, VR, LDVR,
               WORK, LWORK, RWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$W <- W
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

  class(ret) <- "zgeev"
  ret
} # End of qz.zgeev().


### S3 method.
print.zgeev <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("W:\n")
  print(x$W)
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
} # end of print.zgeev().

