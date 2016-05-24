### Generalized eigen decomposition in double precision.
#
# Call "src/R_dggev.c" wrap "src/dggev.f".
# See "src/dggev.f" for details notations.
#
# The right generalized eigenvector v(j) corresponding to the
# generalized eigenvalue lambda(j) of (A,B) satisfies
#
#   A * v(j) = lambda(j) * B * v(j).
#
# The left generalized eigenvector u(j) corresponding to the
# generalized eigenvalues lambda(j) of (A,B) satisfies
#
#   u(j)**T * A = lambda(j) * u(j)**T * B
#
# where u(j)**T is the transpose of u(j).
#

qz.dggev <- function(A, B, vl = TRUE, vr = TRUE, LWORK = NULL){
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

  JOBVL <- ifelse(vl, "V", "N")
  JOBVR <- ifelse(vr, "V", "N")

  N <- as.integer(ncol(A))
  # S <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 
  # T <- B                           # WCC: memory copy, done in C.
  LDB <- as.integer(N)

  ALPHAR <- double(N)
  ALPHAI <- double(N)
  BETA <- double(N)

  if(vl){
    LDVL <- as.integer(N)
    VL <- double(LDVL * N)
    dim(VL) <- c(LDVL, N)
  } else{
    LDVL <- as.integer(1)
    VL <- double(1)
  }

  if(vr){
    LDVR <- as.integer(N)
    VR <- double(LDVR * N)
    dim(VR) <- c(LDVR, N)
  } else{
    LDVR <- as.integer(1)
    VR <- double(1)
  }

  if(is.null(LWORK) || LWORK < 8 * N){
    LWORK <- as.integer(8 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  BWORK <- integer(1)                # WCC: TODO, no effect if SORT = "N".

  INFO <- integer(1)

  ret <- .Call("R_dggev",
               JOBVL, JOBVR, N,
               A, LDA, B, LDB,
               ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,
               WORK, LWORK, BWORK,
               INFO,
               PACKAGE = "QZ")

  ret$ALPHAR <- ALPHAR
  ret$ALPHAI <- ALPHAI
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

  # Extra returns.
  ret$U <- NULL
  ret$V <- NULL
  if(all(ALPHAI == 0)){
    ret$ALPHA <- ALPHAR
    if(vl){
      ret$U <- VL
    }
    if(vr){
      ret$V <- VR
    }
  } else{
    ret$ALPHA <- complex(real = ALPHAR, imaginary = ALPHAI)
    tmp.id <- matrix(which(ALPHAI != 0), nrow = 2)

    if(vl){
      ret$U <- matrix(as.complex(VL), ncol = N)
      for(i in 1:ncol(tmp.id)){
        tmp <- VL[, tmp.id[, i]]
        ret$U[, tmp.id[, i]] <-
          complex(real = cbind(tmp[, 1], tmp[, 1]),
                  imaginary = cbind(tmp[, 2], -tmp[, 2]))
      }
    }

    if(vr){
      ret$V <- matrix(as.complex(VR), ncol = N)
      for(i in 1:ncol(tmp.id)){
        tmp <- VR[, tmp.id[, i]]
        ret$V[, tmp.id[, i]] <-
          complex(real = cbind(tmp[, 1], tmp[, 1]),
                  imaginary = cbind(tmp[, 2], -tmp[, 2]))
      }
    }
  }

  class(ret) <- "dggev"
  ret
} # End of qz.dggev().


### S3 method.
print.dggev <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("ALPHA:\n")
  print(x$ALPHA)
  cat("\nBETA:\n")
  print(x$BETA)
  if(!is.null(x$U)){
    cat("\nU:\n")
    print(x$U)
  }
  if(!is.null(x$V)){
    cat("\nV:\n")
    print(x$V)
  }

  options(op.org)
  invisible()
} # end of print.dggev().
