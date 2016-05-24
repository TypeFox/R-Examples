### Generalized eigen decomposition in double precision.
#
# Call "src/R_dgeev.c" wrap "src/dgeev.f".
# See "src/dgeev.f" for details notations.
#
# The right eigenvector v(j) of A satisfies
#
#                  A * v(j) = lambda(j) * v(j)
#
# where lambda(j) is its eigenvalue.
# The left eigenvector u(j) of A satisfies
#
#               u(j)**T * A = lambda(j) * u(j)**T
#
# where u(j)**T denotes the transpose of u(j).
#

qz.dgeev <- function(A, vl = TRUE, vr = TRUE, LWORK = NULL){
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
  JOBVL <- ifelse(vl, "V", "N")
  JOBVR <- ifelse(vr, "V", "N")

  N <- as.integer(ncol(A))
  # T <- A                           # WCC: memory copy, done in C.
  LDA <- as.integer(N) 

  WR <- double(N)
  WI <- double(N)

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

  if(is.null(LWORK) || LWORK < 4 * N){
    LWORK <- as.integer(4 * N)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_dgeev",
               JOBVL, JOBVR, N,
               A, LDA,
               WR, WI, VL, LDVL, VR, LDVR,
               WORK, LWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$WR <- WR
  ret$WI <- WI
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
  if(all(WI == 0)){
    ret$W <- WR
    if(vl){
      ret$U <- VL
    }
    if(vr){
      ret$V <- VR
    }
  } else{
    ret$W <- complex(real = WR, imaginary = WI)
    tmp.id <- matrix(which(WI != 0), nrow = 2)

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

  class(ret) <- "dgeev"
  ret
} # End of qz.dgeev().


### S3 method.
print.dgeev <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("W:\n")
  print(x$W)
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
} # end of print.dgeev().
