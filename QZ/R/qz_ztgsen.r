### Generalized Schur form in complex*16 precision.
#
# Call "src/R_ztgsen.c" wrap "src/ztgsen.f".
# See "src/ztgsen.f" for details notations.
#
#   (A,B) = Q * (S,T) * (Z)**H = Q_n * (S_n,T_n) * (Z_n)**H
#
# where (Z)**H is the conjugate-transpose of Z.
#

qz.ztgsen <- function(S, T, Q, Z, select, ijob = 4L,
    want.Q = TRUE, want.Z = TRUE, LWORK = NULL, LIWORK = NULL){
  # Check
  if(!(is.complex(S) && is.complex(T) && is.complex(Q) && is.complex(Z))){
    stop("S, T, Q, and Z should be in complex.")
  }
  if(!(is.matrix(S) && is.matrix(T) && is.matrix(Q) && is.matrix(Z))){
    stop("S, T, Q, and Z should be in matrix.")
  }
  if(any((dim(S) != dim(T)) | (dim(S) != dim(Q)) | (dim(S) != dim(Z)))){
    stop("Dimensions of S, T, Q, and Z should be equal.")
  }
  if(dim(S)[1] != dim(S)[2]){
    stop("Squared matrices are required.")
  }
  if(length(select) != ncol(S)){
    stop("The select is incorrect.")
  }
  if(ijob < 0 || ijob > 5){
    stop("The ijob = 0, 1, ..., 5.")
  }

  # Prepare
  IJOB <- as.integer(ijob)
  WANTQ <- as.integer(as.logical(want.Q))
  WANTZ <- as.integer(as.logical(want.Z))
  SELECT <- as.integer(as.logical(select))

  N <- as.integer(ncol(S))
  # S.out <- S                       # WCC: memory copy, done in C.
  LDA <- as.integer(N) 
  # T.out <- T                       # WCC: memory copy, done in C.
  LDB <- as.integer(N)
  # Q.out <- Q                       # WCC: memory copy, done in C.
  LDQ <- as.integer(N)
  # Z.out <- Z                       # WCC: memory copy, done in C.
  LDZ <- as.integer(N)

  ALPHA <- complex(N)
  BETA <- complex(N)

  M <- integer(1)
  PL <- double(1)
  PR <- double(1)
  DIF <- double(2)

  if(is.null(LWORK) || LWORK < N * (N + 1)){
    LWORK <- as.integer(N * (N + 1))
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- complex(LWORK)

  if(is.null(LIWORK) || LIWORK < max(c(N + 2, N * (N + 1) / 2))){
    LIWORK <- as.integer(max(c(N + 2, N * (N + 1) / 2)))
  } else{
    LIWORK <- as.integer(LWORK)
  }
  IWORK <- integer(LIWORK)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_ztgsen",
               IJOB, WANTQ, WANTZ, SELECT, N, S, LDA, T, LDB,
               ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
               WORK, LWORK, IWORK, LIWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$ALPHA <- ALPHA
  ret$BETA <- BETA
  ret$M <- M
  ret$PL <- PL
  ret$PR <- PR
  ret$DIF <- DIF
  ret$WORK <- WORK[1]
  ret$IWORK <- IWORK[1]
  ret$INFO <- INFO

  # For returns.
  if(! want.Q){
    ret$Q <- NULL
  }
  if(! want.Z){
    ret$Z <- NULL
  }

  class(ret) <- "ztgsen"
  ret
} # End of qz.ztgsen().


### S3 method.
print.ztgsen <- function(x, digits = max(4, getOption("digits") - 3), ...){
  op.org <- options()
  options(digits = digits)

  cat("ALPHA:\n")
  print(x$ALPHA)
  cat("\nBETA:\n")
  print(x$BETA)
  cat("\nS:\n")
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
} # end of print.ztgsen().

