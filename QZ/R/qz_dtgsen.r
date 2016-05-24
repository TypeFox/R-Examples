### Generalized Schur form in double precision.
#
# Call "src/R_dtgsen.c" wrap "src/dtgsen.f".
# See "src/dtgsen.f" for details notations.
#
#   (A,B) = Q * (S,T) * (Z)**T = Q_n * (S_n,T_n) * (Z_n)**T
#
# where (Z)**T is the transpose of Z.
#

qz.dtgsen <- function(S, T, Q, Z, select, ijob = 4L,
    want.Q = TRUE, want.Z = TRUE, LWORK = NULL, LIWORK = NULL){
  # Check
  if(!(is.double(S) && is.double(T) && is.double(Q) && is.double(Z))){
    stop("S, T, Q, and Z should be in double")
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

  ALPHAR <- double(N)
  ALPHAI <- double(N)
  BETA <- double(N)

  M <- integer(1)
  PL <- double(1)
  PR <- double(1)
  DIF <- double(2)

  if(is.null(LWORK) || LWORK < max(c(4 * N + 16, N * (N + 1)))){
    LWORK <- as.integer(max(c(4 * N + 16, N * (N + 1))))
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  if(is.null(LIWORK) || LIWORK < max(c(N * (N + 1) / 2, N + 6))){
    LIWORK <- as.integer(max(c(N * (N + 1) / 2, N + 6)))
  } else{
    LIWORK <- as.integer(LIWORK)
  }
  IWORK <- integer(LIWORK)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_dtgsen",
               IJOB, WANTQ, WANTZ, SELECT, N, S, LDA, T, LDB,
               ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
               WORK, LWORK, IWORK, LIWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$ALPHAR <- ALPHAR
  ret$ALPHAI <- ALPHAI
  ret$BETA <- BETA
  ret$M <- M
  ret$PL <- PL
  ret$PR <- PR
  ret$DIF <- DIF
  ret$WORK <- WORK[1]
  ret$IWORK <- IWORK[1]
  ret$INFO <- INFO

  # Extra returns.
  if(all(ALPHAI == 0)){
    ret$ALPHA <- ALPHAR
  } else{
    ret$ALPHA <- complex(real = ALPHAR, imaginary = ALPHAI)
  }
  if(! want.Q){
    ret$Q <- NULL
  }
  if(! want.Z){
    ret$Z <- NULL
  }

  class(ret) <- "dtgsen"
  ret
} # End of qz.dtgsen().


### S3 method.
print.dtgsen <- function(x, digits = max(4, getOption("digits") - 3), ...){
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
} # end of print.dtgsen().

