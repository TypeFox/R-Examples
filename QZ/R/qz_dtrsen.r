### Generalized Schur form in double precision.
#
# Call "src/R_dtrsen.c" wrap "src/dtrsen.f".
# See "src/dtrsen.f" for details notations.
#
#   A = Q * T * (Q)**T = Q_n * T_n * (Q_n)**T
#
# where (Q)**T is the transpose of Q.
#

qz.dtrsen <- function(T, Q, select, job = c("B", "V", "E", "N"),
    want.Q = TRUE, LWORK = NULL, LIWORK = NULL){
  # Check
  if(!(is.double(T) && is.double(Q))){
    stop("T and Q should be in double")
  }
  if(!(is.matrix(T) && is.matrix(Q))){
    stop("T and Q should be in matrix.")
  }
  if(any((dim(T) != dim(Q)))){
    stop("Dimensions of T and Q should be equal.")
  }
  if(dim(T)[1] != dim(T)[2]){
    stop("Squared matrices are required.")
  }
  if(length(select) != ncol(T)){
    stop("The select is incorrect.")
  }
  if(!(job[1] %in% c("B", "V", "E", "N"))){
    stop("The job should be one of \"B\", \"V\", \"E\", or \"N\".")
  }

  # Prepare
  JOB <- as.character(job[1])
  COMPQ <- ifelse(want.Q, "V", "N")
  SELECT <- as.integer(as.logical(select))

  N <- as.integer(ncol(T))
  # T.out <- T                       # WCC: memory copy, done in C.
  LDB <- as.integer(N)
  # Q.out <- Q                       # WCC: memory copy, done in C.
  LDQ <- as.integer(N)

  WR <- double(N)
  WI <- double(N)

  M <- integer(1)
  S <- double(1)
  SEP <- double(1)

  if(is.null(LWORK) || LWORK < N * (N + 1) / 2){
    LWORK <- as.integer(N * (N + 1) / 2)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- double(LWORK)

  if(is.null(LIWORK) || LIWORK < N * (N + 1) / 4){
    LIWORK <- as.integer(N * (N + 1) / 4)
  } else{
    LIWORK <- as.integer(LIWORK)
  }
  IWORK <- integer(LIWORK)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_dtrsen",
               JOB, COMPQ, SELECT, N, T, LDB, Q, LDQ,
               WR, WI, M, S, SEP,
               WORK, LWORK, IWORK, LIWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$WR <- WR
  ret$WI <- WI
  ret$M <- M
  ret$S <- S
  ret$SEP <- SEP
  ret$WORK <- WORK[1]
  ret$IWORK <- IWORK[1]
  ret$INFO <- INFO

  # Extra returns.
  if(all(WI == 0)){
    ret$W <- WR
  } else{
    ret$W <- complex(real = WR, imaginary = WI)
  }
  if(! want.Q){
    ret$Q <- NULL
  }

  class(ret) <- "dtrsen"
  ret
} # End of qz.dtrsen().


### S3 method.
print.dtrsen <- function(x, digits = max(4, getOption("digits") - 3), ...){
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
} # end of print.dtrsen().

