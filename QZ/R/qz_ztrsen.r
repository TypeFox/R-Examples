### Generalized Schur form in complex*16 precision.
#
# Call "src/R_ztrsen.c" wrap "src/ztrsen.f".
# See "src/ztrsen.f" for details notations.
#
#   A = Q * T * (Q)**H = Q_n * T_n * (Q_n)**H
#
# where (Q)**H is the conjugate-transpose of Q.
#

qz.ztrsen <- function(T, Q, select, job = c("B", "V", "E", "N"),
    want.Q = TRUE, LWORK = NULL){
  # Check
  if(!(is.complex(T) && is.complex(Q))){
    stop("T and Q should be in complex.")
  }
  if(!(is.matrix(T) && is.matrix(Q))){
    stop("T and Q should be in matrix.")
  }
  if(any(dim(T) != dim(Q))){
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
  LDA <- as.integer(N) 
  # Q.out <- Q                       # WCC: memory copy, done in C.
  LDQ <- as.integer(N)

  W <- complex(N)

  M <- integer(1)
  S <- double(1)
  SEP <- double(1)

  if(is.null(LWORK) || LWORK < N * (N + 1) / 2){
    LWORK <- as.integer(N * (N + 1) / 2)
  } else{
    LWORK <- as.integer(LWORK)
  }
  WORK <- complex(LWORK)

  INFO <- integer(1)

  # Run
  ret <- .Call("R_ztrsen",
               JOB, COMPQ, SELECT, N, T, LDA,
               Q, LDQ, W, M, S, SEP,
               WORK, LWORK, INFO,
               PACKAGE = "QZ")

  # Return
  ret$W <- W
  ret$M <- M
  ret$S <- S
  ret$SEP <- SEP
  ret$WORK <- WORK[1]
  ret$INFO <- INFO

  # For returns.
  if(! want.Q){
    ret$Q <- NULL
  }

  class(ret) <- "ztrsen"
  ret
} # End of qz.ztrsen().


### S3 method.
print.ztrsen <- function(x, digits = max(4, getOption("digits") - 3), ...){
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
} # end of print.ztrsen().

