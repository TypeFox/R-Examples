SeqimputeEM <- function(data, max.ncomps = 5, max.ssq = 0.99, Init = "mean", 
                     adjmean = FALSE, max.iters = 200, tol = .Machine$double.eps^0.25) {
  dat.import <- data
  if((nrow(na.omit(dat.import)) == nrow(dat.import))) {
    stop("No Missing Values Detected")
  } else {
    dat.import
  }
  Scale <- scale
  dat.num.check <- dat.import[(sapply(dat.import, function(x) is.numeric(x)) == T)]
  if(ncol(dat.num.check) == 0) {
    stop("No continous variables, use rough.fix")
  } else {
    Start.Imputation <- imputeRough(dat.import, Init = Init)
    Pre <- Start.Imputation$Pre.Imputed
    Initials <- Start.Imputation$Initials
    dat <- Start.Imputation$Imputed.Dataframe
  }
  Mindex = is.na(Pre)
  Imputed.DataFrames = NULL
  Dindex <- !is.na(Pre)
  X <- Start.Imputation$Imputed.Dataframe
  nobj <- nrow(X)
  npred <- ncol(X)
  Xm <- colMeans(X)
  Xc <- X - rep(Xm, each = nobj)
  TSSQ <- crossprod(Xc[Dindex])
  SSQ <- 0
  a <- 0
  while (a < max.ncomps  & SSQ < max.ssq) {
    a <- a + 1
    iter <- 1
    ftol <- 2*tol; iter <- 1;
    while (iter < max.iters & ftol > tol){
      SVD.m <- svd(Xc, nu = a, nv = a)
      if(a == 1) {
        D <- diag(as.matrix(SVD.m$d[1:a])) 
      } else {
        D <- diag((SVD.m$d[1:a])) 
      }
      U <- SVD.m$u 
      V <- as.matrix(SVD.m$v)
      Xpredicted <- U %*% D %*% t(V)
      SSold <- crossprod(Xc[Mindex])
      ftol <- abs(crossprod(Xpredicted[Mindex])-SSold)/SSold;
      Xc[Mindex] <- Xpredicted[Mindex]
      if (adjmean) { 
        X <- Xc + rep(Xm, each = nobj)
        Xm <- colMeans(X);
        Xc <- X - rep(Xm, each = nobj)
      }
        iter <- iter + 1
    }
      SSres <- crossprod(Xc[Dindex])-crossprod(Xpredicted[Dindex]);
      SSQ <- 1-SSres/TSSQ
  }
  X <- Xc + rep(Xm, each = nobj)
  output <- list(Imputed.DataFrames = X, ncomps = a)
  class(output) <- "print.seqem"
  output
}
