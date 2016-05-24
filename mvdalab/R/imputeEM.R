imputeEM <- function(data, impute.ncomps = 2, pca.ncomps = 2, CV = TRUE, Init = "mean", 
                  scale = TRUE, iters = 25, tol = .Machine$double.eps^0.25) {
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
  Xa <- dat
  nobj <- nrow(Xa)
  npred <- ncol(Xa)
  Scale <- scale
  Results <- llply(1:impute.ncomps, function(this.pc) {
    Bound <- 1
    Xc <- Xa
    t <- Bound
    Iters <- iters
    Predictions <- matrix(NA, nrow = length(Initials), ncol = Iters)
    Bounds <- matrix(NA, nrow = Iters)
    Predictions[, 1] <- Initials
    Bounds[1, ] <- 1
    for(i in 2:Iters) {
      X <- Xc - rep(colMeans(Xc), each = nobj)
      if (scale) {
        sdscales <- sqrt(colSums((Xc - rep(colMeans(Xc), each = nobj))^2)/(nobj - 1))
        if (any(abs(sdscales) < .Machine$double.eps^0.5)) 
          warning("Scaling with (near) zero standard deviation")
        X <- X/rep(sdscales, each = nobj)
      } else {
        sdscales <- rep(1, npred)
        X <- X
      }
      SVD.m <- svd(X, nu = this.pc, nv = this.pc)
      U <- SVD.m$u 
      V <- as.matrix(SVD.m$v)
      if(this.pc == 1) {
        D <- diag(as.matrix(SVD.m$d[1:this.pc])) 
      } else {
        D <- diag((SVD.m$d[1:this.pc])) 
      }
      dat.pca.a <- U %*% D %*% t(V)
      dat.pca1 <- as.data.frame(sweep(dat.pca.a, 2, sdscales, "*"))
      dat.pca <- as.data.frame(sweep(dat.pca1, 2, colMeans(Xc), "+"))
      names(dat.pca) <- names(X)
      dat.O <- Pre 
      Preds <- as.data.frame(dat.pca)[sapply(dat.O, function(x) is.na(x))]
      dat.O[is.na(dat.O)] <- Preds
      Xc <- dat.O
      Predictions[, i] <- Preds
      Bound <- abs((crossprod(Predictions[, i]) - crossprod(Predictions[, i - 1])) / 
                     crossprod(Predictions[, i - 1]))
      Bounds[i, ] <- Bound
      if(is.na(Bounds[i, ] < tol)) {
        stop("Failure to Converge - No big deal...try again")
      } else if(Bounds[i, ] < tol) break
    }
    if(CV == TRUE) {
      CV.Results <- pcaFit(Xc, ncomp = pca.ncomps, scale = Scale)$GVC
    } else {
      CV.Results <- NULL
    }
    list(Xc, Preds, CV.Results, Initials)
  })
  Imputed.DataFrames <- llply(1:length(Results), function(x) Results[[x]][[1]])
  Imputed.Continous <- llply(1:length(Results), function(x) Results[[x]][[2]])
  if(CV == TRUE) {
    CV.Results <- data.frame(sapply(1:length(Results), function(x) Results[[x]][[3]]))[-1, ]
  } else {
    CV.Results <- NULL
  }
  output <- list(Imputed.DataFrames = Imputed.DataFrames, 
                 Imputed.Continous = Imputed.Continous, 
                 CV.Results = CV.Results, 
                 impute.ncomps = impute.ncomps, 
                 pca.ncomps = pca.ncomps)
  class(output) <- "empca"
  output
}



