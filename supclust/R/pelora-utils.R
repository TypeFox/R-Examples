## The function implementing the sign-flip for Pelora - empirical covariance
sign.change <- function(x,y)
  {
    if(is.null(dx <- dim(x))) stop("'x' must be a numeric matrix")
    signs <- sign(apply(x, 2, cov, y))
    list(x.new = x * rep(signs, each = dx[1]), signs = signs)
  }


## Computing the coefficients for penalized logistic regression
ridge.coef <- function(x,y,lambda)
  {
    X       <- cbind(rep(1,length(y)),x)
    pnlty   <- diag(apply(X,2,var)*lambda*nrow(X),nr=ncol(X))
    th      <- c(log(mean(y)/(1-mean(y))),rep(0,ncol(X)-1))

    for(j in 1:2) {
        p  <- 1 / (1+exp(- drop(X %*% th)))
        W  <- diag(p*(1-p))
        WX <- W %*% X ## <<-- FIXME: make faster (W is diagonal !)
        X.WX <- crossprod(X, WX)
        th <- solve(X.WX + pnlty,
                    crossprod(X, y-p) + X.WX %*% th)
    }
    drop(th)
  }


## The function for the standardization of genes
standardize.genes <- function(exmat)
  {
    means <- apply(exmat,2,mean)
    sdevs <- apply(exmat,2,sd)
    for (i in 1:(dim(exmat)[2]))
      {
        exmat[,i] <- (exmat[,i]-means[i])/sdevs[i]
      }
    list(x=exmat, means=means, sdevs=sdevs)
  }
