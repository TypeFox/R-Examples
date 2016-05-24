dwtest <- function(formula, order.by = NULL, alternative = c("greater", "two.sided", "less"),
  iterations = 15, exact = NULL, tol = 1e-10, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)

  if(!inherits(formula, "formula")) {
    if(!is.null(w <- weights(formula))) {
      if(!isTRUE(all.equal(as.vector(w), rep(1L, length(w)))))
        stop("weighted regressions are not supported")
    }
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  
   
  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  n <- nrow(X)
  if(is.null(exact)) exact <- (n < 100)
  k <- ncol(X)
    
  res <- lm.fit(X,y)$residuals
  dw <- sum(diff(res)^2)/sum(res^2)
  Q1 <- chol2inv(qr.R(qr(X)))
  if(n < 3) {
    warning("not enough observations for computing a p value, set to 1")
    pval <- 1
  } else {
    if(exact)
    {
      A <- diag(c(1,rep(2, n-2), 1))
      A[abs(row(A)-col(A))==1] <- -1
      MA <- diag(rep(1,n)) - X %*% Q1 %*% t(X)
      MA <- MA %*% A
      ev <- eigen(MA)$values[1:(n-k)]
      if(any(Im(ev)>tol)) warning("imaginary parts of eigenvalues discarded")
      ev <- Re(ev)
      ev <- ev[ev>tol]
      pdw <- function(dw) .Fortran("pan", as.double(c(dw,ev)), as.integer(length(ev)),
               as.double(0), as.integer(iterations), x=double(1), PACKAGE = "lmtest")$x
      pval <- switch(alternative,
        "two.sided" = (2*min(pdw(dw), 1-pdw(dw))),
        "less" = (1 - pdw(dw)),
        "greater" = pdw(dw))
  
      if(is.na(pval) || ((pval > 1) | (pval < 0)))
      {
        warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
        exact <- FALSE
      }
    }
    if(!exact)
    {
      if(n < max(5, k)) {
        warning("not enough observations for computing an approximate p value, set to 1")
        pval <- 1        
      } else {
        AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), ncol = k)
        AX[1,] <- X[1,] - X[2,]
        AX[n,] <- X[n,] - X[(n-1),]
        XAXQ <- t(X) %*% AX %*% Q1
        P <- 2*(n-1) - sum(diag(XAXQ))
        Q <- 2*(3*n - 4) - 2* sum(diag(crossprod(AX) %*% Q1)) + sum(diag(XAXQ %*% XAXQ))
        dmean <- P/(n-k)
        dvar <- 2/((n-k)*(n-k+2)) * (Q - P*dmean)
        pval <- switch(alternative,
          "two.sided" = (2*pnorm(abs(dw-dmean), sd=sqrt(dvar), lower.tail = FALSE)),
          "less" = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE),
          "greater" = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
      }
    }
  }
  
  alternative <- switch(alternative,
    "two.sided" = "true autocorrelation is not 0",
    "less" = "true autocorrelation is less than 0",
    "greater" = "true autocorrelation is greater than 0")

  names(dw) <- "DW"
  RVAL <- list(statistic = dw, method = "Durbin-Watson test",
    alternative = alternative, p.value= pval, data.name=dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

