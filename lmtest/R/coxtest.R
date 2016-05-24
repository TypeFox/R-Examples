## see Greene (2003), Section 8.3.4, p.155

coxtest <- function(formula1, formula2, data = list())
{
  ## merge two models (if possible) and compute
  ## response y and regressor matrices X and Z
  ## 1. If formulas are available: build big model first
  if(inherits(formula1, "formula") && inherits(formula2, "formula")) {
    formula <- formula2
    if(length(formula) > 2) formula[[2]] <- NULL
    formula[[2]] <- as.call(list(as.name("+"), as.name("."), formula[[2]]))
    formula <- update(formula1, formula)
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(terms(formula1), mf)
    Z <- model.matrix(terms(formula2), mf)
    m1 <- deparse(formula1)
    m2 <- deparse(formula2)
  } else {
  ## 2. if not, go via row names
    if(!inherits(formula1, "formula")) {
      mf <- model.frame(formula1)
      X <- if(is.matrix(formula1$x)) formula1$x else model.matrix(terms(formula1), mf)
      y <- if(is.vector(formula1$y)) formula1$y else model.response(mf)
      m1 <- deparse(formula(formula1))
    } else {
      mf <- model.frame(formula1, data = data)
      y <- model.response(mf)
      X <- model.matrix(formula1, data = data)
      m1 <- deparse(formula1)
    }  

    if(!inherits(formula2, "formula")) {
      mf <- model.frame(formula2)
      Z <- if(is.matrix(formula2$x)) formula2$x else model.matrix(terms(formula2), mf)
      y2 <- if(is.vector(formula2$y)) formula2$y else model.response(mf)
      m2 <- deparse(formula(formula2))
    } else {
      mf <- model.frame(formula2, data = data)
      y2 <- model.response(mf)
      Z <- model.matrix(formula2, data = data)
      m2 <- deparse(formula2)
    }  

    if(!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% row.names(X))))) {
      warning("models fitted on different subsets")
      allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
      X <- X[allnames,]
      Z <- Z[allnames,]
      y <- y[allnames]
      y2 <- y2[allnames]
    }
    if(!identical(y, y2)) warning("different dependent variables specified")
  }
  ## for pretty printing
  m1 <- paste(m1, collapse = "\n")
  m2 <- paste(m2, collapse = "\n")

  ## Steps in Greene (2003) p.156-7
  # 1.
  x.fit <- lm.fit(X, y)
  s2x <- mean(x.fit$residuals^2)
  
  # 2.
  z.fit <- lm.fit(Z, y)  
  s2z <- mean(z.fit$residuals^2)

  # 3. / 5.
  zx.fit <- lm.fit(Z, x.fit$fitted)
  s2zx <- s2x + mean(zx.fit$residuals^2)
  xz.fit <- lm.fit(X, z.fit$fitted)
  s2xz <- s2z + mean(xz.fit$residuals^2)
  
  # 4.
  xzx.fit <- lm.fit(X, zx.fit$residuals)
  s2xzx <- mean(xzx.fit$residuals^2)
  zxz.fit <- lm.fit(Z, xz.fit$residuals)
  s2zxz <- mean(zxz.fit$residuals^2)
  
  # 6.
  n <- nrow(X)
  c01 <- n/2 * log(s2z/s2zx)
  v01 <- n * (s2x/(s2zx^2)) * s2xzx
  c10 <- n/2 * log(s2x/s2xz)
  v10 <- n * (s2z/(s2xz^2)) * s2zxz

  ## Cox statistic for nonnested models
  stat01 <- c01/sqrt(v01)
  stat10 <- c10/sqrt(v10)
  rval <- cbind(c(c01, c10), sqrt(c(v01, v10)), c(stat01, stat10), 
                2*pnorm(abs(c(stat01, stat10)), lower.tail=FALSE))
  colnames(rval) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(rval) <- c("fitted(M1) ~ M2", "fitted(M2) ~ M1")

  ## put results together
  title <- "Cox test\n"
  topnote <- paste("Model ", 1:2,": ", c(m1, m2), sep="", collapse="\n")
  rval <- structure(as.data.frame(rval), heading = c(title, topnote),
	            class = c("anova", "data.frame"))
  return(rval)
}
