reset <- resettest <- function(formula, power=2:3,
  type = c("fitted", "regressor", "princomp"), data=list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    mf <- model.frame(formula)
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), mf)
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(mf)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }
  y <- scale(y, center = attr(terms(mf), "intercept") > 0L)
  k <- ncol(X)
  n <- nrow(X)
  type <- match.arg(type)
  switch(type,

  "fitted" = {
    y.hat <- lm.fit(X,y)$fitted
    Z <- matrix(t(sapply(y.hat, "^", power)), nrow=n)
  },

  "regressor" = {
    Z <- as.matrix(mf[,which(!sapply(mf,is.factor))[-1]])
    Z <- matrix(as.vector(t(sapply(as.vector(Z), "^", power))), nrow=n)
  },

  "princomp" = {
    Z <- as.matrix(mf[,which(!sapply(mf,is.factor))[-1]])
    pc1 <- as.matrix(eigen(cov(Z))$vectors)[,1]
    pc1 <- as.vector(Z %*% pc1)
    Z <- matrix(t(sapply(pc1, "^", power)), nrow=n)
  })

  XZ <- cbind(X, Z)
  q <- ncol(Z)
  res1 <- lm.fit(X,y)$residuals
  res2 <- lm.fit(XZ,y)$residuals
  res1 <- sum(res1^2)
  res2 <- sum(res2^2)
  df1 <- q
  df2 <- n-(k+q)
  reset <- (df2/df1) * ((res1 - res2) / res2)
  names(reset) <- "RESET"
  df <- c(df1, df2)
  names(df) <- c("df1","df2")
  RVAL <- list(statistic = reset,
	parameter = df,
	method = "RESET test",
	p.value= as.vector(pf(reset, df1, df2, lower.tail = FALSE)),
	data.name=dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
