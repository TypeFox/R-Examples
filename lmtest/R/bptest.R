bptest <- function(formula, varformula = NULL, studentize = TRUE,
  data = list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
	   formula$x
	 else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
	   formula$y
	 else model.response(model.frame(formula))
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
  }  

  ## only use complete cases that are in both models
  if(!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% row.names(X))))) {
    allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
    X <- X[allnames,]
    Z <- Z[allnames,]
    y <- y[allnames]
  }
   
  k <- ncol(X)
  n <- nrow(X)

  resi <- lm.fit(X,y)$residuals
  sigma2 <- sum(resi^2)/n

  if(studentize)
  {
    w <- resi^2 - sigma2
    fv <- lm.fit(Z,w)$fitted
    bp <- n * sum(fv^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
  }
  else
  {
    f <- resi^2/sigma2 -1
    fv <- lm.fit(Z,f)$fitted
    bp <- 0.5 * sum(fv^2)
    method <- "Breusch-Pagan test"
  }

  names(bp) <- "BP"
  df <- ncol(Z)-1
  names(df) <- "df";
  RVAL <- list(statistic = bp,
      parameter = df,
      method = method,
      p.value= pchisq(bp,df,lower.tail=FALSE),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

