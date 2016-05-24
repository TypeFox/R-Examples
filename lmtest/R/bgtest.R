bgtest <- function(formula, order = 1, order.by = NULL, type = c("Chisq", "F"), data = list(), fill = 0)
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
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
  k <- ncol(X)
  order <- 1:order
  m <- length(order)
  resi <- lm.fit(X,y)$residuals

  Z <- sapply(order, function(x) c(rep(fill, length.out = x), resi[1:(n-x)]))
  if(any(na <- !complete.cases(Z))) {
    X <- X[!na, , drop = FALSE]
    Z <- Z[!na, , drop = FALSE]
    y <- y[!na]
    resi <- resi[!na]
    n <- nrow(X)
  }
  auxfit <- lm.fit(cbind(X,Z), resi)

  cf <- auxfit$coefficients
  vc <- chol2inv(auxfit$qr$qr) * sum(auxfit$residuals^2) / auxfit$df.residual
  names(cf) <- colnames(vc) <- rownames(vc) <- c(colnames(X), paste("lag(resid)", order, sep = "_"))

  switch(match.arg(type),

  "Chisq" = {
    bg <- n * sum(auxfit$fitted^2)/sum(resi^2)
    p.val <- pchisq(bg, m, lower.tail = FALSE)
    df <- m
    names(df) <- "df"
  },

  "F" = {
    uresi <- auxfit$residuals
    bg <- ((sum(resi^2) - sum(uresi^2))/m) / (sum(uresi^2) / (n-k-m))
    df <- c(m, n-k-m)
    names(df) <- c("df1", "df2")
    p.val <- pf(bg, df1 = df[1], df2 = df[2], lower.tail = FALSE)
  })

  names(bg) <- "LM test"
  RVAL <- list(statistic = bg, parameter = df,
               method = paste("Breusch-Godfrey test for serial correlation of order up to", max(order)),
               p.value = p.val,
               data.name = dname,
	       coefficients = cf,
	       vcov = vc)
  class(RVAL) <- c("bgtest", "htest")
  return(RVAL)
}

vcov.bgtest <- function(object, ...) object$vcov
df.residual.bgtest <- function(object, ...) if(length(df <- object$parameter) > 1L) df[2] else NULL

