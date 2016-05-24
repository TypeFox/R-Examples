gqtest <- function(formula, point = 0.5, fraction = 0,
  alternative = c("greater", "two.sided", "less"), order.by = NULL, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)

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

  k <- ncol(X)
  n <- nrow(X)
  if(point > 1) {
    if(fraction < 1) fraction <- floor(fraction * n)
    point1 <- point - ceiling(fraction/2)
    point2 <- point + ceiling(fraction/2 + 0.01)
  } else {
    if(fraction >= 1) fraction <- fraction/n
    point1 <- floor((point-fraction/2) * n)
    point2 <- ceiling((point+fraction/2) * n + 0.01)
  }
  if (point2 > n-k+1 | point1 < k) stop("inadmissable breakpoint/too many central observations omitted")

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

  rss1 <- sum(lm.fit(as.matrix(X[1:point1,]),y[1:point1])$residuals^2)
  rss2 <- sum(lm.fit(as.matrix(X[point2:n,]),y[point2:n])$residuals^2)
  mss <- c(rss1/(point1-k), rss2/(n-point2+1-k))

  gq <- mss[2]/mss[1]
  df <- c(n-point2+1-k, point1-k)
  names(df) <- c("df1", "df2")

  PVAL <- switch(alternative,
    "two.sided" = (2*min(pf(gq, df[1], df[2]), pf(gq, df[1], df[2], lower.tail = FALSE))),
    "less" = pf(gq, df[1], df[2]),
    "greater" = pf(gq, df[1], df[2], lower.tail = FALSE))

  alternative <- switch(alternative,
    "two.sided" = "variance changes from segment 1 to 2",
    "less" = "variance decreases from segment 1 to 2",
    "greater" = "variance increases from segment 1 to 2")

  method <- "Goldfeld-Quandt test"
  names(gq) <- "GQ"
  RVAL <- list(statistic = gq,
      parameter = df,
      method = method,
      p.value= PVAL,
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

