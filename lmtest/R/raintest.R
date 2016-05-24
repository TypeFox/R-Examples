raintest <- function(formula, fraction = 0.5, order.by = NULL, center = NULL,
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
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)

  if(is.null(order.by))
  {
    if(is.null(center)) center <- 0.5
    if(center > 1) center <- center/n
    from <- ceiling(quantile(1:n, probs=(center-fraction/2)))
    to <- from + floor(fraction*n) - 1
  }
  else
  if(order.by == "mahalanobis")
  {
    if(is.null(center)) center <- apply(X,2,mean)
    o <- order(mahalanobis(X,center, chol2inv(qr.R(qr(X))), inverted = TRUE))
    X <- as.matrix(X[o,])
    y <- y[o]
    from <- 1
    to <- floor(fraction*n)
  }
  else
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
    if(is.null(center)) center <- 0.5
    if(center > 1) center <- center/n
    from <- ceiling(quantile(1:n, probs=(center-fraction/2)))
    to <- from + floor(fraction*n) - 1
  }

  subX <- as.matrix(X[from:to,])
  suby <- y[from:to]
  n1 <- nrow(subX)
  if(n1 < k) stop("not enough observations in subset")
  resi <- lm.fit(X,y)$residuals
  subresi <- lm.fit(subX, suby)$residuals
  sresi <- sum(resi^2)
  sresi1 <- sum(subresi^2)
  rain <- ((sresi - sresi1)/(n-n1))/(sresi1/(n1-k))
  names(rain) <- "Rain"
  df <- c((n-n1),(n1-k))
  names(df) <- c("df1","df2")
  RVAL <- list(statistic = rain,
      parameter = df,
      method = "Rainbow test",
      p.value= as.vector(pf(rain, df[1], df[2], lower.tail = FALSE)),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

