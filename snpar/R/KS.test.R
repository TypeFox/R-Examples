KS.test <-
function(x, y, ..., 
                   kernel = c("epan","unif","tria","quar","triw","tric","gaus","cos"),
                   hx, hy, alternative = c("two.sided", "less", "greater"))
{
  kernel <- match.arg(kernel)
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  if (any(is.na(x)))
    warning("'x' missing values exist and have been removed")
  x <- x[is.finite(x)]
  m <- length(x)
  if (m < 1L)
    stop("not enough 'x' data")
  if (missing(hx))
    hx <- as.numeric((8*pi/3)^(1/5)*2.0362*(((quantile(x, 0.75) - 
                      quantile(x, 0.25))/1.349)^(2/3))*m^(-1/3))
  Funif <- function(u) (1/2 + u/2)*(abs(u) <= 1) + (u > 1)
  Ftria <- function(u) (1/2 + u + u^2/2 - u^2*(u >= 0))*(abs(u) <= 1) + (u > 1) 
  Fepan <- function(u) (1/2 + 3*u/4 - u^3/4)*(abs(u) <= 1) + (u > 1)
  Fquar <- function(u) (3/16*u^5 - 5/8*u^3 + 15/16*u + 1/2)*(abs(u) <= 1) + (u > 1)
  Ftriw <- function(u) 
  {
    (35*u/32 + 1/2 - 5*u^7/32 + 21*u^5/32 - 35*u^3/32)*(abs(u) <= 1) + (u > 1)
  }
  Ftric <- function(u) 
  {
    (1/2 + 70*u/81 + 7*u^10/81 - 14/81*u^10*(u >= 0) + 
       10/27*u^7 + 35/54*u^4 - 35/27*u^4*(u >= 0))*(abs(u) <= 1) + (u > 1)
  }
  Fgaus <- function(u) pnorm(u)
  Fcos <- function(u) (sin(pi/2)/2 + sin(pi*u/2)/2)*(abs(u) <= 1) + (u > 1)
  ux.x <- outer(x, x, "-")/hx 
  kF.x <- NULL 
  for (i in 1:nrow(ux.x))
    kF.x[i] <- switch(kernel, unif = mean(Funif(ux.x[i,])), tria = mean(Ftria(ux.x[i,])), 
                      epan = mean(Fepan(ux.x[i,])), quar = mean(Fquar(ux.x[i,])), 
                      triw = mean(Ftriw(ux.x[i,])), tric = mean(Ftric(ux.x[i,])),
                      gaus = mean(Fgaus(ux.x[i,])), cos = mean(Fcos(ux.x[i,])))
  if (is.numeric(y))
  {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    METHOD <- "Two-sample kernel Kolmogorov-Smirnov test"
    if (any(is.na(y)))
      warning("'y' missing values exist and have been removed")
    y <- y[is.finite(y)]
    n <- length(y)
    if (n < 1L)
      stop("not enough 'y' data")
    if (missing(hy))
      hy <- as.numeric((8*pi/3)^(1/5)*2.0362*(((quantile(x, 0.75) - 
                        quantile(x, 0.25))/1.349)^(2/3))*n^(-1/3))
    ux.y <- outer(y, x, "-")/hy
    uy.x <- outer(x, y, "-")/hx 
    uy.y <- outer(y, y, "-")/hy
    kG.x <- NULL; kG.y <- NULL; kF.y <- NULL
    for (i in 1:nrow(uy.x))
      kG.x[i]<-switch(kernel, unif = mean(Funif(uy.x[i,])), tria = mean(Ftria(uy.x[i,])), 
                      epan = mean(Fepan(uy.x[i,])), quar = mean(Fquar(uy.x[i,])), 
                      triw = mean(Ftriw(uy.x[i,])), tric = mean(Ftric(uy.x[i,])),
                      gaus = mean(Fgaus(uy.x[i,])), cos = mean(Fcos(uy.x[i,])))  
    for (i in 1:nrow(ux.y))
    {
      kG.y[i]<-switch(kernel, unif = mean(Funif(uy.y[i,])), tria = mean(Ftria(uy.y[i,])), 
                      epan = mean(Fepan(uy.y[i,])), quar = mean(Fquar(uy.y[i,])), 
                      triw = mean(Ftriw(uy.y[i,])), tric = mean(Ftric(uy.y[i,])),
                      gaus = mean(Fgaus(uy.y[i,])), cos = mean(Fcos(uy.y[i,])))
      kF.y[i]<-switch(kernel, unif = mean(Funif(ux.y[i,])), tria = mean(Ftria(ux.y[i,])), 
                      epan = mean(Fepan(ux.y[i,])), quar = mean(Fquar(ux.y[i,])), 
                      triw = mean(Ftriw(ux.y[i,])), tric = mean(Ftric(ux.y[i,])),
                      gaus = mean(Fgaus(ux.y[i,])), cos = mean(Fcos(ux.y[i,])))
    }
      
    DN.twosided <- max(max(abs(kF.x - kG.x)), max(abs(kF.y - kG.y)))
    DN.less <- max(max(kF.x - kG.x), max(kF.y - kG.y))
    DN.greater <- max(max(kG.x - kF.x), max(kG.y - kF.y))
    DN <- switch(alternative, two.sided = DN.twosided, less = DN.less, 
                 greater = DN.greater)
    z <- DN*sqrt(m*n/(m + n))  
  }
  else
  {
    if (is.character(y)) 
      CDFx <- get(y, mode = "function", envir = parent.frame())
    if (!is.function(CDFx)) 
      stop("'y' must be numeric or a function or a string naming a valid function")
    METHOD <- "One-sample kernel Kolmogorov-Smirnov test"
    kF.xorder <- cbind(x, kF.x, CDFx(x, ...))
    kF.xorder <- kF.xorder[order(kF.xorder[,1]),]
    Dm.twosided <- max(abs(kF.x - CDFx(x,...))) 
    Dm.greater <- max(kF.x - CDFx(x,...))
    Dm.less <- max(CDFx(x,...) - kF.x)
    Dm <- switch(alternative, two.sided = Dm.twosided, less = Dm.less, 
                 greater = Dm.greater)
    z <- Dm*sqrt(m)
  }
  p.val <- sum((-1)^(c(1:1000) - 1)*exp(-2*c(1:1000)^2*z^2))
  STATISTIC <- ifelse(is.numeric(y), DN, Dm)
  STATISTIC <- switch(alternative, two.sided = setNames(STATISTIC, "D"),
                       less = setNames(STATISTIC, "D^-"),
                       greater = setNames(STATISTIC, "D^+"))
  ALTERNATIVE <- switch(alternative, two.sided = "two.sided",
                        less = "the CDF of x lies below the null hypothesis",
                        greater = "the CDF of x lies above the null hypothesis")
  RVAL <- list(data.name = DNAME, statistic = STATISTIC, p.value = p.val,
               method = METHOD, alternative = ALTERNATIVE)
  class(RVAL) <- "htest"
  return(RVAL)
}
