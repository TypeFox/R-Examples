kde <-
function(x, h, xgrid, ngrid, kernel = c("epan", "unif","tria","quar",
                                               "triw","tric","gaus","cos"), plot = FALSE)
{
  kernel <- match.arg(kernel)
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric") 
  if (missing(xgrid) && missing(ngrid))
    x.grid <- x[is.finite(x)]
  if (missing(ngrid) && !missing(xgrid))
    x.grid <- xgrid[is.finite(xgrid)]
  if (missing(xgrid) && !missing(ngrid))
    x.grid <- seq(min(x), max(x), length = ngrid)
  if (!missing(xgrid) && !missing(ngrid))
    stop("argument 'xgrid' and 'ngrid' cannot be both specified at the same time")
  if (any(is.na(x)))
    warning("missing values exist and have been removed")
  n <- length(x)
  if (n < 5L)
    stop("not enough 'x' data")
  if (missing(h))
    bw <- as.numeric((8*pi/3)^(1/5)*2.0362*(((quantile(x, 0.75) - 
                      quantile(x, 0.25))/1.349)^(2/3))*n^(-1/5))
  else
  {
    if (!is.numeric(h)) 
      stop("bandwidth 'h' must be numeric")
    else
      bw <- h
  }
  funif <- function(u) 1/2*(abs(u) <= 1)
  ftria <- function(u) (1 - abs(u))*(abs(u) <= 1)
  fepan <- function(u) 3*(1 - u^2)/4*(abs(u) <= 1)
  fquar <- function(u) 15*(1 - u^2)^2/16*(abs(u) <= 1)
  ftriw <- function(u) 35*(1 - u^2)^3/32*(abs(u) <= 1)
  ftric <- function(u) 70*(1 - abs(u)^3)^3/81*(abs(u) <= 1)
  fgaus <- function(u) 1/sqrt(2*pi)*exp(-u^2/2)
  fcos <- function(u) pi*cos(pi*u/2)/4*(abs(u) <= 1)
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
  U <- outer(x.grid, x, "-")/bw
  PDF <- NULL; CDF <- NULL
  for (i in 1:nrow(U))
  {
    PDF[i] <- switch(kernel, unif = mean(funif(U[i,]))/bw, tria = mean(ftria(U[i,]))/bw, 
                     epan = mean(fepan(U[i,]))/bw, quar = mean(fquar(U[i,]))/bw,
                     triw = mean(ftriw(U[i,]))/bw, tric = mean(ftric(U[i,]))/bw, 
                     gaus = mean(fgaus(U[i,]))/bw, cos = mean(fcos(U[i,]))/bw)
    CDF[i] <- switch(kernel, unif = mean(Funif(U[i,])), tria = mean(Ftria(U[i,])), 
                     epan = mean(Fepan(U[i,])), quar = mean(Fquar(U[i,])), 
                     triw = mean(Ftriw(U[i,])), tric = mean(Ftric(U[i,])),
                     gaus = mean(Fgaus(U[i,])), cos = mean(Fcos(U[i,])))
  }                   
  x.PCDF <- cbind(x.grid, PDF, CDF) 
  x.order <- x.PCDF[order(x.PCDF[,1]),]  
  if (plot)
  {
    par(mfrow = c(2,1))
    plot(x.order[,1],x.order[,2], type = "l", xlab = "xgrid",
         ylab = "PDF", main = "kernel CDF estimation")
    plot(x.order[,1],x.order[,3], type = "l", xlab = "xgrid",
         ylab = "CDF", main = "kernel PDF estimation")
    par(mfrow = c(1,1))
  }  
  return(list(data = x, xgrid = x.grid, fhat = PDF, Fhat = CDF, bw = bw)) 
}
