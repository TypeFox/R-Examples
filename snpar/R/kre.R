kre <-
function(x, y, h, kernel = c("epan","unif","tria","quar",
                                     "triw","tric","gaus","cos"), plot = FALSE)
{
  kernel <- match.arg(kernel)
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (!is.numeric(y))
    stop("'y' must be numeric")
  x <- as.vector(x[is.finite(x)])
  y <- as.vector(y[is.finite(x)])
  m <- length(x)
  n <- length(y)
  if (m != n)
    stop("'x' and 'y' must have the same length")
  if (missing(h))
    h <- as.numeric((quantile(x, 0.75) - quantile(x, 0.25))*n^(-1/5)) 
  f.unif <- function(u) 1/2*(abs(u) <= 1)
  f.tria <- function(u) (1 - abs(u))*(abs(u) <= 1)
  f.epan <- function(u) 3*(1 - u^2)/4*(abs(u) <= 1)
  f.quar <- function(u) 15*(1 - u^2)^2/16*(abs(u) <= 1)
  f.triw <- function(u) 35*(1 - u^2)^3/32*(abs(u) <= 1)
  f.tric <- function(u) 70*(1 - abs(u)^3)^3/81*(abs(u) <= 1)
  f.gaus <- function(u) 1/sqrt(2*pi)*exp(-u^2/2)
  f.cos <- function(u) pi*cos(pi*u/2)/4*(abs(u) <= 1)
  YHAT <- NULL
  u <- outer(x,x, "-")/h
  for (i in 1:nrow(u))
    YHAT[i] <- switch(kernel, unif = sum(f.unif(u[i,])*y)/sum(f.unif(u[i,])), 
                      tria = sum(f.tria(u[i,])*y)/sum(f.tria(u[i,])), 
                      epan = sum(f.epan(u[i,])*y)/sum(f.epan(u[i,])), 
                      quar = sum(f.quar(u[i,])*y)/sum(f.quar(u[i,])),
                      triw = sum(f.triw(u[i,])*y)/sum(f.triw(u[i,])), 
                      tric = sum(f.tric(u[i,])*y)/sum(f.tric(u[i,])), 
                      gaus = sum(f.gaus(u[i,])*y)/sum(f.gaus(u[i,])), 
                      cos = sum(f.cos(u[i,])*y)/sum(f.cos(u[i,])))
  RESULTS <- cbind(x, y, YHAT)  
  results.order <- RESULTS[order(RESULTS[,1]),]
  if (plot)
  {
    plot(x, y, xlab = "explanatory", ylab = "response", 
         main = "kernel regression estimation")
    lines(results.order[,1], results.order[,3], col = 2)
  }
  return(list(results = RESULTS, bw = h)) 
}
