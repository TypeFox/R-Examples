#' Student's T Quantile - Quantile Plot
#'
#' Creates emperical QQ-plot of the quantiles of the data set x versus of a t 
#' distribution. The QQ-plot can be used to determine whether the sample in x is
#' drawn from a t distribution with specified number of degrees of freedom.
#'
#' @param Ra Sample data set
#' @param df Number of degrees of freedom of the t distribution
#' 
#' @references Dowd, K. Measuring Market Risk, Wiley, 2007.
#' 
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # t-QQ Plot for randomly generated standard normal data
#'    Ra <- rnorm(100)
#'    TQQPlot(Ra, 20)
#'
#' @export
TQQPlot<- function(Ra, df){
  
  x <- as.vector(Ra)
  
  if(!is.vector(Ra)){
    stop("The first argument should be a vector.")
  }
  
  mu <- mean(x)
  sigma <- sd(x)
  y<- sort(x)
  a <- PlotPos(y)
  x <- a$pp
  n <- a$n
  x <- mu + sigma * qt(x, df) # This is theoretical t-quantile
  xx <- x
  yy <- y
  
  # Dowd's code has following details but since his code does not work for 
  # matrices with rows, columns > 1, it is not efficient to have it here.
  # if(((dim(x)[1] == n) | (dim(x)[1] == 1 & dim(x)[2] == n)) & ~any(is.nan(x))){
  #  xx <- sort(x)
  # } else {
  #   xx <- quantile(x, pvec)[[1]]
  # }
  #
  # if(((dim(y)[1] == n) | (dim(y)[1] == 1 & dim(y)[2] == n)) & ~any(is.nan(y))){
  #   yy <- sort(y)
  # } else {
  #   yy <- quantile(y, pvec)[[1]]
  # }
  xx <- sort(x)
  yy <- sort(y)
  
  q1x <- quantile(x, .25)[[1]]
  q3x <- quantile(x, .75)[[1]]
  q1y <- quantile(y, .25)[[1]]
  q3y <- quantile(y, .75)[[1]]
  qx <- matrix(c(q1x,q3x), 2, length(q1x))
  qy <- matrix(c(q1y,q3y), 2, length(q1y))
  
  dx <- q3x - q1x
  dy <- q3y - q1y
  slope <- dy/dx
  centerx <- (q1x + q3x)/2
  centery <- (q1y + q3y)/2
  maxx <- max(x)
  minx <- min(x)
  maxy <- centery + slope * (maxx - centerx)
  miny <- centery - slope * (centerx - minx)
  
  mx <- matrix(c(minx,maxx), 2, length(minx))
  my <- matrix(c(miny,maxy), 2, length(miny))
  
  xmin <- min(xx, qx, mx)
  xmax <- max(xx, qx, mx)
  ymin <- min(yy, qy, my)
  ymax <- max(yy, qy, my)

  plot(xx, yy, type = "p", pch=3, col="red", xlab = "t-Quantiles", 
       ylab = "Quantiles of Input Sample", 
       main = paste("QQ Plot of Sample Data versus Student-t with ", df,
                    "Degrees of freedom"),
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  par(new = TRUE)
  plot(qx, qy, type = "l", col="blue", xlab = "t-Quantiles", 
       ylab = "Quantiles of Input Sample", 
       main = paste("QQ Plot of Sample Data versus Student-t with ", df,
                    "Degrees of freedom"),
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  par(new = TRUE)
  plot(mx, my, type = "l", xlab = "t-Quantiles", 
       ylab = "Quantiles of Input Sample", 
       main = paste("QQ Plot of Sample Data versus Student-t with ", df,
                    "Degrees of freedom"),
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
} 

# Helper Functions

# Position PLot
PlotPos <- function(Ra){
  # 
  sx <- as.matrix(Ra)
  if (!is.matrix(sx)) {
    stop("Input should be a matrix.")
  }
  n <- dim(sx)[1]
  m <- dim(sx)[2]
  if (n == 1){
    sx <- t(sx)
    n = m
    m = 1
  }
  nvec <- sum(!is.nan(sx))
  pp <- RepMat(as.matrix(1:n), 1, m)
  pp <- (pp - .5)/ RepMat(nvec, n, 1)
  pp[is.nan(sx)] <- NaN
  
  if (nargs() > 1){
    n <- max(nvec)
  }
  return(list("pp" = pp, "n" = n))
  
}

# Implementation of repmat from matlab in R
RepMat <- function(X,m,n){
  X <- as.matrix(X)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  a <- matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
  return(a)
}