#' compute and plot 2 dimensional histogram of ff data
#'
#' function interface modeled after gplots::hist2d
#' @keywords internal
hist2d.ff <- function(x,y=NULL, nbins=100, show=TRUE){
  if (length(nbins) == 1){
    nbins <- rep(nbins, 2)
  }
  
  stopifnot(length(x) == length(y))
  
  x_range <- range(x, na.rm=TRUE)
  y_range <- range(y, na.rm =TRUE)
  
  x_n <- nbins[1]
  x_steps <- seq(from=x_range[1], to=x_range[2], length.out=(x_n+1)) 
  
  y_n <- nbins[2]
  y_steps <- seq(from=y_range[1], to=y_range[2], length.out=(y_n+1))
  
  h <- ffdf(x=x, y=y)
  h <- transform( h 
           , x_bin = findInterval(x, x_steps, rightmost.closed=T)
           , y_bin = findInterval(y, y_steps, rightmost.closed=T)
  )
  
  m <- binned_tabulate(h$y_bin, h$x_bin, nbins=x_n, nlevels=y_n)
  m <- m[,-1] # drop na column
  if (show){
    image(x_steps, y_steps, m)
  }
  structure(
    list( counts=m
        , x.breaks = x_steps
        , y.breaks = y_steps
        , nobs = length(x)
        , call = match.call()
        ), class="hist2d")
}

# # testing
# iris_f <- as.ffdf(iris)
# hist2d.ff(iris_f$Sepal.Length, iris_f$Sepal.Width, nbins=10)
