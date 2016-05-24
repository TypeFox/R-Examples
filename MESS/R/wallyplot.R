wallyplot.default <- function(x, y=x, FUN=residualplot,
                              hide=TRUE,
                              simulateFunction=rnorm,
                              ...) {

  simulateFunction <- match.fun(simulateFunction)

  if (is.vector(y) && length(x)==length(y))
    y <- cbind(y, sapply(1:8, function(k) {simulateFunction(length(x))}))
  
  if (!is.numeric(x) || !is.matrix(y))
    stop("x and y input must be a vector and matrix, respectively")

  if (length(x) != nrow(y))
    stop("x and y does not conform")

  if (ncol(y) != 9)
    stop("y must be a matrix with 9 columns")
  
  if (!is.numeric(x) && !is.numeric(y)) {
    stop("Must have a pair of numeric vectors or an lm object as input")
  }
  
  cc <- complete.cases(x, y)
  x <- x[cc]
  y <- y[cc,]

  plot3x3(x,y, FUN=FUN, hide=hide, ...)

  invisible(NULL)
}


wallyplot.lm <- function(x, y=x, FUN=residualplot,
	                 hide=TRUE,
                         simulateFunction=rnorm,
                         ...) {

  # Extract information from model fit 
  y <- rstandard(x)
  x <- predict(x)

  wallyplot.default(x, y, FUN=FUN, hide=hide, simulateFunction=simulateFunction, ...)
}


wallyplot <- function(x, y=x, FUN=residualplot,
                      hide=TRUE,
                      simulateFunction=rnorm,
                       ...) {
  UseMethod("wallyplot")
}


# qqnorm.wally <- function(x, y, ...) {qqnorm(y, ...) ; abline(a=0, b=1)}


plot3x3 <- function(x, y, FUN=plot, hide=TRUE, ylim=range(y), mar=c(4, 4, .1, .1)+.1, ...) {
  # Input check

  if (!is.numeric(x) || !is.matrix(y))
    stop("x and y input must be a vector and matrix, respectively")

  if (length(x) != nrow(y))
    stop("x and y does not conform")

  if (ncol(y) != 9)
    stop("y must be a matrix with 9 columns")
         
  oldpar <- par(no.readonly = TRUE)
  par(mar=mar)

  FUN <- match.fun(FUN)
  
  pos <- c(1:9)
  if (hide)
    pos <- sample(pos)

  par(mfrow=c(3,3))
  
  for (i in 1:length(pos)) {
    FUN(x, y[,pos[i]], ylim=ylim, ...)
  }
  
  if (hide) {
    readline("Hit <Enter> to show the original plot. ")
  }
  
  figpos <- order(pos)[1]
  par(mfg=c(figpos %/% 3.1 + 1, figpos - (figpos %/% 3.1)*3 ))
  box(col="red", lwd=2)
  
  par(oldpar)
}
