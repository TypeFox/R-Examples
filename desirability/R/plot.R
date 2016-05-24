plot.dMax <- function(x, add = FALSE, nonInform = TRUE, ...)
{
   xRange <- extendrange(c(x$low, x$high))
   if(!add) plot(xRange, c(0, 1), type = "n", xlab = "Input", ylab = "Desirability")
   segments(min(xRange), 0, x$low, 0, ...)
   segments(x$high, 1, max(xRange), 1, ...)
   input <- seq(from = x$low, to = x$high, length = 100)
   output <- predict.dMax(x, input)
   points(input, output, type = "l", ...)
   
   if(nonInform) abline(h = x$missing, lty = 2, ...)
   invisible(x)
}



plot.dBox <- function(x, add = FALSE, nonInform = TRUE, ...)
{
   xRange <- extendrange(c(x$low, x$high))
   if(!add) plot(xRange, c(0, 1), type = "n", xlab = "Input", ylab = "Desirability")
   segments(min(xRange), 0, x$low, 0, ...)
   segments(max(xRange), 0, x$high, 0, ...)   
   segments(x$low, 1, x$low, 0, ...)
   segments(x$high, 1,x$high, 0, ...)
   segments(x$low, 1, x$high, 1, ...)
   
   if(nonInform) abline(h = x$missing, lty = 2, ...)
   invisible(x)
}


plot.dMin <- function(x, add = FALSE, nonInform = TRUE, ...)
{
   xRange <- extendrange(c(x$low, x$high))
   if(!add) plot(xRange, c(0, 1), type = "n", xlab = "Input", ylab = "Desirability")
   segments(min(xRange), 1, x$low, 1, ...)
   segments(x$high, 0, max(xRange), 0, ...)
   input <- seq(from = x$low, to = x$high, length = 100)
   output <- predict.dMin(x, input)
   points(input, output, type = "l", ...)
   if(nonInform) abline(h = x$missing, lty = 2, ...)
   invisible(x)   
}



plot.dTarget <- function(x, add = FALSE, nonInform = TRUE, ...)
{
   xRange <- extendrange(c(x$low, x$high))
   if(!add) plot(xRange, c(0, 1), type = "n", xlab = "Input", ylab = "Desirability")
   segments(min(xRange), 0, x$low, 0, ...)
   segments(x$high, 0, max(xRange), 0, ...)
   input <- seq(from = x$low, to = x$high, length = 100)
   output <- predict.dTarget(x, input)
   points(input, output, type = "l", ...)
   if(nonInform) abline(h = x$missing, lty = 2, ...)
   invisible(x)   
}


plot.dCategorical <- function(x, nonInform = TRUE, ...)
{
   barplot(x$values, ylab = "Desirability", ...)
   if(nonInform) abline(h = x$missing, lty = 2,)
   invisible(x)   
}



plot.dArb <- function(x, add = FALSE, nonInform = TRUE, ...)
{
   xRange <- extendrange(x$x)
   if(!add) plot(xRange, c(0, 1), type = "n", xlab = "Input", ylab = "Desirability")
   input <- seq(from = xRange[1], to = xRange[2], length = 100)
   output <- predict(x, input)
   points(input, output, type = "l", ...)
   if(nonInform) abline(h = x$missing, lty = 2, ...)
   invisible(x)   
   
}


