plot_fitted <- function(x, y, object, xlab, ylab, main, ...) 
{
  name_x <-deparse(substitute(x))
  name_y <-deparse(substitute(y))
  
  if(missing(xlab))
    xlab <- name_x
  if(missing(ylab))
    ylab <- name_y 
  if(missing(main))
    main <- "Spatial Stochastic Frontier Analysis" 

    data=data.frame(x,y,fitted.ssfa(object))
    names(data)=c("x","y","fitted")
    plot(x=data$x, y=data$y, xlab=xlab, ylab=ylab, main=main, ...)
    lines(sort(data$x), data$fitted[order(data$x)], ...)

}

