plot.geofd <-
function(x, emp.trace.vari=x$emp.trace.vari,
                      trace.vari.array=x$trace.vari.array,
                      colors=rainbow(length(trace.vari.array)), ... )
{

  if(missing(x)) x <- list(emp.trace.vari=emp.trace.vari, trace.vari.array=trace.vari.array)
  if(!is.list(trace.vari.array)) stop("the trace.vari.array parameter must be a list containing Variogram Model elements")
  if(!is.list(emp.trace.vari) || class(emp.trace.vari)!="variogram") stop("the emp.trace.vari parameter must be a list object of the class \"variogram\"")
  if(length(colors)!=length(trace.vari.array)) stop("the arguments \"colors\" and \"trace.vari.array\" must have the same length")

  # The empirical trace variogram is plotted
  plot(emp.trace.vari, xlab="Distance", ylab="Trace-Variogram", main="Trace-Variogram Cloud")

  # The argument legend parameter for the legend is loaded
  legend <- c("empirical trace variogram")

  cont <- 1
  # Each calculated trace variogram is plotted
  for(trace.vari in trace.vari.array){
    lines(trace.vari,col=colors[cont],lwd=2)
    legend <- c(legend, trace.vari$cov.model)
    cont <- cont+1
  }

  # Other arguments for the legend function are loaded
  colors <- c(1, colors)
  pch <- c(21, array(-1,length(trace.vari.array)) )
  lty <- c(0, array(1,length(trace.vari.array)) )
  pt.cex <- c(1, array(1,length(trace.vari.array)) )

  # The legend is plotted
  legend("topleft", "(x,y)", legend, col=colors, pch=pch, lty=lty, pt.cex=pt.cex)
}
