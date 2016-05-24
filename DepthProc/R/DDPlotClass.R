
#' @rdname plot-methods
#' @aliases plot,DDPlot
#' @export
setMethod("plot", signature = c(x = "DDPlot"), function(x){
  p = getPlot(x)
  print(p)
})

#' @rdname getPlot-methods
#' @aliases getPlot,DDPlot
#' @export
setMethod("getPlot", "DDPlot", function(object){
  a_est = data.frame(x = object@X, y = object@Y)
  p = ggplot()
  # eval(as.name("x")) - small hack to fix:
  # getPlot,DDPlot: no visible binding for global variable 'x'
  # getPlot,DDPlot: no visible binding for global variable 'y'
  # i cannot use aes(x,y)
  p = p + geom_point(data = a_est, aes(eval(as.name("x")),eval(as.name("y"))), color = "blue", 
                     shape = 1, size = 3)
  p = p + theme_bw() + .depTheme()
  p = p + ggtitle(object@title)
  p = p + xlab("X depth")
  p = p + ylab("Y depth")
  p = p + ylim(c(0, max(a_est$y)))
  p = p + xlim(c(0, max(a_est$x)))
  p = p + geom_abline(color = "grey")
  return(p)
})

setMethod("show", "DDPlot", function(object){
  cat("DDPlot\n")
  plot(object)
  cat("\nDepth Metohod:\n\t", object@X@method)
})

setGeneric("indexLiu", function(ddplot, gamma) standardGeneric("indexLiu"))
setMethod("indexLiu", signature(ddplot = "DDPlot", gamma = "numeric"),
function(ddplot, gamma)
{
  tmp = abs(as.vector(ddplot@X - ddplot@Y))
  indLiu = sapply(gamma, function(x) sum(tmp>x))
  return(indLiu)
})
