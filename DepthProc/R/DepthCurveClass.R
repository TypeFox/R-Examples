


#' @rdname plot-methods
#' @aliases plot,DepthCurve
#' @export
setMethod("plot", signature = c(x = "DepthCurve"), function(x)
{
 plot(new(paste0(class(x),"List"),x))
})

#' @rdname plot-methods
#' @aliases plot,DepthCurveList
#' @export
setMethod("plot", signature = c(x = "DepthCurveList"), function(x)
{
  p = getPlot(x)
  print(p)
})

#############################################################

setMethod("initialize","DepthCurveList", function(.Object, ...)
{
  tmp = list(...)
  n = length(tmp)
  if(n>0)
  {
    .Object[[1]] = tmp[[1]]
    if(n > 1) for(i in 2:length(tmp)) .Object = .Object %+% tmp[[i]]
  }
  return(.Object)
})

#######################################################################


#######################################################################
#' @rdname grapes-plus-grapes-methods
#' @aliases grapes-plus-grapes-methods,DepthCurveList,DepthCurve
#' @export
setMethod("%+%", signature(e1 = "DepthCurveList", e2 = "DepthCurve"), function(e1, e2)
{
  names = sapply(e1,function(x) x@depth@name)
  new_name = e2@depth@name
  if(any(new_name == names))
  {
    warning("Names in DepthCurveList are not unique!")
    k = 1
    new_name_tmp = paste0(new_name)
    while(any(new_name_tmp == names)) { new_name_tmp = paste0(new_name,k); k = k+1}
    e2@depth@name = new_name_tmp
  }
  
  n = length(e1)
  e1[[n+1]] = e2
  return(e1)
})

#' @rdname grapes-plus-grapes-methods
#' @aliases grapes-plus-grapes-methods,DepthCurve,DepthCurveList
#' @export
setMethod("%+%", signature(e1 = "DepthCurve", e2 = "DepthCurveList"), function(e1, e2)
{
  return(e2 %+% e1)
})

#' @rdname grapes-plus-grapes-methods
#' @aliases grapes-plus-grapes-methods,DepthCurve,DepthCurve
#' @export
setMethod("%+%", signature(e1 = "DepthCurve", e2 = "DepthCurve"), function(e1, e2)
{
  return(new(paste0(class(e1),"List"), e1, e2))
})

########################################################################################

setMethod(".getPlot", "DepthCurveList", function(object)
{
  value = unlist(object)
  alpha = as.vector(sapply(object, function(x) x@alpha))
  len_alpha = sapply(object, function(x) length(x@alpha))
  names = sapply(object,function(x) x@depth@name)
  names = rep(names,len_alpha)
  
  data = data.frame(value, alpha, names)
  
  p = ggplot()
  p = p + geom_line(data = data, aes(x = alpha,y = value, col = names), size = 1.5)
  p = p + theme_bw() + .depTheme()
  p = p + ylim(c(0, max(data$value, na.rm=TRUE)))
  p = p + xlim(c(0, max(data$alpha, na.rm=TRUE)))
  return(p)
})

########################

#' @rdname as.matrix-methods
#' @aliases as.matrix,DepthCurveList
#' @export
setMethod("as.matrix", signature(x = "DepthCurveList"), function(x)
{
  names = sapply(x,function(x) x@depth@name)
  tmp = matrix(unlist(x), ncol = length(x))
  colnames(tmp) = names
  tmp
})

setMethod("show", "DepthCurve", function(object)
{
  cat("Object of class:", class(object))
  plot(object)
})

setMethod("show", "DepthCurveList", function(object)
{
  cat("Object of class:", class(object))
  print(getPlot(object))
})
