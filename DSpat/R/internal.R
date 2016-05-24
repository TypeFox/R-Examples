print.dspat=function(x,...)
##################################################################################
# Prints relevant parts of a dspat object
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
  cat("\nStudy area: ")
  print(x$study.area)
  cat("\nLines: ")
  print(x$lines.psp)
  cat("\nObservations: ")
  print(x$model$Q$data)
  cat("\nFitted model: ")
  print(x$model)
  return(NULL)
}
summary.dspat=function(object,...)
##################################################################################
# Summarizes ppm model fit
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
  summary(object$model)
}
coef.dspat=function(object,...)
##################################################################################
# Extracts coefficients and separates into intensity and detection parameters
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
  coeff=coef(object$model)
  detection=grep("distance",names(coeff))
  if(length(detection)!=0)
    return(list(intensity=coeff[-detection],detection=coeff[detection]))
  else
    return(list(intensity=coeff,detection=NULL))
}
vcov.dspat=function(object,...)
##################################################################################
# Extracts variance-covariance matrix for coefficients of fitted model
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
  return(vcov.ppm(object$model,gamaction="silent"))
}
AIC.dspat=function(object,...,k)
##################################################################################
# Extracts AIC for fitted model
#
# Jeff Laake
# 9 April 2008
##################################################################################
{
  if(missing(k)) k=2
  return(-2*object$model$maxlogpl+k*length(coef(object$model)))
}

im.clipped=function(x,window)
##################################################################################
# Fills an image with values x into a clipped non-rectangular window and
# returns an image (class im).  An im must be rectagular, so this function
# creates a rectangular image using the bounding box for window with a 0 value
# and then fills in the values (x) that are in window.
#
# Arguments:
#    x      -a vector of image values in order defined by spatstat with
#             y increasing fastest and then x increasing (see rev_val)
#
#   window  - a polygonal window of class owin
#
# Value: an image with NA outside window and image values in window of the
#        image
#
#
#
# Jeff Laake
# 20 April 2008
##################################################################################
{
   x.im=as.im(0,W=window)
   x.im[window]=x
   return(x.im)
}

rev_val=function(x,y,val)
##################################################################################
# Reverses order of vector val such that it matches order needed for
# an image with y increasing within x increasing.  If val aleady in that
# order it will remain unchanged.
#
# Arguments:
#    x      - x coordinates
#    y      - y coordinates
#    val    - values at x,y
#
# Value: reordered vector of values
#
# Jeff Laake
# 20 April 2008
##################################################################################
{
   X=cbind(x,y,val)
   return(X[order(x,y),3])
}

Ops.psp=function(e1,e2)
##################################################################################
# Allows syntax of x==y or x!=y where x and y are 2 psp objects.
# Tests whether 2 line segment objects are = or !=
#
# Arguments:  e1,e2 - psp objects
# Value:  TRUE or FALSE
#
# Jeff Laake
# 14 April 2008
##################################################################################
{
   ok <- switch(.Generic, "==" = , "!=" = TRUE, FALSE)
   if (!ok) {
       warning(.Generic, " not meaningful for psp")
       return(rep.int(NA, max(length(e1), if (!missing(e2)) length(e2))))
   }
   if(!class(e1)[1]=="psp" | !class(e2)[1]=="psp") stop("\nOne or more arguments is not of class psp")
   x.end=endpoints.psp(e1)
   y.end=endpoints.psp(e2)
   if(.Generic == "==")
     return(x.end$n==y.end$n & all(x.end$x==y.end$x) & all(x.end$y==y.end$y) )
   else
     return(!(x.end$n==y.end$n & all(x.end$x==y.end$x) & all(x.end$y==y.end$y)) )
}

owin.gpc.poly=function(window)
##################################################################################
# Converts an owin class composed of a single polygon to a gpc.poly
#
# Arguments:  window  - an owin class
#
# Value    :  gpc.poly from first polygon in owin
#
# Jeff Laake
# 18 April 2008
##################################################################################
{
if(is.null(window$bdry))
  return(as(cbind(c(window$xrange,rev(window$xrange)),
             c(rep(window$yrange[1],2),rep(window$yrange[2],2))),
             "gpc.poly"))
else
  return(as(cbind(window$bdry[[1]]$x,window$bdry[[1]]$y),"gpc.poly"))
}

