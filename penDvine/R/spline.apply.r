spline.apply <- function(k,cond.cop.apply,y,xout) {
  return(spline(x=cond.cop.apply[,k],y=y,xout=xout[k])$y)
}
