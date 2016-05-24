plot.MICsplines <-
function(x, ...)
{
  mat=x$mat
  xvec=x$x
  matplot(xvec[order(xvec)],mat[order(xvec),],type="l",xlab="",ylab="",main=paste(x$type,"Splines Basis"))
  
}
