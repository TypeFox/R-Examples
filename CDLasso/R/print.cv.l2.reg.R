print.cv.l2.reg <-
function (x, ...)
{
  cat ("\n Call: \n")
  print (x$call)
  
  cat ("\n Selected Coefficient Estimates: \n")
  print (cbind(Lambda=x$lam.vec,"# Selected"=x$num.pred, "CV Error"=x$mean.error),justify="centre")
  
  cat ("\n Optimal Lambda: \n")
  print (x$lam.opt)

}
