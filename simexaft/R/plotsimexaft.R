plotsimexaft <-
function(obj,var,extrapolation = c("linear", "quadratic", "both"),ylimit) {

  extrapolation= match.arg(extrapolation)

  x = obj$lambda

  y = apply(obj$theta[[var]],2,mean)

  plot(x,y,type ="p",xlim=c(-1.5,max(x)),ylim=ylimit, xlab=substitute(paste(lambda)), ylab="Estimated Coefficient", main = substitute(paste("Extrapolation Effect on ", var)),col ="green",pch =19)

  new = data.frame(x = seq(-1.2,2,0.2))
  
  if(extrapolation =="linear"){
   
  regl = lm(y~x)
    
  fitted.l = as.numeric(predict(regl, new))
  
  lines(new$x,fitted.l, col="blue")

  points(new$x[2],fitted.l[2],col="red",pch =19)
  
  legend("topright", "Linear Extrapolation", col="blue", lty=1)
  
  }
  
  if( extrapolation=="quadratic"){
  
  regq = lm(y~x+I(x^2))
  
  fitted.q = as.numeric(predict(regq, new))
 
  lines(new$x,fitted.q, col="red")

  points(new$x[2],fitted.q[2],col="blue",pch =19)
  
  legend("topright", "Quadratric Extrapolation", col="red", lty=1)

}

 if(extrapolation =="both")
 {  
  regl = lm(y~x)
    
  fitted.l = as.numeric(predict(regl, new))
  
  lines(new$x,fitted.l, col="blue")

  points(new$x[2],fitted.l[2],col="red",pch =19)
  
  regq = lm(y~x+I(x^2))
  
  fitted.q = as.numeric(predict(regq, new))
 
  lines(new$x,fitted.q, col="red")

  points(new$x[2],fitted.q[2],col="blue",pch =19)
  
  legend("topright", c("Linear Extrapolation", "Quadratric Extrapolation"), col=c("blue","red"), lty=1)
 
 }

}
