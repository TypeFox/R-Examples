##############################################################################################################
##########Obtain the variance estimation and CI for t.minus, t.plus and youden under normality assumption
##########if Youden index is estimated by "Normal"
##Note:if also use the function for TN, var on youden is fine but the variance on t.minus and t.plus are for t.minus.TN and t.plus.TN on the after Box-Cox trasnformation scale not directly for t.minus and t.plus on the orignal data scale
##############################################################################################################

Youden3Grp.Variance.Normal <- function(x,y,z,alpha=0.05)
{
  

  ########Changed from name"Normal.Youden.Cutoff.CI" to "Youden3Grp.Variance.Normal"
  ####x,y,z: samples from D-,D0,D+ groups
  ####alpha: (1-alpha)*100% CI
  
  
  temp.res <- Youden3Grp.PointEst(x,y,z, method="Normal",randomStart.N=1)$est

  
  ###obtain estimated mu and sd from the above function, for TN, these are on the after Box-Cox transformation scale, for Normal, on original scale

  n.minus <- length(na.exclude(x))
  n0 <- length(na.exclude(y))
  n.plus <- length(na.exclude(z))
  
  mu.minus <- temp.res$mu.minus
  mu0 <- temp.res$mu0
  mu.plus <- temp.res$mu.plus

  s.minus <- temp.res$s.minus
  s0 <- temp.res$s0
  s.plus <- temp.res$s.plus

  t.minus <- temp.res$t.minus
  t.plus <- temp.res$t.plus
  youden <- temp.res$youden

  #if(method=="Normal")
  #  {
  
  #t.minus <- temp.res$t.minus
  #t.plus <- temp.res$t.plus
  #youden <- temp.res$youden
  
  #  }
  #else if(method=="TN")
  #  {
  #    t.minus <- temp.res$t.minus.TN
  #    t.plus <- temp.res$t.plus.TN
  #    youden <- temp.res$youden
  #  }
  #else stop("method can only be Normal/TN!")

  var.t.minus <- Var.t.minus(mu.minus,mu0,s.minus,s0,n.minus,n0)
  var.t.plus <- Var.t.minus(mu0,mu.plus,s0,s.plus,n0,n.plus)

  youden.var.res <- Var.Youden(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,n.minus,n0,n.plus,t.minus,t.plus)

  var.youden <- youden.var.res$var0##get variance on youden estimate
  partialDeriv <- youden.var.res$partialDeriv##get partial derivatives and other terms needed on testing markers
  
  ####calculate Fisher's Z transformation on youden and the variance of z
  youden.z <- FisherZ(youden)
  var.youden.z <- FisherZ.Var(youden,var.youden)
  ######
    
  res0 <- CI.normal(t.minus,var.t.minus,alpha)
  t.minus.CI <- c(res0$lower,res0$upper)
  names(t.minus.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
  
  res0 <- CI.normal(t.plus,var.t.plus,alpha)
  t.plus.CI <- c(res0$lower,res0$upper)
  names(t.plus.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
  
  res0 <- CI.normal(youden,var.youden,alpha)
  youden.CI <- c(res0$lower,res0$upper)
  names(youden.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
  
  res0 <- CI.normal(youden.z,var.youden.z,alpha)
  youden.z.CI <- c(res0$lower,res0$upper)
  names(youden.z.CI) <- c(paste(alpha/2*100,"%",sep=""),paste(100-alpha/2*100,"%",sep=""))
 
  
  #out <- data.frame(t.minus=t.minus,t.plus=t.plus,youden=youden,youden.z=youden.z,var.t.minus=var.t.minus,var.t.plus=var.t.plus,var.youden=var.youden,var.youden.z=var.youden.z,t.minus.lower=t.minus.CI[1],t.minus.upper=t.minus.CI[2],t.plus.lower=t.plus.CI[1],t.plus.upper=t.plus.CI[2],youden.lower=youden.CI[1],youden.upper=youden.CI[2],youden.z.lower=youden.z.CI[1],youden.z.upper=youden.z.CI[2])
  out <- list(var.youden=var.youden,var.t.minus=var.t.minus,var.t.plus=var.t.plus,var.youden.z=var.youden.z,youden.CI=youden.CI,t.minus.CI=t.minus.CI,t.plus.CI=t.plus.CI,youden.z.CI=youden.z.CI,partialDeriv=partialDeriv)
  
  return(out)
  
}

