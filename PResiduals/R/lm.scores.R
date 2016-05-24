#' @importFrom stats lm dnorm density
lm.scores = function(y, X, emp=FALSE){
  N = length(y)  
  mod = lm(y~X)
  smod = summary(mod)
  resid = smod$residuals
  ## bread = [1/N sum (- partial phi)]^-1 
  ##       = [- 1/N  d2l. dtheta. dtheta]^-1
  #d2l.dtheta.dtheta = - solve(bread(mod))*N
  d2l.dtheta.dtheta = -crossprod(cbind(1, X))
  #dl.dtheta = estfun(mod)
  dl.dtheta <- resid*cbind(1, X)
  dresid.dtheta = t(cbind(-1, -X))

  if(emp == TRUE) {
      presid = 2*pnorm((y - mod$fitted.values)/smod$sigma) -1
      dpresid.dtheta = t(cbind(-2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma,
                               -2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma *
                               X))
  }
  else {
      f.y<-density(resid)
      fy.ry <- NULL
      presid <- NULL
      for (i in 1:length(resid)){
          fy.ry[i] <- f.y$y[which(abs(f.y$x-resid[i])==min(abs(f.y$x-resid[i])))]
          presid[i] <- sum(resid<resid[i])/length(resid) - sum(resid>resid[i])/length(resid)
      }
      dpresid.dtheta <- t(cbind(-2*fy.ry,
                                  -2*fy.ry*X))
  }
  
  list(mod = mod, 
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = resid,
       dresid.dtheta = dresid.dtheta,
       presid = presid,
       dpresid.dtheta = dpresid.dtheta)
}
