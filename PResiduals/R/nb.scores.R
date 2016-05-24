#' @importFrom stats dnbinom
#' @importFrom MASS glm.nb
nb.scores <- function(y, X, ...){
  N <-  length(y)

  mod <- glm.nb(y~X,  ...)

  #### score for beta's
  dl.dbeta <- (y - mod$fitted.values)/(1 + mod$fitted.values/mod$theta) * cbind(1, X)

  #### score for 1/size  (size = theta in glm.nb)
  dl.dsize <- unlist(lapply(y, FUN=function(y) ifelse(y==0, 0, sum(1/(mod$theta + 0:(y-1)))))) - log( 1 + mod$fitted.values/mod$theta) + (mod$fitted.values - y)/(mod$theta + mod$fitted.values)

  dl.dtheta <- cbind(dl.dbeta, dl.dsize)
  
  #### fisher information for beta's
  d2l.dbeta.dbeta <- crossprod(cbind(1, X)*sqrt(mod$fitted.values/(1+mod$fitted.values/mod$theta)))

  #### fisher information for theta
  d2l.dsize.dsize <- sum( unlist(lapply(y, FUN=function(y) ifelse(y==0, 0, sum(1/(mod$theta + 0:(y-1))^2)))) - mod$fitted.values/mod$theta/(mod$fitted.values + mod$theta))
  d2l.dtheta.dtheta <- rbind(cbind(d2l.dbeta.dbeta,0),0)
  d2l.dtheta.dtheta[dim(d2l.dtheta.dtheta)[1] , dim(d2l.dtheta.dtheta)[2] ] <- d2l.dsize.dsize
  presid <- pnbinom(y-1, mu=mod$fitted.values, size=mod$theta ) + pnbinom(y, mu=mod$fitted.values, size=mod$theta) -1
  pearson.resid <- (y - mod$fitted.values)/sqrt(mod$fitted.values + mod$fitted.values^2/mod$theta)
  dmu.dbeta <- mod$fitted.values * cbind(1, X)
  dpresid.dmu <- dnbinom(y-1, size=mod$theta, mu=mod$fitted.values)/(mod$fitted.value + mod$theta) - (dnbinom(y-1, size=mod$theta, mu=mod$fitted.values) + dnbinom(y, size=mod$theta, mu=mod$fitted.values))*(mod$theta + y)/(mod$theta + mod$fitted.values)
  dpresid.dbeta <- t(dpresid.dmu * dmu.dbeta)
  df.dsize <- function(k, mu, size) {
    dnbinom(k, size=size, mu=mu)*((mu-k)/(mu+size) +
                                     log(size) -log(mu+size) +
                                     unlist(lapply(k, FUN=function(k) ifelse(k==0, 0, sum(1/(size + (0:(k-1))))))) )                                  
  }
  
  ### dpresid.dsize = 2 dF(y-1).dsize + df(y).dsize = 2 sum(df(y-1).dsize) + df(y).dsize
  #tmp <-  cbind(y, mod$fitted.values)
  #tmp2 <- apply(tmp, 1, FUN=function(x) ifelse(x[1]==0, 0, 2*sum(df.dsize(k=0:(x[1]-1), mu=x[2], size=mod$theta))))
 
  dpresid.dsize <- df.dsize(y, mod$fitted, mod$theta)  +
    apply(cbind(y, mod$fitted.values), 1, FUN=function(x) ifelse(x[1]==0, 0, 2*sum(df.dsize(k=0:(x[1]-1), mu=x[2], size=mod$theta))))
  dpresid.dtheta <- rbind(dpresid.dbeta, dpresid.dsize)
  
  dpearson.resid.dmu <- -0.5 * (y * (mod$fitted.values + mod$fitted.values^2/mod$theta)^(-3/2)*(1+2*mod$fitted.values/mod$theta) +
                                  (1/mod$fitted.values + 1/mod$theta)^(-3/2)/mod$fitted.values^2)
  dpearson.resid.dbeta <- t(dpearson.resid.dmu * dmu.dbeta)
  dpearson.resid.dsize <- 0.5*(y-mod$fitted.values)*(mod$fitted.values + mod$fitted.values^2/mod$theta)^(-3/2)*mod$fitted.values^2/mod$theta^2
  dpearson.resid.dtheta <- rbind(dpearson.resid.dbeta, dpearson.resid.dsize)
  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       presid = presid,
       dpresid.dtheta = dpresid.dtheta,
       pearson.resid = pearson.resid,
       dpearson.resid.dtheta = dpearson.resid.dtheta)
  

}
