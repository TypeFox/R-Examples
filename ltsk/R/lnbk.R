require(sp)
require(gstat)
lnbk <-
function(query,obs,model,xcoord='x',ycoord='y',wcoord="w",bcoord='b',zcoord='z',verbose=T,epsilon=1e-5)
{
  ## query: query location data frame with associated block IDs;
  ## obs : observed location and data frame with values (assumed on point support);
  ## model: a variogram model for kriging
  ## xcoord: coordinate for easting /longitude
  ## ycoord: coordinate for northing / latitude
  ## bcoord: coordinate for block in query data frame
  ## wcoord: coordinate for query location weight (e.g. area)
  ## zcoord: coordinate for observation
  ## verbose: debug information
  ## epsilon: small constant to reduce small number inflation
  
  ## Values:
  ## for each block a estimated average with the mean square prediction error
  
  # potentially check for coordinates, for positive observed values, something else
  l.query <- as.matrix(query[,c(xcoord,ycoord,bcoord,wcoord)])
  l.obs <- as.matrix(obs[,c(xcoord,ycoord,zcoord)])
  chk1 <- all(l.obs[,3]>=0)
  if(!chk1){
    stop(zcoord," contains negatives...")
  }
  chk2 <- mean(l.obs[,3]==0)
  if(chk2>0){
    cat(round(chk2,2)," data values are zero, added small constant\n")
    l.obs[,3] <- l.obs[,3]+epsilon
  }
  l.obs[,3]  <- log(l.obs[,3]) ## transform to Gaussian scale
  ## detect zero distance
  tmp <- data.frame(l.obs[,1:2])
  coordinates(tmp) <- as.formula(paste("~",xcoord,"+",ycoord,sep=""))
  ii  <- zerodist(tmp)
  if(nrow(ii)>0){
    cat("zero distance detected\n")
    l.obs[unique(ii[,1]),1] <- l.obs[unique(ii[,1]),1]+1000*runif(length(unique(ii[,1])))
    l.obs[unique(ii[,2]),2] <- l.obs[unique(ii[,2]),2]+1000*runif(length(unique(ii[,2])))
  }
  
  s <- l.obs[,1:2]
  Sigma.Y <- variogramLine(model, dist_vector = spDists(s), covariance = TRUE) ## covariance among observed process
  l.Sigma.Y <- chol(Sigma.Y) ## cholesky factor Sigma.Y
  
  tmp <- rep(1,nrow(l.obs))
  mu.Y.xx <- sum(backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,tmp,upper.tri=T, transpose=T))) ## 1^T Sigma.y^(-1) 1
  mu.Y.xy <-  sum(backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,l.obs[,3],upper.tri=T, transpose=T))) ## 1^T Sigma.y^(-1) Y
  mu.Y <- mu.Y.xy / mu.Y.xx  ## GLS esimtate of the unknown mean
  
  s0 <- l.query[,1:2]
  sigma.Ys0 <- variogramLine(model, dist_vector = spDists(s0,s), covariance = TRUE) ## covariance among observed process
  res.xy <- backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,l.obs[,3]-mu.Y,upper.tri=T, transpose=T)) ## Sigma.y^(-1)(Y-mu.Y)
  res.xy <- as.vector(sigma.Ys0 %*% res.xy)  ## sigma.ys0^T %*% Sigma.y^(-1)(Y-mu.Y)
  
  yhat <- mu.Y + res.xy ## ordinary Kriging estimate

  m.num <- backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,t(sigma.Ys0),upper.tri=T, transpose=T)) ## Sigma.y^(-1) %*% {C_j^T}
  m.num.sum <- colSums(m.num) ## 1^T %*% Sigma.y^(-1) %*% {C_j^T}
  m.num.out <- 1 - m.num.sum
  m <- m.num.out / mu.Y.xx ## Lagrange multiplier estimate
  
#   tmp <- rep(0,nrow(l.query))
#   for(i in 1:length(tmp)){
#     c1 <- sigma.Ys0[i,]
#     c2 <- solve(Sigma.Y,c1)
#     tmp[i] <- sum(c2)
#   }
  #browser()
  c.adj <- outer(rep(1,nrow(l.obs)),m,"*")
  c.adjusted <- t(sigma.Ys0) + c.adj
  lambda <-  backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,c.adjusted,upper.tri=T, transpose=T)) ## kriging weights nxm
  
#   yhat2 <- as.vector(crossprod(lambda,l.obs[,3]))
#   summary(yhat-yhat2)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -1.954e-14 -7.105e-15  0.000e+00 -2.006e-15  2.665e-15  1.421e-14 

  totvar <- variogramLine(model, dist_vector = rep(0,nrow(l.query)), covariance = TRUE) ## C_Y(s_0,s_0)
  var.part2 <- rowSums(t(lambda) * sigma.Ys0) ## lambda(s0)^T cY(s0)
  sigma.s0 <-  totvar$gamma - var.part2 + m
  
  zhat <- exp(yhat + 0.5*sigma.s0 -m) ## kriging estimates
  
  zhat.w <- zhat*l.query[,4] ## weighted estimates
  zB.w <- tapply(zhat.w,l.query[,3],sum) ## aggregate by block
  zB.d <- tapply(l.query[,4],l.query[,3],sum) ## total weight by block
  zBhat <- zB.w / zB.d ## weighted average by block
  #tmp <- tapply(zhat,l.query[,3],mean)
  #summary(zBhat-tmp)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -4.441e-16  0.000e+00  0.000e+00 -1.542e-17  0.000e+00  4.441e-16 
  #browser()
  lst.sigma.Ys0 <- rowSplit(sigma.Ys0,l.query[,3]) ## list of cY(u) for each Block
  lst.m <- split(m,l.query[,3]) ## list of m(u) Lagrange multiplier for each Block
  lst.loc <- rowSplit(l.query[,1:2],l.query[,3]) ## list of coordinates x(u), y(u) for each Block
  lst.w <- split(l.query[,4],l.query[,3]) ## list of areal weights
  
  fac <- exp(2*mu.Y+totvar$gamma[1]) ## exp{mu.Y+1/2 C_Y(u,u)} * exp{mu.Y+1/2 C_Y(v,v)}
  do.mspe <- function(loc,w,sigma.Ys0,m){
    ## loc: coordinates
    ## w : areal weights
    ## sigma.Ys0: covariance with observed data
    ## m : Lagrange multiplier
    ## return mean squared prediction error
    #browser()
    c.adj <- outer(m,rep(1,ncol(sigma.Ys0)),"*") # fill Lagrange multiplier by each data location
    c.adjusted <- sigma.Ys0 + c.adj # add Lagrange multiplier c'Y(u)= cY(u)+1 m(u)
    a.mat <- variogramLine(model, dist_vector = spDists(loc), covariance = TRUE) ## a=exp{C_Y(u,v)}
    a.mat <- exp(a.mat)
    b.sigmainv <- backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,t(sigma.Ys0),upper.tri=T, transpose=T)) ## Sigma_y^(-1) cY(v)
    b.mat <- exp(c.adjusted %*% b.sigmainv) ## exp{c'Y(u)^T Sigma_y^(-1) cY(v)}
    c.mat <- t(b.mat) ## exp{c'Y(v)^T Sigma_y^(-1) cY(u)}
    d.sigmainv <- backsolve(l.Sigma.Y,forwardsolve(l.Sigma.Y,t(c.adjusted),upper.tri=T, transpose=T)) ## exp{c'Y(u)^T Sigma_y^(-1) cY(v)}
    d.mat <- exp(c.adjusted %*% d.sigmainv)
    w.mat <- outer(w,w,"*")
    
    tmp <- fac*(a.mat-b.mat-c.mat+d.mat) ## integrand in (24) Cressie 2006
    sum(tmp * w.mat) / sum(w.mat)  ## weighted average
  }
  mspe <- mapply(do.mspe,lst.loc,lst.w,lst.sigma.Ys0,lst.m)
  #browser()
  r <- cbind(sort(unique(l.query[,3])),tapply(l.query[,4],l.query[,3],sum),zBhat,mspe,sqrt(mspe))
  colnames(r) <- c(bcoord,wcoord,"krig","mspe","rmspe")
  return(r)
}
