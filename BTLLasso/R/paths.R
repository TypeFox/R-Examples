paths <- function(model){
  
  coef <- model$coefs
  covar <- colnames(model$X)
  
  acoefs <- model$acoefs
  norm <- rowSums(abs(coef%*%acoefs))
  norm <- norm/max(norm)
  
  m <- model$m
  n.theta <- model$n.theta
  
  deviances <- model$deviances
  
 coef <- model$coefs.repar


  if(n.theta>0){
  theta <- coef[,1:n.theta,drop=FALSE]
  }
  
  intercepts <- coef[,(n.theta+1):(n.theta+m)]
  
  gamma <- coef[,(n.theta+m+1):ncol(coef)]
  
  p <- ncol(gamma)/(m)

  if(!is.null(deviances)){
    x.axis.min <- norm[which.min(deviances)]
  }

pos <- combn(m,2)
diff_mat <- matrix(0,nrow=m,ncol=choose(m,2))
for(o in 1:ncol(pos)){
  diff_mat[pos[1,o],o] <- 1
  diff_mat[pos[2,o],o] <- -1
}

cov_norms <- matrix(0,nrow=nrow(coef),ncol=p)
index = 1
for(o in 1:p){
  gamma_p <- gamma[,index:(index+m-1)]
  cov_norms[,o] <- rowSums(abs(gamma_p%*%diff_mat))
  index <- index + m
}


plot(norm,cov_norms[,1],type="l",ylim=range(cov_norms),ylab="",
     xlab=expression(sum(sum(abs(beta["rj"]-beta["sj"]),"r<s"),"j")
                     /max(sum(sum(abs(beta["rj"]-beta["sj"]),"r<s"),"j"))),
     las=1)
for(o in 2:p){
  lines(norm,cov_norms[,o])
}
if(!is.null(deviances)){
abline(v=x.axis.min,lty=2,col=2)
}  

axis(4,at=cov_norms[nrow(cov_norms),],labels=covar,las=2)



}