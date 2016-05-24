singlepaths <- function(model, colors = NULL, equal.ranges = FALSE, subs = NULL){

  coefs <- model$coefs
  covar <- colnames(model$X)
  
  acoefs <- model$acoefs
  norm <- rowSums(abs(coefs%*%acoefs))
  norm <- norm/max(norm)
  
  m <- model$m
  n.theta <- model$n.theta
  
  labels <- model$labels
  
  coefs <- model$coefs.repar

  if(n.theta>0){
    theta <- coefs[,1:n.theta,drop=FALSE]
  }
  
  intercepts <- coefs[,(n.theta+1):(n.theta+m)]
  
  gamma <- coefs[,(n.theta+m+1):ncol(coefs)]
  
  p <- ncol(gamma)/(m)
  
  cols <- floor(sqrt(p))
  
  rows <- ceiling((p)/cols)
  
  layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))
  
  y.range <- range(gamma)
  
  if(!is.null(model$deviances)){
    x.axis.min <-norm[which.min(model$deviances)]
  }
  
  if(is.null(colors)){colors <- rep(1,m)}
  
  index <- 1
  for(i in 1:p){
    if(!equal.ranges){y.range <- range(gamma[,index:(index+m-1)])}
    plot(norm, gamma[,index],ylim=y.range,type="l",main=covar[i],ylab="",
         xlab=expression(sum(sum(abs(beta["rj"]-beta["sj"]),"r<s"),"j")/max(sum(sum(abs(beta["rj"]-beta["sj"]),"r<s"),"j"))),
         col=colors[1])

    for(u in 1:(m-1)){
      lines(norm,gamma[,index+u],col=colors[u+1])
    }
    axis(4,at=gamma[nrow(gamma),index:(index+m-1)],labels = labels, las=2)

    mtext(subs[i],side=3,line=0.5,cex=par()$cex)
    if(!is.null(model$deviances)){
    abline(v=x.axis.min,lty=2,col=2)
    }
    index <- index+m
  }
  
  layout(1)
}