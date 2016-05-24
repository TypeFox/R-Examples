ci.BTLLasso <- function(object,  subs=NULL, rows= NULL, subset = NULL){
  
model <- object$cv.model
epsilon <- model$epsilon
accuracy <- -log10(epsilon)
covariates <- colnames(model$X)
conf.ints <- object$conf.ints.repar

  m <- model$m
  labels <- model$labels
  n.theta <- model$n.theta
  estimates <- model$coefs.repar[which.min(model$deviances),]
  estimates <- round(estimates,accuracy)
  conf.ints <- round(conf.ints,accuracy)
  
  unpen <- estimates[1:(n.theta+m)]
  
  unpen.ci <- conf.ints[,1:(n.theta+m)]
  
  gamma.ci <- conf.ints[,(n.theta+m+1):ncol(conf.ints)]
  gamma <- estimates[(n.theta+m+1):ncol(conf.ints)]
  
  p <- length(gamma)/(m)
  
  if(is.null(subset)){subset <- 1:p}  
  
  
  if(!is.null(rows)){
    cols <- ceiling(length(subset)/rows)
  }else{
    cols <- floor(sqrt(length(subset)))
    
    rows <- ceiling((length(subset))/cols)
  }
  
  
    layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))

  index <- 1
  for(i in 1:(p)){
    if(i %in% subset){
    plot(gamma[index:(index+m-1)],1:m,xlim=range(gamma.ci[,index:(index+m-1)]),pch=16,yaxt="n",xlab="",ylab="",
         main=covariates[i])  
    segments(y0=1:m,x0=gamma.ci[1,index:(index+m-1)],x1=gamma.ci[2,index:(index+m-1)])
    axis(2,at =1:m, labels=labels,las=2)
    mtext(subs[i],side=3,line=0.5,cex=par()$cex)
    abline(v=0,lty=2,col="lightgray")
    }
    index <- index + m
  }

layout(1)

}

