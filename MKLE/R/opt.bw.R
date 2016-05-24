"opt.bw" <- function(data, bws=c(sd(data),4*sd(data)), B=1000, gridsize=2^14){
  "mse"<-function(bw,data,B){
    n<-length(data)
    mkles<-apply(matrix(sample(data,n*B,replace=TRUE),ncol=B),2,mkle,bw=bw,gridsize=gridsize)
    (mean(mkles)-mean(data))^2+var(mkles)
  }

  res<-optimize(mse,interval=bws,data=data,B=B)

  mkles<-apply(matrix(sample(data,length(data)*B,replace=TRUE),ncol=B),2,mean)
  mse0<-(mean(mkles)-mean(data))^2+var(mkles)

  if(mse0<res$objective)return(0)

  return(res$minimum)

}
