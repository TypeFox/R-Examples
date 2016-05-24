plot.gmm <- function(x,...){
  if(class(x)!="gmm"){stop("This object is not of class gmm")}
  ncomp<-ncol(x$mu)
  cols<-c("red","blue","green","paleturquoise","yellow","violet","orange","slategrey")
  if(ncomp>8){cols<-rep("black",ncomp)}
  
  r<-range(x$x)
  se<-seq(r[1],r[2],length.out=length(x$x))
  d<-list()
  for(i in 1:ncomp){
  d[[i]]<-dnorm(x = se,mean = x$mu[i],sd = x$sigma[i])*x$lambda[i]
  }
  
  h<-hist(x$x, plot=F)
  br<-mean(diff(h$breaks))
  h$counts <- h$counts / sum(h$counts)
  ylim<-c(0,max(c(max(unlist(d)*br),max(h$counts))))
  plot(h, freq=TRUE, ylab="Relative Frequency",xlab="Data",main="",ylim=ylim)
  for(i in 1:length(d)){
  lines(se,d[[i]]*br,col=cols[i],lwd=5)
    }
}


