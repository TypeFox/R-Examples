histStack<-function(x,...) {
 UseMethod("histStack")
}

histStack.formula<-function(x,data,breaks="Sturges",col="rainbow",right=TRUE,
 main="",xlab=NULL,legend.pos=NULL,cex.legend=0.75,...) {

 mf<-model.frame(x,data=data)
 if(is.null(xlab)) xlab<-names(mf)[1]
 histStack.default(mf[,1],mf[,2],breaks=breaks,col=col,right=right,
  main=main,xlab=xlab,legend.pos=legend.pos,cex.legend=cex.legend,...)
}

histStack.default<-function(x,z,breaks="Sturges",col="rainbow",right=TRUE,
 main="",xlab=NULL,legend.pos=NULL,cex.legend=0.75,...) {

 if(!is.factor(z)) {
  z<-factor(z)
  warning("z was converted to a factor")
 }
 seps=levels(z)
 numseps<-length(seps)
 if(length(col) == 1) col<-do.call(col,list(n=numseps))
 if(length(col) < numseps) col<-rep(col,length.out=numseps)
 if(!is.numeric(x)) stop("x must be numeric",call.=FALSE)
 # plot the histogram of all x
 hS<-hist(x,breaks=breaks,col=col[1],right=right,main=main,xlab=xlab,...)
 # plot the remaining 
 for(i in 1:(numseps-1))
  hist(x[z %in% seps[-(1:i)]],breaks=hS$breaks,col=col[i+1],
  right=right,add=TRUE)
 if(!is.null(legend.pos)) {
  if(length(legend.pos > 1))
   legend(legend.pos[1],legend.pos[2],seps,fill=col,cex=cex.legend)
  else legend(legend.pos,seps,fill=col,cex=cex.legend)
 }
}
