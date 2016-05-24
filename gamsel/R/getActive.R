getActive<-
  function(object, index=NULL, type=c("nonzero","linear","nonlinear"),EPS=0) {
    if(is.null(index))index=seq(along=object$interecept)
    type=match.arg(type)
    anya=abs(object$alpha[,index,drop=FALSE])>EPS
    degrees=object$degrees
    p=length(degrees)
    counter=rep(1:p,degrees)
    counter=diag(p)[counter,]
    anyb=t(counter)%*%abs(object$beta[,index,drop=FALSE]) > EPS
    anytype=switch(type,
      linear=anya,
      nonlinear=anyb,
      nonzero=anya|anyb
      )
    out=lapply(seq(along=index),function(j){which=anytype[,j]; if(any(which))seq(p)[which]else NULL})
    names(out)=paste("l",index,sep="")
    out
 }
