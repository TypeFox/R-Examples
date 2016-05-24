predict.CGP <-
function(object,newdata=NULL,PI=FALSE,...){
  
  UU<-newdata
  DD<-object$X
  yobs<-object$yobs
  n<-nrow(DD)
  p<-ncol(DD)
  one<-rep(1,n)
  
  lambda<-object$lambda
  theta<-as.vector(object$theta)
  alpha<-as.vector(object$alpha)
  bw<-object$bandwidth
  Sig<-object$Sig_matrix
  sf<-object$sf
  res2<-object$res2
  temp<-object$temp_matrix
  invQ<-object$invQ
  tau2<-object$tau2
  beta<-object$mu

  if(is.null(UU)){
    Yp<-NULL
    gp<-NULL
    lp<-NULL
    v<-NULL
    Y_low<-NULL
    Y_up<-NULL
  }
  
  if(!is.null(UU)){
    UU<-as.matrix(UU)
    N<-nrow(UU)
    if(ncol(UU)!=p) print("Predictive location input UU is of wrong dimension!")
    ppp<-rep(0,N)
    g<-rep(0,n)
    gbw<-rep(0,n)
    l<-rep(0,n)
    Yp<-rep(0,N)
    gp<-rep(0,N)
    lp<-rep(0,N)
    v<-rep(0,N)
    for (k in 1:N){
      for (r in 1:n){
        g[r]<-exp(-(DD[r,]-UU[k,])^2%*%(theta))
        gbw[r]<-exp(-(DD[r,]-UU[k,])^2%*%(theta*bw))
        l[r]<-exp(-(DD[r,]-UU[k,])^2%*%(alpha))
      }
      v[k]<-(t(gbw)%*%(res2)/(t(gbw)%*%one))/sf
      q<-g+lambda*sqrt(v[k])*Sig^(1/2)%*%l
      Yp[k]<-beta+t(q)%*%temp
      gp[k]<-beta+t(g)%*%temp
      if(PI){
        lp[k]<-lambda*sqrt(v[k])*t(l)%*%Sig^(1/2)%*%temp
        ppp[k]<-1+lambda*v[k]-t(q)%*%invQ%*%q+(1-t(q)%*%invQ%*%one)^2/(one%*%invQ%*%one)
      }
    }
    if(PI){
      ppp[ppp<0]<-0
      ka<-1.96
      Y_up<-Yp+ka*sqrt(tau2*ppp)
      Y_low<-Yp-ka*sqrt(tau2*ppp)
    }
    if(!PI){
      Y_low<-NULL
      Y_up<-NULL
    }
  }
  
  val<-list(Yp=Yp,gp=gp,lp=lp,v=v,Y_low=Y_low,Y_up=Y_up)
  return(val)
  
}
