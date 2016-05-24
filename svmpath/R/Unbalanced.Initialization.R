Unbalanced.Initialization<-function(K,y,nplus,nminus,eps=1e-8){
pos.init<-function(K,y,nplus,nminus,eps){
  
### Build up initial Left and dominant class indices, 
  Iplus<-seq(y)[y>0]
  Iminus<-seq(y)[y<0]

  Left<-Iminus
  Kscript <- K * outer(y, y)	
  Rmat<-Kscript[Iplus,Iplus]
  c.objective<-Kscript[Iplus,Iminus]%*%rep(1,nminus)
  alpha.opt=InitsvmPath(Rmat,c.objective,nminus)$alpha
### alpha.opt<-OptInit.alpha(nplus,nminus,Rmat,c.objective)$alpha# LSSOL wierd bug
###alpha.opt<-new.init(nplus,nminus,Rmat,c.objective,eps*eps)#quadprog - unsatisfactory
 
  alpha<-rep(1,length(y))
  alpha[Iplus]<-alpha.opt
  
  zeros<-alpha.opt<eps
  Right<-if(any(zeros))Iplus[zeros] else NULL
  elbows<-(!zeros)&(alpha.opt<1-eps)
  Elbow<-if(any(elbows))Iplus[elbows] else NULL
  Iplus.Left<-setdiff(Iplus,c(Right,Elbow))
  Left<-c(Left,Iplus.Left)
  
  f<-K%*%(y*alpha)

  fmin<-min(f[Iminus])
    
  if(length(Elbow)>0)
    fmax<-max(f[Elbow])
  else
    {
      fmax<- max(f[Iplus.Left])
      Elbow <-Iplus.Left[abs(f[Iplus.Left]-fmax)<eps]
      Left<-setdiff(Left,Elbow)
    }
  iminus<-Iminus[abs(f[Iminus]-fmin)<eps]
  Elbow<-c(Elbow,iminus)
  Left<-setdiff(Left,iminus)
  lambda<- (fmax-fmin)/2
  beta0<-1-fmax/lambda
  list(Elbow=Elbow,Right=Right,Left=Left,lambda=lambda,alpha=alpha,beta0=beta0)
}

###Find out which class is dominant
more.y <- +1
if(nminus>nplus){
  more.y<- -1
  ### We like to work with + class dominant
  init<-pos.init(K,-y,nminus,nplus,eps)
  ###Fix up some of the entries
  init$beta0<- -init$beta0
}
else init<- pos.init(K,y,nplus,nminus,eps)
alpha0<-init$beta0*init$lambda
 
### package parameters for the left of the start
### alpha0 is linear in lambda, so just solve a system of equations
alpha00<-c(slope = more.y, intercept=init$lambda*(init$beta0-more.y))
c(init,list(alpha0=alpha0,alpha00=alpha00))
}

      
