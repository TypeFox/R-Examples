betabreaker <-
function(object){
  frac.step=function(x){
    ##  step=x[1]-1
    ##  start=x[2]
    ##  end=x[3]
    ##  step -start/(end-start)
    x[1]-1 -x[2]/(x[3]-x[2])
  }
  beta=object$beta
  sbeta=sign(beta)
  kp=dim(beta)
  k=kp[1];p=kp[2]
  dsbeta=abs(sbeta[-1,]-sbeta[-k,])
  if(any(dsbeta==2)){
  bbeta=matrix(cbind(
    step.end=rep(1:(k-1),p),
    beta.start=as.vector(beta[-k,]),
    beta.end=as.vector(beta[-1,])
    )[dsbeta==2],ncol=3)
  fsteps=apply(bbeta,1,frac.step)
  new.beta=predict(object,type="coefficient",s=fsteps+1,mode="step")$coef
 new.beta=rbind(beta,new.beta)
  fo=c(seq(k)-1,fsteps)
  beta= new.beta[order(fo),]
  dimnames(beta)[[1]]=format(round(sort(fo),2))
}
  beta
  
}

