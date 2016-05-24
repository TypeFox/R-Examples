optim.step1.2 <-
function(dat, coef, random.eff){
  
  lam=dat$lam
  lam0=which(lam==0)
  lam1=which(lam==1)
  n.lam0=length(lam0)
  n.lam1=length(lam1)
  dfs=dat$dfs
  n.coef=length(coef)
  A=coef[1]
  B=coef[2]
  beta=coef[3:(n.coef-4)]
  ss0=coef[n.coef-3]
  ss1=coef[n.coef-2]
  rho=coef[n.coef-1]
  sigma=coef[n.coef]
  b=random.eff[, 1]
  d=random.eff[, 2]
  
  optimdat=dat$dat
  dimdat=dim(optimdat)
  nn=dimdat[1]
  ids=unique(optimdat[,1])
  mm=length(ids)
  D=as.matrix(optimdat[,4:(sum(dfs)+4)])
  Z=matrix(0, nrow=nn, ncol=mm)  
  
  for(i in 1:mm){
    id.x=(optimdat[,1]==ids[i])
    ni=sum(id.x)
    t1=rep(1, ni)
    Z[id.x,i]=t1
  }
  
  Sigma=matrix(c(ss0^2,rho*ss0*ss1,rho*ss0*ss1,ss1^2),2,2)/(sigma^2)
  
  Gi=solve(t(chol(Sigma)))  ##cholosky decomposition; t(the upper matrix)
  G=kronecker(diag(mm), Gi)
  delta=0.1   # the judgement of the stop of the while loop
  w=rep(0, 2*mm)
  
  while(delta>0.00001) {
    
    DDt=D%*%as.matrix(beta)
    f0.1=-as.matrix(A*exp(Z%*%b))/(1+exp(-log(DDt)/(B*exp(Z%*%d))))
    GG=c(seq(1, mm, 1)*2-1)
    GG2=c(seq(1, mm, 1)*2)
    w[GG]=b
    w[GG2]=d
    f0.2=G%*%w
    #w is the random effect vector with size of 72 w=(b, d)
    f0=as.matrix(rbind(f0.1, f0.2))
    
    y=as.matrix(rbind(as.matrix(optimdat[,3]), as.matrix(rep(0, 2*mm))))
    nn.new=dim(y)[1]
    
    
    
    f0.A.1=-as.matrix(exp(Z%*%b))/(1+exp(-log(DDt)/(B*exp(Z%*%d))))
    f0.A=as.matrix(rbind(f0.A.1, as.matrix(rep(0, 2*mm))))
    
    f0.B.1=-f0.1*(1+exp(-log(DDt)/(B*exp(Z%*%d))))^(-1)*exp(-log(DDt)/(B*exp(Z%*%d)))*(log(DDt)/exp(Z%*%d))*B^(-2)
    f0.B=as.matrix(rbind(f0.B.1, as.matrix(rep(0, 2*mm))))
    
    
    f0.beta.miss=f0.1*(1+exp(-log(DDt)/(B*exp(Z%*%d))))^(-1)*exp(-log(DDt)/(B*exp(Z%*%d)))*((DDt)*(B*exp(Z%*%d)))^(-1)
    applyftn1=function(x){x*f0.beta.miss}
    f0.beta.1=apply(D, 2, applyftn1)
    f0.beta.2=matrix(0, nrow=2*mm, ncol=(dim(D)[2]-n.lam0))
    
    f0.beta.0=as.matrix(f0.beta.1[,lam0])
    f0.beta.c=f0.beta.1[,lam1]
    
    f0.beta.0=rbind(f0.beta.0, matrix(0, nrow=2*mm, ncol=n.lam0))
    f0.beta.c=rbind(f0.beta.c, f0.beta.2)
    
    
    
    applyftn2=function(x){x*f0.1}
    f0.b.1=apply(Z, 2, applyftn2)
    
    f0.b.2=G[,GG]
    f0.b=rbind(f0.b.1, f0.b.2)
    
    f0.d.1=f0.1*(1+exp(-log(DDt)/(B*exp(Z%*%d))))^(-1)*exp(-log(DDt)/(B*exp(Z%*%d)))*(log(DDt)/(B*exp(Z%*%d)))
    applyftn3=function(x){x*f0.d.1}
    f0.d.1=apply(-Z, 2, applyftn3)
    
    f0.d.2=G[,GG2]
    f0.d=rbind(f0.d.1, f0.d.2)
    
    
    
    y.new=y-f0+f0.beta.c%*%beta[lam1]
    X1=cbind(f0.A, f0.B, f0.beta.0, f0.b, f0.d)
    X2=f0.beta.c
    I=diag(nn.new)
    
    beta.new=cls(y.new, (I-Px(X1))%*%X2)
    y.part.hat=X2%*%beta.new$betahat
    y.part1.hat=y.new-y.part.hat
    increment=as.numeric(solve(t(X1)%*%X1)%*%t(X1)%*%y.part1.hat)   #parameter increment
    
    A=increment[1]+A
    B=increment[2]+B
    temp=beta
    if(n.lam0!=0){
      beta[lam0]=increment[3:(2+n.lam0)]+beta[lam0]
    }    
    beta[lam1]=beta.new$betahat
    b=b+increment[(3+n.lam0):(length(increment)-mm)]
    d=d+increment[(length(increment)-mm+1):length(increment)]
    diff.beta=abs(temp-beta)
    delta=max(abs(increment), diff.beta)
   # print(c("delta=", delta))
  }
  
  y.hat=-as.matrix(A*exp(Z%*%b))/(1+exp(-log(D%*%beta)/(B*exp(Z%*%d))))
  Zhat=matrix(0, nrow=nn, ncol=2*mm)
  applyftn2=function(x){x*y.hat}
  Zhat[,GG]=apply(Z, 2, applyftn2)
  f0.d.1=y.hat*(1+exp(-log(D%*%beta)/(B*exp(Z%*%d))))^(-1)*exp(-log(D%*%beta)/(B*exp(Z%*%d)))*(log(D%*%beta)/(B*exp(Z%*%d)))
  applyftn3=function(x){x*f0.d.1}
  Zhat[,GG2]=apply(-Z, 2, applyftn3)
  
  time=optimdat[,2]
  ww=optimdat[,3]-y.hat
  Dt=round(log(D%*%beta), 8)
  colnames(Dt)="Dt"
  
  id=as.matrix(Z%*%c(1:mm))
  colnames(id)="id"
  DAMAGE_Y=as.matrix(optimdat[,3])
  res.dat=cbind(id=id, TIME=optimdat[,2], Dt=Dt , DAMAGE_Y=DAMAGE_Y)
  colnames(res.dat)=c("id", "TIME", "Dt", "DAMAGE_Y")
  newdata=as.data.frame(res.dat)
  random.eff=cbind(b, d)
  coef=c(A, B, beta, ss0, ss1, rho, sigma)
  res=list(coef=coef, random.eff=random.eff, yhat=y.hat, Zhat=Zhat,  time=time, residual=ww, dat=dat, newdata=newdata)
  return(res)
}
