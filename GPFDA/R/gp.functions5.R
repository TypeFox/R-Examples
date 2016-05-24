######################## Covariance functions ############################
cov.linear=function(hyper,Data,Data.new=NULL){
  #data should have the form that, each column is a variable
  #k(x^1,x^2)=sum_{q=1}^{Q} a_q * x^1_q * x^2_q
  #hyper is a list of hyper-parameters
  if(is.null(Data.new)) Data.new=Data
  hyper=lapply(hyper,exp);n.hyper=length(hyper$linear.a)
  cov.lin=xixj(Data.new,Data,hyper$linear.a)
  return(cov.lin)
}

cov.pow.ex=function(hyper,Data,Data.new=NULL,gamma=1){
  #hyper is a list of hyper-parameters
  if(is.null(gamma)) gamma=1
  hyper=lapply(hyper,exp);
  datadim=dim(Data)
  #Data=t(t(Data)*hyper$w^(1/hyper$gamma))
  v.power=xixj_sta(Data,Data.new,hyper$pow.ex.w,power=gamma)
  
  exp.v.power=hyper$pow.ex.v*exp(-v.power/2)
}

cov.rat.qu=function(hyper,Data,Data.new=NULL){
  #hyper is a list of hyper-parameters
  hyper=lapply(hyper,exp);
  datadim=dim(Data)
  
  v.power=xixj_sta(Data,Data.new,hyper$rat.qu.w,power=1)
  
  rat.qu.v=(1+log(hyper$rat.qu.s)*v.power)^(-hyper$rat.qu.a)
  return(rat.qu.v)
}



# 
# cov.matern=function(hyper,Data,Data.new=NA){
#   #hyper is a list of hyper-parameters
#   n.sample=dim(Data)[1];n.var=dim(Data)[2]
#   hyper=lapply(hyper,exp);
#   if(!is.matrix(Data.new)){
#     cov.=sapply(1:n.sample,function(i) matrix(rep(Data[i,],n.sample),nrow=n.sample,byrow=T)-Data)#v=x_q-x_q' (for each variable)
#     cov.=abs(cov.)^2*matrix(rep(hyper$matern.w,each=n.sample*n.sample),ncol=n.sample,byrow=T)#w*v^gamma (for each vraiable)
#     cov..=matrix(0,ncol=n.sample,nrow=n.sample)
#     if(n.var>1)
#       for(i in 1:(n.var-1)){
#         cov..=cov..+cov.[1:n.sample,];cov.=cov.[-(1:n.sample),]}
#     cov.=cov..+cov. # sum all variable up 
#     cov.matern=cov.}
#   
#   
#   if(is.matrix(Data.new)){
#     nn.sample=dim(Data.new)[1]
#     cov.=sapply(1:n.sample,function(i) matrix(rep(Data[i,],nn.sample),nrow=nn.sample,byrow=T)-Data.new)
#     cov.=abs(cov.)^2*matrix(rep(hyper$matern.w,each=nn.sample*n.sample),ncol=n.sample,byrow=T)
#     cov..=matrix(0,ncol=n.sample,nrow=nn.sample)
#     if(n.var>1)
#       for(i in 1:(n.var-1)){
#         cov..=cov..+cov.[1:nn.sample,];cov.=cov.[-(1:nn.sample),]}
#     cov.=cov..+cov.
#     cov.matern=t(cov.matern)
#   }
#   ## modefied bessel function problem
#   
#   return(cov.matern)
# }

######################## Predict ############################
gppredict=function(train=NULL,Data.new=NULL,hyper=NULL, Data=NULL, Y=NULL, Cov=NULL,gamma=NULL,lrm=NULL,mean=0){
  Data.new=as.matrix(Data.new)
  if(mean==1){
    mean=mean(Y)
  }
  if(!is.null(Data)) Data=as.matrix(Data)
  if(!is.null(Y)) Y=as.matrix(Y-mean)
  if(class(train)=='gpr'){
    hyper=train$hyper
    Data=train$train.x
    Y=train$train.y
    Cov=train$CovFun
    gamma=train$gamma
    mean=train$mean
    lrm=train$lrm
  }
  if(is.null(train)){
    train=gpr(Data=Data,response=Y,Cov=Cov,hyper=hyper,gamma=gamma)
  }
  if(is.null(Data.new)) Data.new=Data
  n=dim(Data)[1];nn=dim(Data.new)[1];n.var=dim(Data)[2]
  
  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  
  CovL1=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(hyper,Data,Data.new,gamma=gamma))
    else
      return(f(hyper,Data,Data.new))
  }  )
  Q1=Reduce('+',CovL1)
  
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(hyper,Data,gamma=gamma))
    if(j!='cov.pow.ex')
      return(f(hyper,Data))
  }  )
  Q=Reduce('+',CovL)
  Q=Q+diag(exp(hyper$vv),dim(Q)[1])
  
  for(i in 1:n) CovList[i]=list(paste0('diag.',Cov[i]))
  CovLn=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(hyper,Data.new,gamma=gamma))
    if(j!='cov.pow.ex')
      return(as.matrix(f(hyper,Data.new)))  }  )
  Qstar=Reduce('+',CovLn)

  # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
  # invQ=pseudoinverse(Q)
  invQ=mymatrix2(Q)$res
  QQ1=invQ%*%t(Q1)
  QR=invQ%*%Y
  
  # if(!exists('lrm')) lrm=NULL

  if(class(lrm)=='lm'){
    newtrend=data.frame(xxx=Data.new[,1])
    mean=predict(lrm,newdata=newtrend)
  } 
  mu=Q1%*%QR+mean
  sigma2=(Qstar-as.matrix(diag(Q1%*%invQ%*%t(Q1))))+exp(hyper$vv)
#   sigma2=abs(Qstar-as.matrix(colSums(t(Q1)*QQ1)))+exp(hyper$vv)
  #sigma=sqrt(sigma2)
  result=c(list('pred.mean'=mu[,1],'pred.sd'=sqrt(sigma2[,1]),'newdata'=Data.new),unclass(train),nsigma=any(Qstar<as.matrix(colSums(t(Q1)*QQ1))))
  class(result)='gpr'
  return(result)
}

gpr=function(Data, response, Cov=c('linear','pow.ex'), hyper=NULL, NewHyper=NULL, mean=0, gamma=1,itermax=100,reltol=8e-10,trace=0){#,Xprior,Xprior2){
#   set.seed(60);
  Data=as.matrix(Data)
  y.original=response
  response=as.matrix(response)
  if(is.null(hyper)){
    hyper=list()
    if(any(Cov=='linear'))
    #  hyper$linear.a=rep(log(0.01),dim(Data)[2])
      hyper$linear.a=rnorm(dim(Data)[2],sd=0.01)
    if(any(Cov=='pow.ex')){
    #  hyper$pow.ex.v=log(1)
    #  hyper$pow.ex.w=rep(log(10),dim(Data)[2])
      hyper$pow.ex.v=rnorm(1,sd=0.01)
      hyper$pow.ex.w=-(abs(rnorm(dim(Data)[2],sd=0.01)))
    }
    if(any(Cov=='rat.qu')){
    #  hyper$rat.qu.w=rep(log(1),dim(Data)[2])
    #  hyper$rat.qu.s=log(1.5)
    #  hyper$rat.qu.a=log(1.5)
      hyper$rat.qu.w=abs(rnorm(dim(Data)[2],sd=0.01))
      hyper$rat.qu.s=runif(1,0.01,0.5)
      hyper$rat.qu.a=runif(1,0.01,0.5)
    }
    hyper$vv=sample(x=c(0.2,0.5,1,1.5),1)
    hyper.nam=names(hyper)
    
    if(!is.null(NewHyper)){
      hyper.nam=c(hyper.nam,NewHyper)
      nh.length=length(NewHyper)
      for(i in 1:nh.length){
        hyper=c(hyper,runif(1,-1,1))
      }
      names(hyper)=hyper.nam
    }
  }
  if(!is.null(hyper)){
    hyper=hyper[substr(names(hyper),1,6)%in%c(Cov,'vv')]
  }  
  hp.name=names(unlist(hyper))
  
  if(mean==0) {response=response; mean=0;lrm=0}
  if(mean==1) {mean=mean(response);response=as.matrix(response-mean);lrm=1}
  if(mean=='t') {
    trend=data.frame(yyy=response,xxx=Data[,1])
    lrm=lm(yyy~xxx,data=trend); 
    response=as.matrix(resid(lrm));
    mean=fitted(lrm)
  }
  
  trace=round(trace)
  if(trace>0)
    cat(c('\n','title: -likelihood:',hp.name,'\n'),sep='     ')
  CG0 <- nlminb(unlist(hyper), gp.loglikelihood2, gp.Dlikelihood2,Data=Data,response=response,Cov=Cov,gamma=gamma,control=list(iter.max=itermax,rel.tol=reltol,trace=trace))
  # CG0 <- optim(unlist(hyper), gp.loglikelihood2, gp.Dlikelihood2,Data=Data,response=response,Cov=Cov,gamma=gamma,method='CG',control=list(maxiter=itermax,reltol=reltol,trace=trace))
  # if(trace!=F&CG0$convergence==0)
  #   cat('\n','    optimization finished. Converged.','\n')
  # if(trace!=F&CG0$convergence==1)
  #   cat('\n','    optimization finished. Failed Converge.','\n')
  if(trace>0)
    cat('\n','    optimization finished.','\n')
  CG=CG0[[1]]
  names(CG)=hp.name
  CG.df=data.frame(CG=CG,CG.N=substr(hp.name,1,8))
  names(CG.df)=c('CG','CG.N')
  hyper.cg=split(CG.df$CG,CG.df$CG.N)

  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(hyper.cg,Data,Data,gamma=gamma))
    if(j!='cov.pow.ex')
      return(f(hyper.cg,Data,Data))
  }  )
  if(length(CovL)==1)
    Q=CovL[[1]]
  if(length(CovL)>1)
    Q=Reduce('+',CovL)

  response=as.matrix(response)
  Q=Q+diag(exp(hyper.cg$vv),dim(Q)[1])
  # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
  # invQ=pseudoinverse(Q)
  invQ=mymatrix2(Q)$res
  QR=invQ%*%response
  AlphaQ=QR%*%t(QR)-invQ
  
  D2fx=lapply(seq_along(hyper.cg),function(i){
    Dp=hyper.cg[i]
    name.Dp=names(Dp)
    f=get(paste0('D2',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      D2para=f(hyper.cg,Data,gamma=gamma,inv.Q=invQ,Alpha.Q=AlphaQ)
    if(!name.Dp%in%c('pow.ex.w','pow.ex.v'))
      D2para=f(hyper.cg,Data,inv.Q=invQ,Alpha.Q=AlphaQ)
    return(D2para)
  })
  names(D2fx)=names(hyper.cg)
  II=(-1/(unlist(D2fx)*dim(Data)[1]))
  
  fitted=(Q-diag(exp(hyper.cg$vv),dim(Q)[1]))%*%invQ%*%(response)+mean
  fitted.var=exp(hyper.cg$vv)*rowSums((Q-diag(exp(hyper.cg$vv),dim(Q)[1]))*t(invQ))
  result=list('hyper'=hyper.cg,'I'=II,'fitted.mean'=fitted[,1],fitted.sd=sqrt(fitted.var),'train.x'=Data,'train.y'=response,'train.yOri'=y.original, 'CovFun'=Cov,'gamma'=gamma,'Q'=Q,'inv'=invQ,'mean'=mean,'lrm'=lrm,conv=CG0$convergence,'hyper0'=hyper)
  class(result)='gpr'
  return(result)
}                                    


########################### likelihood ######################################
gp.loglikelihood2=function(hyper.p,Data, response,Cov,gamma=1){
  #this function doesn't return anything, it's for the conjugate gradian
  #hyper is a list of hyper-parameters
  #Data should have the form that, each column is a variable
  #response is the given response vector
  #Cov is a function contains all the covariance matrix, defult is:
  ###cov.linear(hyper,Data)+cov.pow.ex(hyper,Data), but 
  #####it could be other forms.
  
  Data=as.matrix(Data)
  datadim=dim(Data)

  hp.class=substr(names(hyper.p),1,8)
  kernel.class=unique(substr(names(hyper.p),1,6))
  hp.class=data.frame(class=hp.class,hp=hyper.p)
  names(hp.class)=c('class','hp')
  hp.list=split(hp.class$hp,hp.class$class)
  hyper.p=hp.list

  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    if(j=='cov.pow.ex')
      return(f(hyper.p,Data,Data,gamma=gamma))
    if(j!='cov.pow.ex')
      return(f(hyper.p,Data,Data))
  }  )
  Q=Reduce('+',CovL)
  Q=Q+diag(exp(hyper.p$vv),dim(Q)[1])

  response=as.matrix(response)
  # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
  # invQ=pseudoinverse(Q)
  invQ=mymatrix2(Q)$res
  invQ.response=invQ%*%response
  logdetQ=sum(determinant(Q,logarithm=T)$modulus)
  
  # fX=c(0.5*logdetQ + 0.5*t(response)%*%invQ%*%response + 0.5*dim(Data)[1]*log(2*pi))
  fX=0.5*logdetQ + 0.5*t(response)%*%invQ.response + 0.5*nrow(Data)*log(2*pi)
  
  # temp=0
  # if(any(is.na(Xprior2))==F){
  # if(any(Cov=='linear')){
  #   for (d in 1:n.hyper){
  #     temp=temp+hyper$linear.a[d]+0.5*((hyper$linear.a[d]-Xprior2$mua[d])/Xprior2$sigma[d])^2}}
  #     # updating (mu, sigma, log(a))
  # if(any(Cov=='pow.ex')){
  # 	for (d in 1:n.hyper){
  # 		temp=temp+(Xprior2$linear.alpha[d]+1)*hyper$w[d]+Xprior2$mu[d]/(Xprior2$linear.alpha[d]*exp(hyper$w[d]))}
  # 		# updating (alpha, mu, log(w))
  #     temp=temp+hyper$v1[1]+0.5*((hyper$v1[1]-Xprior2$muv1[1])/Xprior2$sigmav1[1])^2
  # 	# updating (muv1[1],sigmav1[1], v1[1])
  #     temp=temp+hyper$v0+0.5*((hyper$v0-Xprior2$muv0[1])/Xprior2$sigmav0[1])^2
  # 	# updating (muv0[1],sigmav0[1], v0[1])
  # }}
  # 
  # fX=fX+temp
  return(fX)
}

gp.Dlikelihood2=function(hyper.p,  Data, response,Cov,gamma){
  #this function doesn't return anything, it's for the conjugate gradian
  #hyper is a list of hyper-parameters
  #Data should have the form that, each column is a variable
  #response is the given response vector
  #Cov is a function contains all the covariance matrix, defult is:
  ###cov.linear(hyper,Data)+cov.pow.ex(hyper,Data), but 
  #####it could be other forms.
  
  Data=as.matrix(Data)
  datadim=dim(Data);
  
  hp.class=substr(names(hyper.p),1,8)
  kernel.class=unique(substr(names(hyper.p),1,6))
  hp.class=data.frame(class=hp.class,hp=hyper.p)
  names(hp.class)=c('class','hp')
  hyper.p=split(hp.class$hp,hp.class$class)
  
  n=length(Cov)
  CovList=vector('list',n)
  for(i in 1:n) CovList[i]=list(paste0('cov.',Cov[i]))
  CovL=lapply(CovList,function(j){
    f=get(j)
    f(hyper.p,Data,Data)
  }  )
  Q=Reduce('+',CovL)
  Q=Q+diag(exp(hyper.p$vv),dim(Q)[1])

  response=as.matrix(response)
  # invQ=pseudoinverse(Q+diag(1e-9,ncol=ncol(Q),nrow=nrow(Q)))
  # invQ=pseudoinverse(Q)
  invQ=mymatrix2(Q)$res
  Alpha=invQ%*%response
  #   Alpha2=t(Alpha)%*%Alpha
  AlphaQ=Alpha%*%t(Alpha)-invQ

  Dfx=lapply(seq_along(hyper.p),function(i){
    Dp=hyper.p[i];
    name.Dp=names(Dp)
    f=get(paste0('DCov.',name.Dp))
    if(name.Dp%in%c('pow.ex.w','pow.ex.v') )
      Dpara=f(hyper.p,Data,AlphaQ,gamma=gamma)
    if(name.Dp=='vv')
      Dpara=f(hyper.p,Alpha,invQ)
    if(!name.Dp%in%c('pow.ex.w','pow.ex.v','vv'))
      Dpara=f(hyper.p,Data,AlphaQ)
    return(Dpara)
  })
  
  names(Dfx)=names(hyper.p)
  Dfx=-0.5*unlist(Dfx)
  Dfx
}






################################# tools ###########################################
rmse=function(t,a){ 
#compute the root mean squar error between two vectors
y = sqrt(sum((a-t)^2)/length(t))
return(y)
}

mymatrix=function(matrix,log=T){
#singular decomposition 
	m=matrix
	a=svd(m)
	U=a$u; V=a$v; D=a$d
	l=length(D)
	idx=which(D<2e-13)
	if(length(idx)>0)
		D=D+1e-12
	else D=D
	inv=V%*%diag(1/D)%*%t(U)
	if(log==T)
		det=log(prod(D))
	else
		det=prod(D)
	return(list("inv"=inv,"det"=det))}


mymatrix2=function(smatrix,sB='sB',det=F,log=T,jitter=1e-10){
  mat=smatrix+diag(jitter,dim(smatrix)[1])
  smatrix=as.spam(mat,eps=1e-8)
  if(is.character(sB)) sB=diag(1,dim(mat)[1])
  else sB=as.matrix(sB)
  sB=as.spam(sB,eps=1e-8)
  x=solve.spam(smatrix,sB)
  d=NULL
  if(det==T){
    L=chol(smatrix)
    if(log==T)
      d=2*log(prod(diag(L)))
    if(log==F)
      d=prod(diag(L))
  }
  return(list('res'=as.matrix(x),'det'=d))
}

xixj=function(mat,mat.new=NULL,a=NULL){
  mat=as.matrix(mat)
  mdim=dim(mat)
  #   err=1
  
  if(is.null(mat.new)){
    #     err=0
    mat.new=mat
  }
  
  if(is.null(a))  a=rep(1,mdim[2])
  if(length(a)<mdim[2]) {
    a1=rep(1,mdim[2])
    a1[1:length(a)]=a
    a=a1;rm(a1)
    warning('number of "a" is less than the number of columns, use 1 as the missing "a"')
  }
  if(length(a)>mdim[2]) {
    a=a[1:mdim[2]]
    warning('number of "a" is more than the number of columns, omit the extra "a"')
  }
  
  aa=matrix(rep(a,mdim[1]),ncol=mdim[2],byrow=T)
  out=(aa*mat)%*%t(mat.new)
  return(out)
}

xixj_sta=function(mat,mat.new=NULL,w=NULL,power=NULL){
  mat=as.matrix(mat)
  if(is.null(mat.new)) mat.new=mat
  mdim=dim(mat);mdim.new=dim(mat.new)
  cov.=matrix(sapply(1:mdim[1],function(i) matrix(rep(mat[i,],mdim.new[1]),nrow=mdim.new[1],byrow=T)-mat.new),ncol=mdim[1])
  if(is.null(power)) power=1
  cov.=((cov.)^2)^power;
  if(is.null(w)) {
    w=rep(1,mdim[2])
    warning('missing "weight", use 1 instead')
  }
  if(length(w)==1&mdim[2]>1){
    w=rep(w,mdim[2])
    warning('only one "weight" found, applied to all columns')
  }

  if(length(w)>1&length(w)<mdim[2]){
    w1=rep(1,mdim[2])
    w1[1:length(w)]=w
    w=w1;rm(w1)
    warning('number of "weight" is less than the number of columns, use 1 as the missing "weight"')
  }
  if(length(w)>mdim[2]){
    w=w[1:mdim[2]]
    warning('number of "weight" is more than the number of columns, omit the extra "weight"')
  }
  
  wmat=matrix(rep(w,each=dim(cov.)[1]*dim(cov.)[2]/mdim[2]),ncol=dim(cov.)[2],byrow=T)

  cov.=wmat*cov.

  cov..=matrix(0,ncol=mdim[1],nrow=mdim.new[1])
  if(mdim[2]>1){
    for(i in 1:(mdim[2]-1)){
      cov..=cov..+cov.[1:mdim.new[1],];cov.=cov.[-(1:mdim.new[1]),]}
    cov.=cov..+cov.
  }
  return(cov.)  
}

Dpow.ex=function(vec,Data,hyper,Q=NULL,gamma){
  DQ=cov.pow.ex(hyper,Data,gamma=gamma);
  DQ=-0.5*DQ*xixj_sta(as.matrix(vec),w=exp(hyper$pow.ex.w[which(apply(Data,2,mean)==mean(vec) & apply(Data,2,max)==max(vec) & apply(Data,2,min)==  min(vec))]),power=gamma)
  return(DQ)
}

Drat.qu=function(vec,Data,hyper,Q=NULL){
  DQ=cov.rat.qu(hyper,Data)
  power=exp(hyper$rat.qu.a)
  DQ=-power*hyper$rat.qu.s*DQ^((power+1)/power)
  DQ=DQ%*%xixj_sta(vec,w=exp(hyper$rat.qu.w[which(apply(Data,2,mean)==mean(vec) & apply(Data,2,max)==max(vec) & apply(Data,2,min)==  min(vec))]),power=1)
  return(DQ)
}


D2pow.ex=function(vec,Data,hyper,Q=NULL,gamma){
  DQ=cov.pow.ex(hyper,Data,gamma=gamma)
  extra=xixj_sta(as.matrix(vec),w=exp(hyper$pow.ex.w[which(apply(Data,2,mean)==mean(vec) & apply(Data,2,max)==max(vec) & apply(Data,2,min)==  min(vec))]),power=gamma)
  D2Q=0.25*DQ*extra^2-0.5*DQ*extra
  return(D2Q)
}

D2rat.qu=function(vec,Data,hyper,Q=NULL){
  Q=cov.rat.qu(hyper,Data)
  power=exp(hyper$rat.qu.a)
  extra=hyper$rat.qu.s*xixj_sta(vec,w=exp(hyper$rat.qu.w[which(apply(Data,2,mean)==mean(vec) & apply(Data,2,max)==max(vec) & apply(Data,2,min)==  min(vec))]),power=1)
  DQ=(-power)*((-power-1)*(extra^2)*Q^((power+2)/power)+extra*Q^((power+1)/power))
  return(DQ)
}

DCov.linear.a=function(hyper,data,AlphaQ){
  Dlinear.aj=apply(data,2,function(i) sum(AlphaQ*exp(hyper$linear.a[which(data[1,]==i[1])])*xixj(as.matrix(i),a=1)) )
  return(Dlinear.aj)
}


DCov.pow.ex.w=function(hyper,data,gamma=1,AlphaQ){
  Dpow.ex.wj=apply(data,2,function(i) sum(AlphaQ*Dpow.ex(as.matrix(i),data,hyper,gamma=gamma)) )
  return(Dpow.ex.wj)
}


DCov.pow.ex.v=function(hyper,data,gamma,AlphaQ){
  DDpow.ex.v=cov.pow.ex(hyper,data)
  Dpow.ex.v=sum(AlphaQ*DDpow.ex.v)
  return(Dpow.ex.v)
}

DCov.rat.qu.w=function(hyper,data,AlphaQ){
  Drat.qu.wj=apply(data,2,function(i) sum(AlphaQ*Drat.qu(i,data,hyper)) )
  return(Drat.qu.wj)
}

DCov.rat.qu.s=function(hyper,data,AlphaQ){
  hyper=lapply(hyper,exp);
  
  v.power=xixj_sta(data,data,hyper$rat.qu.w,power=1)
  Drat.qu.s=(-hyper$rat.qu.a)*v.power*(1+log(hyper$rat.qu.s)*v.power)^(-hyper$rat.qu.a-1)
  Drat.qu.s=sum(AlphaQ*Drat.qu.s)
  return(Drat.qu.s)
}

DCov.rat.qu.a=function(hyper,data,AlphaQ){
  DDrat.qu.a=cov.rat.qu(hyper,data)
  DDrat.qu.a=log(DDrat.qu.a)%*%DDrat.qu.a
  Drat.qu.a=sum(AlphaQ*DDrat.qu.a)
  return(Drat.qu.a)
}

DCov.vv=function(hyper,Alpha,invQ){
  Dfvv=-sum(diag(invQ))*exp(hyper$vv) + t(Alpha)%*%Alpha*exp(hyper$vv)
  return(Dfvv)
}

D2linear.a=function(hyper,data,inv.Q,Alpha.Q){
  D2linear.aj=apply(data,2,function(i) D2(exp(hyper$linear.a[which(data[1,]==i[1])])*xixj(as.matrix(i),a=1),exp(hyper$linear.a[which(data[1,]==i[1])])*xixj(as.matrix(i),a=1),inv.Q,Alpha.Q))
  return(D2linear.aj)
}

D2pow.ex.w=function(hyper,data,gamma,inv.Q,Alpha.Q){
  D2pow.ex.wj=apply(data,2,function(i) D2(Dpow.ex(i,data,hyper,gamma=gamma),D2pow.ex(i,data,hyper,gamma=gamma),inv.Q,Alpha.Q))
  return(D2pow.ex.wj)
}

D2pow.ex.v=function(hyper,data,gamma,inv.Q,Alpha.Q){
  DDpow.ex.v=cov.pow.ex(hyper,data)
  D2pow.ex.v=D2(DDpow.ex.v,DDpow.ex.v,inv.Q,Alpha.Q)  
  return(D2pow.ex.v)
}

D2rat.qu.w=function(hyper,data,inv.Q,Alpha.Q){
  D2rat.qu.wj=apply(data,2,function(i) D2(Drat.qu(i,data,hyper),D2rat.qu(i,data,hyper),inv.Q,Alpha.Q))
  return(D2rat.qu.wj)
}

D2rat.qu.s=function(hyper,data,inv.Q,Alpha.Q){
  hyper=lapply(hyper,exp)
  v.power=xixj_sta(data,data,hyper$rat.qu.w,power=1)
  power=hyper$rat.qu.a
  sD1=(-hyper$rat.qu.a)*v.power*(1+log(hyper$rat.qu.s)*v.power)^(-hyper$rat.qu.a-1)
  sD2=power*(power+1)*v.power^2%*%(sD1/(-power))^((power+1)/(power+2))
  D2rat.qu.s=apply(data,2,function(i) D2(sD1,sD2,inv.Q,Alpha.Q))
  return(D2rat.qu.s)
}
D2rat.qu.a=function(hyper,data,inv.Q,Alpha.Q){
  Q=cov.rat.qu(hyper,data)
  aD1=log(Q)%*%Q
  aD2=log(Q)%*%(log(Q)%*%Q-Q*exp(hyper$rat.qu.a))/exp(hyper$rat.qu.a)
  D2rat.qu.a=D2(aD1,aD2,inv.Q,Alpha.Q)
  return(D2rat.qu.a)
}


D2vv=function(hyper,data,inv.Q,Alpha.Q){
  D2fvv=D2(diag(exp(hyper$vv),dim(data)[1]),diag(exp(hyper$vv),dim(data)[1]),inv.Q,Alpha.Q)
  return(D2fvv)
}


D2=function(d1,d2,inv.Q,Alpha.Q){
  Aii=t(d1)%*%inv.Q%*%d1
  al=Alpha.Q+inv.Q
  return(0.5*(sum(Alpha.Q*(d2-Aii))-sum(al*Aii)))
}

diag.linear=function(hyper,data){
  Qstar=data^2%*%matrix(exp(hyper$linear.a))
  return(Qstar)
}

diag.pow.ex=function(hyper,data){
  Qstar=rep(exp(hyper$pow.ex.v),dim(data)[1])
  return(Qstar)
}

diag.rat.qu=function(hyper,data){
  Qstar=rep(1,dim(data)[1])
  return(Qstar)
}

##################### plot ##########################
plot.gpr=function(x,...,fitted=F,col.no=1){
  obj=x
  if(fitted==T){
    if(is.null(obj$fitted.mean)){
      warning('fitted values not found, ploting predicted values')
      type='Prediction'
      mu=obj$pred.mean
      sd=obj$pred.sd
      x=obj$newdata
      X=obj$train.x
      Y=obj$train.yOri
    }
    if(!is.null(obj$fitted.mean)){
      type='Fitted values'
      mu=obj$fitted.mean
      sd=obj$fitted.sd
      X=obj$train.x
      Y=obj$train.yOri
      x=X
    }
  }
  else{
    if(is.null(obj$pred.mean)){
      warning('predicted values not found, ploting fitted values')
      type='Fitted values'
      mu=obj$fitted.mean
      sd=obj$fitted.sd
      X=obj$train.x
      Y=obj$train.yOri
      x=X
    }
    if(!is.null(obj$pred.mean)){
      type='Prediction'
      mu=obj$pred.mean
      sd=obj$pred.sd
      x=obj$newdata
      X=obj$train.x
      Y=obj$train.yOri
    }
  }
  if(dim(X)[1]<=150|length(X)<=150){
    pchType=4
    PcexNo=0.8
    LcexNo=1.5
  }
  else{
    pchType=20
    PcexNo=0.1
    LcexNo=0.8
  }
  upper=mu+1.96*(sd);
  lower=mu-1.96*(sd);
  plot(-100,-100,col=0,xlim=range(X[,col.no],x[,col.no]),ylim=range(upper,lower,Y),main=type, xlab="input ",ylab="response",...)
  #
  polygon(c(x[,col.no], rev(x[,col.no])), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  #
  points(X[,col.no],Y,pch=pchType,col=2,cex=PcexNo)
  # lines(X[,1],Y)
  lines(x[,col.no],mu,col=4,lwd=LcexNo)  
}