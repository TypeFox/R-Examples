BBRecap=function(
  data,
  last.column.count=FALSE,
  neval=1000,
  by.incr=1,
  mbc.function=c("standard","markov","counts","integer","counts.integer"),
  mod=c("linear.logistic","M0","Mb","Mc","Mcb","Mt","Msubjective.cut","Msubjective"),
  nsim=5000,
  burnin=round(nsim/10),
  nsim.ML=1000,
  burnin.ML=round(nsim.ML/10),
  num.t=50,
  markov.ord=NULL,
  prior.N=c("Rissanen","Uniform","one.over.N","one.over.N2"),
  meaningful.mat.subjective=NULL,
  meaningful.mat.new.value.subjective=NULL,
  z.cut=NULL,
  output=c("base","complete","complete.ML")){

          log.marg.likelihood=NULL
  
  mod=match.arg(mod)
  prior.N=match.arg(prior.N)
  output=match.arg(output)
  mbc.function=match.arg(mbc.function)
  
  model=switch(mod,
               linear.logistic="mod.linear.logistic",
               M0="mod.M0",
               Mb="mod.Mb",
               Mc="mod.Mc",
               Mcb="mod.Mcb",
               Mt="mod.Mt",
               Msubjective.cut="mod.Msubjective.cut",
               Msubjective="mod.Msubjective")  
  
  prior=switch(prior.N,
               Rissanen="Rissanen.prior",
               Uniform="Uniform.prior",
               one.over.N="1overN.prior",
               one.over.N2="1overN^2.prior")	
 
  qb.function=switch(mbc.function,
  					standard="qb.standard",
  					markov="qb.markov",
  					integer="qb.integer",
  					counts="qb.count",
  					counts.integer="qb.count.integer")
  
  if((mod!="Msubjective.cut") & length(z.cut)>0){
 	warning("z.cut without mod='Msubjective.cut' will be IGNORED!")
}
   if(!(class(data)=="data.frame" | class(data)=="matrix" | class(data)=="array" | class(data)=="table")){
       stop("input data must be a data.frame or a matrix object or an array")
   }
  					
 if(class(data)=="table" | class(data)=="array"){
   	
   	n.occ=length(dim(data))
mm=matrix(ncol=n.occ,nrow=2^n.occ)

for(i in 1:2^n.occ){
mm[i,]=as.numeric(intToBits(i-1)>0)[1:n.occ]
}

data.vec=c(data)
data[1]=0
dd=c()
for(i in 1:length(data)){
		
	dd=c(dd,rep(mm[i,],data.vec[i]))
	
		}

data.matrix=matrix(dd,ncol=length(dim(data)),byrow=T)
  	
   }
   
   if(class(data.matrix)!="matrix"){
       data.matrix=as.matrix(data) 
       }
   
  if(last.column.count){

      if (any(data[,ncol(data)]<0)){
    stop("Last column must contain non negative frequencies/counts")
}

  	data=as.matrix(data)
  	data.matrix=matrix(ncol=(ncol(data)-1))

for(i in 1:nrow(data)){

d=rep(data[i,1:(ncol(data)-1)],(data[i,ncol(data)]))
dd=matrix(d,ncol=(ncol(data)-1),byrow=T)
data.matrix=rbind(dd,data.matrix)

}

data.matrix=data.matrix[1:(nrow(data.matrix)-1),]
  	  	
  }

  if(any(data.matrix!=0 & data.matrix!=1)) stop("data must be a binary matrix")
  
   if(sum(apply(data.matrix,1,sum)==0)){
 	warning("input data argument contains rows with all zeros: these rows will be removed and ignored")

        data.matrix=data.matrix[apply(data.matrix,1,sum)!=0,]
    }

  t=ncol(data.matrix)
  M=nrow(data.matrix)
  logistic=0
    
  if(model=="mod.Mc" | model=="mod.Mcb"){
    
    if (!is.numeric(markov.ord) || markov.ord != round(markov.ord) || markov.ord <= 0 || markov.ord > t) 
      stop("The Markov order must be a non-negative integer smaller than t!")
    
    mbc.function="markov"
    z.cut=seq(0,1,by=1/(2^markov.ord))
    z.cut[1]=-0.1
    if(model=="mod.Mcb"){
      z.cut=c(z.cut[1],0,z.cut[-1])
    }
  }

if(qb.function=="qb.standard") recodebinary=quant.binary
if(qb.function=="qb.integer") recodebinary=quant.binary.integer
if(qb.function=="qb.markov") recodebinary=quant.binary.markov
if(qb.function=="qb.count") recodebinary=quant.binary.counts
if(qb.function=="qb.count.integer") recodebinary=quant.binary.counts.integer

#####          
  v=c()

 meaningful.mat=matrix(NA,ncol=ncol(data.matrix),nrow=M)
  
  for(i in 1:ncol(data.matrix)){
    for(j in 1:M){
      v=data.matrix[j,1:i]
      v=paste(v,collapse="")
      if(qb.function=="qb.markov"){
      		meaningful.mat[j,i]=recodebinary(v,markov.ord)}
      else{meaningful.mat[j,i]=recodebinary(v)}	
    }
  }
  
  meaningful.mat=cbind(rep(0,M),meaningful.mat[,1:t-1])
  meaningful.mat.new.value=rep(0,t*by.incr)
  
  if(model=="mod.M0"){
          z.cut=c(-0.1,1)
        }
        
  if(model=="mod.Mb"){
          z.cut=c(-0.1,0,1)
        }
        
  if(model=="mod.Mt"){
          z.cut=seq(0,t)
          meaningful.mat=matrix(rep(1:t),nrow=M,ncol=t,byrow=T)
          meaningful.mat.new.value=1:t        
          }
  
  if(model=="mod.Msubjective"){
    if(is.matrix(meaningful.mat.subjective)){
      if(ncol(meaningful.mat.subjective)!=ncol(data.matrix)||
           nrow(meaningful.mat.subjective)!=nrow(data.matrix))
        stop("meaningful.mat.subjective has different dimensions from data matrix")
      meaningful.mat=meaningful.mat.subjective  
    }
    else{ stop("meaningful.mat.subjective must be a numerical matrix") }
    if(is.numeric(meaningful.mat.new.value.subjective)){
      meaningful.mat.vec.new=meaningful.mat.new.value.subjective
    } 
  }
  
  meaningful.mat.vec=c(meaningful.mat)
  y=c(data.matrix)
  
  nn.1=c()
  nn.0=c()
  
  coeff.vec1=c()
  coeff.vec2=c()
  
  log.post.N=c()
  post.N=c()
  vv=c()
  
  prior.inv.const=0
  
  prior.distr.N=function(x){
    
    if(prior=="Rissanen.prior") {out=(rissanen(x))}
    if(prior=="Uniform.prior") {out=(1/(x^0))}
    if(prior=="1overN.prior") {out=(1/(x^1))}
    if(prior=="1overN^2.prior") {out=(1/(x^2))}
    
    return(out)
  }
    
  if(model=="mod.linear.logistic" | model=="mod.Msubjective"){logistic=1}
  
  if(logistic==0){

            if(output=="complete.ML"){
        output="complete"
       } 
    
    for(k in 1:neval){
 meaningful.mat.new=rbind(meaningful.mat,matrix(rep(meaningful.mat.new.value,(k>1)),nrow=((k-1) *by.incr),ncol=t,byrow=T))
 dati.new=rbind(data.matrix,matrix(0,nrow=((k-1)*by.incr),ncol=t))
 meaningful.mat.vec.new=c(meaningful.mat.new)
 y.new=c(dati.new)

        for(j in 1:(length(z.cut)-1)){
          
          if(j==1){nn.0[j]=sum(meaningful.mat.vec.new>=z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==0)}          
          else{nn.0[j]=sum(meaningful.mat.vec.new>z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==0)} 
          if(j==1){nn.1[j]=sum(meaningful.mat.vec.new>=z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==1)}          
          else{nn.1[j]=sum(meaningful.mat.vec.new>z.cut[j] &
          	meaningful.mat.vec.new<=z.cut[j+1] & y.new==1)} 
        }
                
        prior.inv.const=prior.inv.const+ prior.distr.N((M+k-1))
        
        log.post.N[k]=lchoose((sum(nn.1)+sum(nn.0))/t,M)+sum(lbeta((nn.1+1),(nn.0+1)))+log(prior.distr.N((M+k-1)))
        
      }
      
     # closed k cycle

      l.max=max(log.post.N)
      
      for(k in 1:neval){
        vv[k]<-exp(log.post.N[k]-l.max)
      }
      
      ss=sum(vv)
      post.N=vv/ss

	  if(output=="complete"){
      if(model=="mod.Mt"){
      	Exp.beta=matrix(ncol=(length(z.cut)-1),nrow=neval)   
      	for(i in 1:neval){     	
      	Exp.beta[i,]=(nn.1+1)/(nn.1+nn.0-(neval-i)+2)	
 	 	} 
      	pH.post.mean=apply((post.N*Exp.beta),2,sum)      	
      }
      else{
      	Exp.beta=matrix(ncol=(length(z.cut)-1),nrow=neval)   
      	for(i in 1:neval){     	
      	Exp.beta[i,]=(nn.1+1)/(nn.1+nn.0+2)	
 	 	} 
  	 	rr=sort(0:(neval-1)*t,decreasing=T)
 	 	Exp.beta[,1]=Exp.beta[,1]*(nn.1[1]+nn.0[1]+2)/(nn.1[1]+nn.0[1]-rr+2)
      	pH.post.mean=apply((post.N*Exp.beta),2,sum)            	
      	}
    
    names(pH.post.mean)=paste("pH",1:(length(z.cut)-1),sep="")
      }
		 
    val=seq(M,((M+(neval-1)*by.incr)),by=by.incr)
    
    ord=order(post.N,decreasing=T)
    p.max=ord[1]
    mode.N=val[p.max]
    
    mean.N<-round(sum(post.N*val))
    
    median.N=M
    
    g<-c()				
    for(k in 1:neval)	{
      g=c(g,post.N[k])
      if (sum(g)<=0.5) median.N=median.N+1
    }
    
    funzioneRMSE=function(x){
      sum((((x/val)-1)^2)*post.N)
    }		
    
    estimate.N=round(optimize(funzioneRMSE,c(min(val),max(val)))$minimum)	
    
    ### Credible Set ###
    
    alpha<-0.05	
    
    g=0
    d=0
    
    aa=c()
    
    ordine<-order(post.N,decreasing=T)
    
    w<-val
    w<-w[ordine]
    p<-post.N
    p<-p[ordine]
    
    for(k in 1:neval)		{
      if (g<(1-alpha)) 	{g=g+p[k]  
                         d=d+1}
    }
    aa<-w[1:d]
    inf.lim<-min(aa)
    sup.lim<-max(aa)
    
    log.marg.likelihood=log(sum(exp(log.post.N-max(log.post.N)-log(prior.inv.const))))+max(log.post.N)
    
  }
  
  ## ARMS
  
  if(model=="mod.linear.logistic" | model=="mod.Msubjective"){logistic=1}
  
  if(logistic==1){
    
    print("implementation via ARMS")
    
    ### LOG-PRIOR N
    
    log.prior.N=function(x){
      
      if(prior=="Rissanen.prior") {out=log(rissanen(x))}
      if(prior=="Uniform.prior") {out=0}
      if(prior=="1overN.prior") {out=log(1/(x^1))}
      if(prior=="1overN^2.prior") {out=log(1/(x^2))}
      
      return(out)
    }
    
    ### LOG-LIKELIHOOD
    
    if(logistic==1){
      
      log.likelihood=function(NN,aa0,aa1){
        
        z=matrix(0,ncol=t,nrow=NN)
        x=matrix(0,ncol=t,nrow=NN)
        
        z[(1:M),]=meaningful.mat
        x[(1:M),]=data.matrix
        l = x*log(exp(aa0+aa1*z)/(1+exp(aa0+aa1*z)))-(1-x)*log((1+exp(aa0+aa1*z)))
        out=sum(l)	+lchoose(NN,M)	
        return(out)
      }
      ###
      
      N.vec=c(M,rep(0,(nsim-1)))
      a0.vec=c(0,rep(0,(nsim-1)))
      a1.vec=c(0,rep(0,(nsim-1)))
      
      ### MODIFIED-LOG-LIKELIHOOD 
      
      log.likelihood.t=function(NN,aa0,aa1,tt){
        
        z=matrix(0,ncol=t,nrow=NN)
        x=matrix(0,ncol=t,nrow=NN)
        
        z[(1:M),]=meaningful.mat
        x[(1:M),]=data.matrix
        l = x*log(exp(aa0+aa1*z)/(1+exp(aa0+aa1*z)))-(1-x)*log((1+exp(aa0+aa1*z)))
        out=tt*(sum(l)	+lchoose(NN,M))	
        return(out)
      }	
      ###
      
      for(g in 2:(nsim+1)){
        
        N.vec[g]=arms(
          N.vec[g-1],
          function(N) log.prior.N(N)+log.likelihood(N,a0.vec[g-1],a1.vec[g-1]),
          function(N) (N>=M)*(N<=M+neval),
          1
        )
        N.vec[g]=round(N.vec[g])
        
        a0.vec[g]=arms(
          a0.vec[g-1],
          function(a0) log.likelihood(N.vec[g],a0,a1.vec[g-1]),
          function(a0) (a0>-5)*(a0<=5),
          1
        )
        
        a1.vec[g]=arms(
          a1.vec[g-1],
          function(a1) log.likelihood(N.vec[g],a0.vec[g],a1),
          function(a1) (a1>-5)*(a1<=5),
          1
        )
        
      }
      
      mean.N=round(mean(N.vec[-(1:burnin)]))
      mean.a0=mean(a0.vec[-(1:burnin)])
      mean.a1=mean(a1.vec[-(1:burnin)])
      median.N=median(N.vec[-(1:burnin)])
      
      mode.N=as.numeric(names(sort(table(N.vec[-(1:burnin)]))[length(table(N.vec[-(1:burnin)]))]))
      
      val=seq(M,(M+neval))
      
      post.N=(table(c(N.vec,M:(M+neval)))-1)/sum(N.vec)
      
      funzioneRMSE=function(x){
        sum((((x/val)-1)^2)*post.N)
      }		
      
      estimate.N=round(optimize(funzioneRMSE,c(min(val),max(val)))$minimum)	    
      HPD.interval=function(x.vec){
        obj=density(x=x.vec,from=min(x.vec),to=max(x.vec),n=diff(range(x.vec))+1,kernel="rectangular")
        density.estimation.N=obj$y
        names(density.estimation.N)=obj$x
        ordered.x.vec=x.vec[order(density.estimation.N[as.character(x.vec)],decreasing=T)]
        HPD.interval=range(ordered.x.vec[1:round(quantile(1:length(x.vec),prob=c(0.95)))])
        return(HPD.interval)
      }
      
      HPD.N=HPD.interval(N.vec[-(1:burnin)])
      HPD.a0=HPD.interval(a0.vec[-(1:burnin)])
      HPD.a1=HPD.interval(a1.vec[-(1:burnin)])
      
      if(output=="complete.ML"){

        output="complete"
        
        a=seq(0.0016,1,length=num.t)
        
        a=a^4	
        
        l=matrix(ncol=length(a),nrow=nsim.ML)
        
        for(j in 1:length(a)){

  print(paste("ML evaluation via power-posterior method of Friel & Pettit (2008) ; t =",j,"of",length(a)))
          
          N.vec.ML=c(M,rep(0,(nsim-1)))
          a0.vec.ML=c(0,rep(0,(nsim-1)))
          a1.vec.ML=c(0,rep(0,(nsim-1)))
          
          for(g in 2:(nsim.ML+1)){
            
            N.vec.ML[g]=arms(
              N.vec.ML[g-1],
              function(N) log.prior.N(N)+log.likelihood.t(N,a0.vec.ML[g-1],a1.vec.ML[g-1],a[j]),
              function(N) (N>=M)*(N<=M+neval),
              1
            )
            
            N.vec.ML[g]=round(N.vec.ML[g])
            
            a0.vec.ML[g]=arms(
              a0.vec.ML[g-1],
              function(a0) log.likelihood.t(N.vec.ML[g],a0,a1.vec.ML[g-1],a[j]),
              function(a0) (a0>-5)*(a0<=5),
              1
            )
            
            a1.vec.ML[g]=arms(
              a1.vec.ML[g-1],
              function(a1) log.likelihood.t(N.vec.ML[g],a0.vec.ML[g],a1,a[j]),
              function(a1) (a1>-5)*(a1<=5),
              1
            )
            
            l[g-1,j]=log.likelihood(N.vec.ML[g],a0.vec.ML[g],a1.vec.ML[g])
            
          }	
        }	
        
        ll=apply(l,2,mean)	
        
        v=c()
        for(k in 1:(length(a)-1)){	
          v[k]=(a[k+1]-a[k])*(ll[k+1]+ll[k])/2	
        }
        
        log.marg.likelihood=sum(v)
        
      }
    }	
  }
  
  if(logistic==0){
    if(mod=="M0" | mod=="Mb" | mod=="Mt" | mod=="Mc" | mod=="Mcb"){
    out=switch(output,
               base=
                 list(Model=model,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim)),
                 complete=
                 list(Model=model,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),pH.post.mean=pH.post.mean,                    N.range=val,posterior.N=post.N,z.matrix=meaningful.mat,vec.cut=z.cut,log.marginal.likelihood=log.marg.likelihood)
                 ) 
  }
 else{
      out=switch(output,
               base=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim)),
                 complete=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=c(inf.lim,sup.lim),pH.post.mean=pH.post.mean, N.range=val,posterior.N=post.N,z.matrix=meaningful.mat,vec.cut=z.cut,log.marginal.likelihood=log.marg.likelihood
                 ) )}
}  
  
  if(logistic==1) {
   if(mod=="Msubjective"){ 
    out=switch(output,
               base=
                 list(Model=model,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N),
               complete=
                 list(Model=model,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat,log.marginal.likelihood=log.marg.likelihood))
    }
  
    else{  
    	out=switch(output,
               base=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N),
               model=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat),
               complete=
                 list(Model=model,mbc.function=mbc.function,prior.N=prior.N,N.hat.mean=mean.N,N.hat.median=median.N,N.hat.mode=mode.N,N.hat.RMSE=estimate.N,HPD.N=HPD.N,N.vec=N.vec,mean.a0=mean.a0,HPD.a0=HPD.a0,a0.vec=a0.vec,mean.a1=mean.a1,HPD.a1=HPD.a1,a1.vec=a1.vec,z.matrix=meaningful.mat,log.marginal.likelihood=log.marg.likelihood))}
  
  }
  
  return(out)
  
}
