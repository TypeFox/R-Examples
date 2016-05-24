supcluster<-function(data,outcome="outcome",features=1:(dim(data)[2]-1),log.transform=TRUE,maxclusters=10,nstart=100,
                     n=500,shape=1,scale=1,
                     alpha=1,betaP=1,fixj="random",
                     fbeta=FALSE,starting.value=NULL,nchains=1){
  subcluster.env=environment()
  
  #if(exists("jprimej)")){remove(jprimej)}
  #we assume that the outcome has the name outcome
  #and is the last variable, the rest are genes.
  if (!fbeta){
    x=data[,outcome]
    x=x-mean(x)}
  y=data[,features]
  if (log.transform) y=log(y)
 
  dimY=dim(y)
  ngenes=dimY[2]
  dimData=c((dimY[1]),ngenes+1)
  maxmatrix=matrix(1:maxclusters,maxclusters,1)
 #extract the columns of the matrix of genes
  prod.dimY=prod(dimY)
  #y=y-matrix(rowMeans(y),dimData[1],dimY[2])
 y=y-matrix(colMeans(y),dimData[1],dimY[2],byrow=TRUE)
  #yxmaxclusters=matrix(y,dimY)??
  y2=y*y
  ty=t(y)
  ty2=t(y2)
  y2g=colSums(y2)
  yg=colSums(y)
 
#diagTM=cmpfun(diagTM1)
 diagTM=function(a,b,s=TRUE){
   mb=dim(b)
   if(s){
     return(sum(b*rep(a,mb[2])))
   } else{
     return(b*matrix(a,mb[1],mb[2]))}
 }
sigma<-function(xx){
  J=xx$jmatrix
  mat2=xx$tau*J%*%xx$beta
  Sigma=rbind(cbind(diag(xx$sig,dimData[2]-1,dimData[2]-1)+
                      tcrossprod(xx$tau*J,J),
                    mat2),cbind(t(mat2),
                                xx$tau*sum(xx$beta^2)+xx$gamma))
  return(Sigma)
}
#OLD function too slow 
loglike1=function(parms){
  #calculates the log-likelihood using 3.
  #generate the variance covariance matrix.
  #note that cluster membership is define dy the J matrix
  #jmatrix[j,i]=1 iff gene j is in cluster i.
  loglik=0
  for (i in 1:dimData[1])
    loglik=loglik+dmvnorm(cbind(y[i,],x[i]),sigma=sigma(parms),log=TRUE)
  return(loglik)
} 
#end loglike

  jprimej=ydot=y2dot=ydot2=jpd=ex1=NULL

updatej=function(parms){
    assign("jprimej",value=colSums(parms$jmatrix),envir=subcluster.env)
    assign("ydot",value=crossprod(parms$jmatrix,ty),envir=subcluster.env)
    assign("y2dot",value=rowSums(crossprod(parms$jmatrix,ty2)),envir=subcluster.env)
    assign("ydot2",value=rowSums(ydot^2),envir=subcluster.env)
  }
  updateb1=function(parms){
    denom=1/(parms$sig+parms$tau*jprimej)
    if(!fbeta){
    assign("jpd",value=sum(jprimej*parms$beta^2*denom),envir=subcluster.env)
    assign("ex1",value=colSums(diagTM(c(parms$beta)*denom,ydot,s=FALSE)),envir=subcluster.env)}
  }
  
  loglike=function(parms,updateJ=FALSE,updateb=FALSE){
    #jprimej=tpj%*%(parms$jmatrix)
  if(updateJ) {updatej(parms);updateb1(parms)}
  if(!updateJ&updateb) updateb1(parms)
    denom=1/(parms$sig+parms$tau*jprimej)
    if(!fbeta){
    ex=parms$tau*ex1
    var.x=parms$tau*sum(parms$beta^2)+parms$gamma-(parms$tau^2)*jpd
    
    loglik=-.5*(prod.dimY*log(2*pi)+prod.dimY*log(parms$sig)+
                  dimY[1]*sum(log(1+jprimej*parms$tau/parms$sig))+ #1
                  #diagTM(denom,y2dot)+
                  sum(denom*y2dot)+ #2
                  parms$tau*(sum(jprimej*denom*y2dot)-#3
                  sum(denom*ydot2))/parms$sig+#4
                  dimY[1]*(log(2*pi)+log(var.x))+ #5
                  sum((x-ex)^2)/var.x)} #6
                  else {loglik=-.5*(prod.dimY*log(2*pi)+prod.dimY*log(parms$sig)+
                                       dimY[1]*sum(log(1+jprimej*parms$tau/parms$sig))+ #1
                                       #diagTM(denom,y2dot)+
                                       sum(denom*y2dot)+ #2
                                       parms$tau*(sum(jprimej*denom*y2dot)-#3
                                                    sum(denom*ydot2))/parms$sig)}
    return(loglik)
  }
 
  #loglike=cmpfun(loglike1)
 dloglike=function(parms,j,i1,i2){
   # new-old loglikelihood when you move gene j to cluster i1 to i2
   if (i1==i2) return(0)
   jprimej.i1=jprimej[i1]
   jprimej.i2=jprimej[i2]
   y2dot.i1=y2dot[i1]
   y2dot.i2=y2dot[i2]
   ##reduce hash look up
   sig=parms$sig
   tau=parms$tau
   if(!fbeta){
   beta=parms$beta
   beta.i1=beta[i1]
   beta.i2=beta[i2]}
   idenom.i1=( sig+ tau*jprimej.i1)
   idenom.i2=( sig+ tau*jprimej.i2)
   new.idenom.i1=( sig+ tau*(jprimej.i1-1))
   new.idenom.i2=( sig+ tau*(jprimej.i2+1))
   if (!fbeta){
   var.x1= tau*sum( beta^2)+parms$gamma
   var.x=var.x1-( tau^2)*jpd
   new.ex= tau*(ex1-( beta.i1/idenom.i1)*ydot[i1,]-( beta.i2/idenom.i2)*ydot[i2,] 
                     +( beta.i1/new.idenom.i1)*(ydot[i1,]-ty[j,])+
                                      ( beta.i2/new.idenom.i2)*(ydot[i2,]+ty[j,]))
   new.varx=var.x1-( tau^2)*(jpd-jprimej.i1* beta.i1^2/( sig+ tau*jprimej.i1)+
                                     (jprimej.i1-1)* beta.i1^2/( sig+ tau*(jprimej.i1-1))-
                                     jprimej.i2* beta.i2^2/( sig+ tau*jprimej.i2)+
                                     (jprimej.i2+1)* beta.i2^2/( sig+ tau*(jprimej.i2+1)))
   }
   if (!fbeta){
   dloglike= -.5*(
              dimY[1]*(log(1+(jprimej.i1-1)* tau/ sig)+log(1+(jprimej.i2+1)* tau/ sig)-
              (log(1+(jprimej.i1)* tau/ sig)+log(1+(jprimej.i2)* tau/ sig)) )+ #1
        (y2dot.i1-y2g[j])/new.idenom.i1-y2dot.i1/idenom.i1+(y2dot.i2+y2g[j])/new.idenom.i2-y2dot.i2/idenom.i2+#2
         tau*((jprimej.i1-1)*(y2dot.i1-y2g[j])/new.idenom.i1-jprimej.i1*y2dot.i1/idenom.i1+
                     (jprimej.i2+1)*(y2dot.i2+y2g[j])/new.idenom.i2-jprimej.i2*y2dot.i2/idenom.i2-#3
          (sum((ydot[i1,]-ty[j,])^2)/new.idenom.i1-ydot2[i1]/idenom.i1+
           sum((ydot[i2,]+ty[j,])^2)/new.idenom.i2-ydot2[i2]/idenom.i2))/ sig+#4
          dimY[1]*(log(new.varx)-log(var.x))+#5
          sum((x-new.ex)^2)/new.varx-sum((x-tau*ex1)^2)/var.x)}  #6
         else {
           dloglike= -.5*(
             dimY[1]*(log(1+(jprimej.i1-1)* tau/ sig)+log(1+(jprimej.i2+1)* tau/ sig)-
                        (log(1+(jprimej.i1)* tau/ sig)+log(1+(jprimej.i2)* tau/ sig)) )+ #1
               (y2dot.i1-y2g[j])/new.idenom.i1-y2dot.i1/idenom.i1+(y2dot.i2+y2g[j])/new.idenom.i2-y2dot.i2/idenom.i2+#2
               tau*((jprimej.i1-1)*(y2dot.i1-y2g[j])/new.idenom.i1-jprimej.i1*y2dot.i1/idenom.i1+
                      (jprimej.i2+1)*(y2dot.i2+y2g[j])/new.idenom.i2-jprimej.i2*y2dot.i2/idenom.i2-#3
                      (sum((ydot[i1,]-ty[j,])^2)/new.idenom.i1-ydot2[i1]/idenom.i1+
                         sum((ydot[i2,]+ty[j,])^2)/new.idenom.i2-ydot2[i2]/idenom.i2))/ sig)#4
              }  
               return(dloglike)
         }
 
#not been updated


#iteration function
  itgamma<-function(parms,name){ 
    #save old parameter value
    start=parms
    sig=1/rgamma(1,shape=shape,scale=(1/start[[name]])/scale)#proposal
    parms[name]=sig #set new value
    lb2=loglike(start)
    lb1=loglike(parms,updateb=TRUE)
    vs=min(1,exp(lb1-lb2-
                   (dgamma(1/sig,shape=shape,scale=(1/start[[name]])/scale,log=TRUE)-
                      dgamma(1/start[[name]],shape=shape,scale=(1/sig)/scale,log=TRUE))))
    prior=((sig)/(start[[name]]))^(3/2) #I don;t understand this
    if (runif(1)<vs*prior) return(parms) else {updateb1(start);return(start)}
  }
  
  itbeta<-function(parms,clusters=NULL){
    updatej(parms)
    if(!fbeta){
      if(is.null(clusters)){clusters=1:maxclusters}
      for (i in clusters){
      betan=rnorm(n=1,parms$beta[i],1/betaP)#proposal distribution beta
      parmsnew=parms
      parmsnew$beta[i]=betan
      lb2=loglike(parms,updateJ=TRUE)
      lb1=loglike(parmsnew,updateJ=TRUE)
      vs=min(1,exp(lb1-lb2))
      dec=(runif(1)<vs)
      if (dec) {parms=parmsnew} else {updatej(parms);updateb1(parms)}
      }
      }
    return(parms)
  }
  
  #updates the jmatrix and the beta's 
  itj<-function(parms){
    if(!is.numeric(fixj)){
      
      for (j in 1:(dimData[2]-1)){
        iold=match(1,parms$jmatrix[j,]) 
        logprobs=rep(0,maxclusters)        
        for (i in 1:maxclusters){
          new.jprimej=jprimej
          new.jprimej[i]=jprimej[i]+1
          new.jprimej[iold]=jprimej[i]-1
          #conceptual issue here do I weight each gene or use the multivariate
          #probabiity for the new determination
          logprobs[i]=dloglike(parms,j,iold,i)+dmultinom(new.jprimej,prob=parms$wts,log=TRUE)}
 
        vs=logprobs-max(logprobs)
        probs=exp(vs)/sum(exp(vs))
        #select new clustering for gene j
        js=sample.int(maxclusters,size=1,prob=probs)
        #update jmatrix for gene j
        parms$jmatrix[j,]=rep(0,maxclusters)
        parms$jmatrix[j,js]=1
        #conceptual what is the prior for Beta or can I just ignore this.
        #conceptual issue how do you propose for beta, do you assume that beta has 
        #dimension maxclusters
        #conceptual issue I do this separately for each gene. 
    if (iold!=js) parms=itbeta(parms,c(iold,js))
      }
      ww1=colSums(parms$jmatrix)
      parms$wts=rdirichlet(1,ww1+rep(alpha,maxclusters))
    }
    return(parms)  
  } 
  iter<-function(parms,name){  
    if (name != 'jmatrix'){
      return(itgamma(parms,name))
      
    } else {
      #update jmatrix
      return(itj(parms))    
    }
  }
  #jmatrix,beta,parm=nue,sig,tau,beta,gamma
  #starting values
  #conceptual issue what were her starting values
  output=list()
  for (i in 1:nchains) output[[i]]=list(inp=c(0,0),parms=data.frame())
  nreps=n-nstart+1
  ss=rep(0,nreps)
  outpt=data.frame(sig=ss,tau=ss,gamma=ss,beta=matrix(0,nreps,maxclusters),genes=matrix(0,nreps,ngenes))
  for (stream in 1:nchains){
  output[[stream]]$parms=outpt
  output[[stream]]$inp=c(maxclusters=maxclusters,ngenes=ngenes)
  #create starting values
  if(fixj=='kmeans'&&(nchains==1)){
    clus=kmeans(t(y),maxclusters-1,
                algorithm="Hartigan-Wong")
    jmatrix=matrix(0,dimData[2]-1,maxclusters)
    for (i in 1:(dimData[2]-1)) jmatrix[i,clus$cluster[i]]=1
    } else
      { if(is.numeric(fixj)) jmatrix=fixj else {
        jmatrix=matrix(0,dimData[2]-1,maxclusters)
        for (i in 1:(dimData[2]-1)){
          
          jmatrix[i,sample.int(maxclusters,1)]=1}
        }
    }
  if(!is.null(starting.value)){
    jmatrix=matrix(0,dimData[2]-1,maxclusters)
    for (i in 1:(dimData[2]-1)){
         jmatrix[i,starting.value[3+maxclusters+i]]=1
    }
  sig=starting.value[1]
  tau=starting.value[2]
  gamma=starting.value[3]
  beta=matrix(starting.value[3+1:(maxclusters)],maxclusters,1)
  parms=list(jmatrix=jmatrix,sig=sig,tau=tau,gamma=gamma,beta=beta,wts=
               colMeans(jmatrix)/sum(colMeans(jmatrix)))
  } else {
  parms=list(jmatrix=jmatrix,sig=1,
             tau=1,gamma=1,beta=matrix(0,maxclusters,1),wts=
               rep(1/maxclusters,maxclusters))}
  
  #iterate ds
  updatej(parms)
  updateb1(parms)
  for (ii in 1:n){
    #print(parms$sig)
    parms=iter(parms,"sig")
    parms=iter(parms,"tau")
    if(!fbeta) parms=iter(parms,"gamma")
    parms=iter(parms,'jmatrix')
    parms=itbeta(parms)
    if (!(ii<nstart)){
    jj=ii-nstart+1
    output[[stream]]$parms[jj,]=c(parms$sig,parms$tau,parms$gamma,parms$beta,
                                      t(parms$jmatrix%*%matrix(1:maxclusters,maxclusters,1)))}
      }
  }
  return(output)
}

