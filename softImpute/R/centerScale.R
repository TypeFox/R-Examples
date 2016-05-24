cs.default <-
  function(x,maxit,thresh,row.center,row.scale,col.center,col.scale,trace,m,n,alpha,beta,tau,gamma){
    centerscaleD=function(x,tau,gamma,alpha,beta){
      xhat=scale(scale(x,beta,FALSE)-alpha,FALSE,gamma)/tau
      attr(xhat,"scaled:center")=NULL
      attr(xhat,"scaled:scale")=NULL
      xhat
    }
    critmat=NULL;
    crit=1;iter=0
    while((crit>thresh)&(iter<maxit)){
      iter=iter+1
### centering
###column mean
      if(col.center){
        xc=x-alpha
        dbeta=beta
        beta=apply(xc,2,weighted.mean,w=1/tau,na.rm=TRUE)
        dbeta=beta-dbeta
      }
      else dbeta=0
### row mean
      if(row.center){
        xc=scale(x,beta,FALSE)
        dalpha=alpha
        alpha=apply(xc,1,weighted.mean,w=1/gamma,na.rm=TRUE)
        dalpha=alpha-dalpha
      }
      else dalpha=0
###Scaling
      dgamma=1
      dtau=1
      if(row.scale|col.scale){
        xc=scale(x-alpha,beta,FALSE)
### Column scales
        if(col.scale){
          dgamma=gamma
          gamma=sqrt(apply( (xc/tau)^2,2,mean,na.rm=TRUE))
          dgamma=gamma/dgamma
        }
### Row scales
      if(row.scale){
        dtau=tau
        tau=sqrt(apply(scale(xc,FALSE,gamma)^2,1,mean,na.rm=TRUE))
        dtau=tau/dtau
      }

###Check
      }
        xhat=centerscaleD(x,tau,gamma,alpha,beta)
        crit=sum(dalpha^2)+sum(dbeta^2)+sum(log(dtau)^2)+sum(log(dgamma)^2)      
      
        critmat=rbind(critmat,c(iter=iter,crit=crit))
        if(trace)cat("Iter",iter,"Total Changes",crit,"\n")
      }

      alist=attributes(xhat)
     alist=c(alist,list("biScale:row"=list(center=alpha,scale=tau),"biScale:column"=list(center=beta,scale=gamma),critmat=critmat))
      attributes(xhat)=alist
xhat
  }

cs.sparseMatrix <-
  function(x,maxit,thresh,row.center,row.scale,col.center,col.scale,trace,m,n,alpha,beta,tau,gamma){
      if(!inherits(x,"sparseMatrix"))stop("x should be in 'sparseMatrix' format")
      if(!inherits(x,"dgCMatrix"))x=as(x,"dgCMatrix")
      centerscaleC=function(x,tau,gamma,alpha,beta,m,n){
        irow=x@i
        pcol=x@p
        sfac=suvC(as.matrix(tau),as.matrix(gamma),irow,pcol)
        x@x=x@x/sfac
        alpha=cbind(alpha,rep(1,m))/tau
        colnames(alpha)=NULL
        beta=cbind(rep(1,n),beta)/gamma
        colnames(beta)=NULL
        splr(x,alpha,-beta)
      }
      xhat=centerscaleC(x,tau,gamma,alpha,beta,m,n)
      critmat=NULL;crit=1;iter=0
      while((crit>thresh)&(iter<maxit)){
        iter=iter+1
### centering
###column mean
      if(col.center){
        dbeta=gamma*colSums(xhat)/sum(1/tau)
        beta=beta+dbeta
        xhat=centerscaleC(x,tau,gamma,alpha,beta,m,n)
      }
      else dbeta=0
### row mean
      if(row.center){
        dalpha=tau*rowSums(xhat)/sum(1/gamma)
        alpha=alpha+dalpha
        xhat=centerscaleC(x,tau,gamma,alpha,beta,m,n)
      }
      else dalpha=0
###Scaling
### Column scales
      if(col.scale){
        dgamma=sqrt(colsum2.splr(xhat)/m)
        gamma=gamma*dgamma
        xhat=centerscaleC(x,tau,gamma,alpha,beta,m,n)
        }
      else dgamma=1

### Row scales
      if(row.scale){
        dtau=sqrt(rowsum2.splr(xhat)/n)
        tau=tau*dtau
        xhat=centerscaleC(x,tau,gamma,alpha,beta,m,n)
      }
      else dtau=1
###Check
   crit=sum(dalpha^2)+sum(dbeta^2)+sum(log(dtau)^2)+sum(log(dgamma)^2)      
      
      critmat=rbind(critmat,c(iter=iter,crit=crit))
      if(trace)cat("Iter",iter,"Total Changes",crit,"\n")
    }
      alist=attributes(xhat)
     alist=c(alist,list("biScale:row"=list(center=alpha,scale=tau),"biScale:column"=list(center=beta,scale=gamma),critmat=critmat))
      attributes(xhat)=alist
xhat
  }


cs.Incomplete <-
  function(x,maxit,thresh,row.center,row.scale,col.center,col.scale,trace,m,n,alpha,beta,tau,gamma){
### Incomplete is sparseMatrix, where the 0's correspond to NAs
      if(!inherits(x,"sparseMatrix"))stop("x should be in 'sparseMatrix' format")
      if(!inherits(x,"dgCMatrix"))x=as(x,"dgCMatrix")
      irow=x@i
      pcol=x@p

      rowsum.along <-function(b,x){
        ##b is a n vector
        x@x[]=1
        drop(x%*%b)
      }
      colsum.along=function(a,x){
        ##a is a m vector
        x@x[]=1
        drop(t(a)%*%x)
      }

      centerscaleI=function(x,tau,gamma,alpha,beta,m,n,irow=x@i,pcol=x@p){
        alpha=cbind(alpha,rep(1,m))
        beta=cbind(rep(1,n),beta)
        x@x=(x@x-suvC(alpha,beta,irow,pcol))/suvC(as.matrix(tau),as.matrix(gamma),irow,pcol)
        x
      }
      xhat=centerscaleI(x,tau,gamma,alpha,beta,m,n,irow,pcol)
    critmat=NULL;crit=1;iter=0
    while((crit>thresh)&(iter<maxit)){
      iter=iter+1
### centering
###column mean
      if(col.center){
        dbeta=gamma*colSums(xhat)/colsum.along(1/tau,x)
        beta=beta+dbeta
        xhat=centerscaleI(x,tau,gamma,alpha,beta,m,n,irow,pcol)
      }
      else dbeta=0
### row mean
      if(row.center){
        dalpha=tau*rowSums(xhat)/rowsum.along(1/gamma,x)
        alpha=alpha+dalpha
        xhat=centerscaleI(x,tau,gamma,alpha,beta,m,n,irow,pcol)
      }
      else dalpha=0
###Scaling
### Column scales

      if(col.scale){
        dgamma=sqrt(colSums(xhat^2)/colsum.along(rep(1,m),x))
        gamma=gamma*dgamma
        xhat=centerscaleI(x,tau,gamma,alpha,beta,m,n,irow,pcol)
        }
      else dgamma=1

### Row scales
###    a=1/sqrt(apply(scale(xc,FALSE,1/b)^2,1,mean,na.rm=TRUE))
      if(row.scale){
        dtau=sqrt(rowSums(xhat^2)/rowsum.along(rep(1,n),x))
        tau=tau*dtau
        xhat=centerscaleI(x,tau,gamma,alpha,beta,m,n,irow,pcol)
      }
      else dtau=1
###Check
   crit=sum(dalpha^2)+sum(dbeta^2)+sum(log(dtau)^2)+sum(log(dgamma)^2)      
      
      critmat=rbind(critmat,c(iter=iter,crit=crit))
      if(trace)cat("Iter",iter,"Total Changes",crit,"\n")
    }
      alist=attributes(xhat)
     alist=c(alist,list("biScale:row"=list(center=alpha,scale=tau),"biScale:column"=list(center=beta,scale=gamma),critmat=critmat))
      attributes(xhat)=alist
xhat
  }


setGeneric("centerScale",cs.default)
setMethod("centerScale","sparseMatrix",cs.sparseMatrix)
setMethod("centerScale","Incomplete",cs.Incomplete)


