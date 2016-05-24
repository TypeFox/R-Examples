
solqp=function(H=NULL,f,A,B,X,neqtl,nctl,numVar){
  f=as.matrix(f)
  B=as.matrix(B)
  numVar=length(X)
  flag = 1
  iter = 0
  iterflag = 0
  if (is.null(H)){
    qpflag=0
  } else if (norm(H,'I')==0){
    qpflag=0
  }else{
    qpflag=1
  }
  
  imax = 10*numVar
  tolcon  = 10e-5
  Astan = matrix(sqrt(apply(A^2,1,sum)),ncol=1)
  B=B/Astan
  A=A/matrix(rep(Astan,ncol(A)),ncol=ncol(A))
  tolerr = 0.01*sqrt(.Machine$double.eps)
  tolf = 100*numVar*.Machine$double.eps
  
  orignctl = nctl;
  lambda = matrix(0,orignctl,1)
  if(neqtl>0){
    indxact=ieq = 1:neqtl
  }else{
    indxact=ieq = NULL
  }
  
  xieq = matrix(0,nctl,1)
  xieq[indxact] = 1
  iact = length(indxact)
  Act = A[indxact,,drop=F]
  indall = 1:nctl
  idel = NULL
  if (iact > 0 ){
    tolcons = 1e-10
    Z=NULL
    s=qr(A[ieq,])
    Qa=qr.Q(s)
    Ra=qr.R(s)
    idxdep =which(abs(diag(Ra))<tolf)
    if (neqtl > numVar){
      idxdep = rbind(idxdep, matrix((numVar+1):neqtl,ncol=1))
    }
    if (length(idxdep)>0){
      flag = 1
      bidxdep =  abs(t(Qa[,idxdep])%*%B[ieq]) >= tolf ;
      if (any( bidxdep )){  
        flag = -2;
      } else {
        st=qr(t(A[ieq,]))
        Qat=qr.Q(st)
        Rat=qr.R(st)
        idel=st$pivot[idxdep]
        numDepend = length(idel)
        A=A[-idel,]
        B=B[-idel,,drop=F]
        neqtl = neqtl - numDepend;
        nctl = nctl - numDepend;
        ieq = 1:neqtl;
        xieq=xieq[-idel,] 
        indxact=indxact[-(1:numDepend) ]
        indxact = indxact - numDepend;      
        Act = A[indxact,,drop=F]
        iact = iact - numDepend;       
      }    
    }
    
    if (iact >= numVar){
      iact = max(neqtl, numVar-1);
      indxact = indxact[1:iact];
      Act = A[indxact,,drop=F];
      xieq = matrix(0,nctl,1);
      xieq[indxact] = 1;
    }
    
    if (iact > neqtl){
      sa=qr(t(A[ieq,]))
      Qat=qr.Q(sa)
      Rat=qr.R(sa)
      idxdep =which(abs(diag(Ra))<tolf)
      if (length(idxdep)>0){       
        idel2=st$pivot[idxdep]
        idelEq   = idel2[which(idel2 <= neqtl)]
        idelIneq = idel2[which(idel2 > neqtl)]
        if (length(idelEq)>0){
          if(neqtl){
            indxact = 1:neqtl
          }else{indxact=NULL}
        } else{
          indxact=indxact[-idelIneq]
        }
        xieq = matrix(0,nctl,1);
        xieq[indxact] = 1;
        Act = A[indxact,,drop=F]
        iact = length(indxact); 
      }
    }
    
    s=qr(t(Act))
    Q=qr.Q(s,complete=T)
    R=qr.R(s,complete=T)
    Z = Q[,(iact+1):numVar]
    
    if(flag != -2 & iact > 0){       
      deltax = Q[,1:iact]%*%(solve(t(R[1:iact,1:iact]), (B[indxact] - Act%*%X)))
      X = X + deltax
      err = A%*%X - B;
      err[ieq] = abs(err[ieq])
      if (any(err > .Machine$double.eps)){
        Xb = ginv(Act)%*%B[indxact,,drop=F];
        errb = A%*%Xb - B;
        errb[ieq] = abs(errb[ieq])
        if (max(errb) < max(err)){
          X = Xb
        }
      }
    }
    
    if (length(idel)>0){
      indall=indall[-idel]
      Astan = Astan[indall,,drop=F]
    }
    
    if (flag == -2){   
      indxact = indall[indxact]
      flag = -2
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
      
    }
    
    err = 0
    if (neqtl >= numVar){
      err = max(abs(A[ieq,]%*%X-B[ieq]))
      if(err > 1e-8){
        flag = -2
        indxact = indall[indxact]
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        
        
      } else {   
        if(max(A%*%X-B) > 1e-8){
          flag = -2
        }     
      }
      
      if(qpflag){
        lambdact = -ginv(R)%*%(t(Q)%*%(H%*%X+f))
      } else{
        lambdact = -ginv(R)%*%(t(Q)%*%f)
      }
      
      lambda[indall[indxact],,drop=F] = (lambdact/Astan[indxact]);
      indxact = indall[indxact];
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))    
    }      
  } else {
    if(iact == 0){
      Q = diag(1,numVar)
      R = NULL
      Z = 1
    }else{
      s=qr(t(Act))
      Q=qr.Q(s,complete=T)
      R=qr.R(s,complete=T)
      Z = Q[,(iact+1):numVar,drop=F]
    }   
  }
  
  const = A%*%X-B;
  constdf = 0
  tolconnorm = .Machine$double.eps
  
  if(nctl > neqtl){
    constdf=max(const[(neqtl+1):nctl]) 
    id=which.max(const[(neqtl+1):nctl])
    tolconnorm = tolcon/Astan[neqtl+id]
  }
  
  if (constdf > tolconnorm){   
    if(neqtl){
      indxact2 = 1:neqtl
    }else{indxact2=NULL}
    A2=cbind(rbind(A,matrix(0,1,numVar)),rbind(matrix(0,neqtl,1),-matrix(1,nctl+1-neqtl,1)))
    st=solqp(NULL,matrix(c(rep(0,numVar),1),ncol=1),A2,rbind(B,1e-5),rbind(X,constdf+1),
             neqtl,nrow(A2),numVar+1)
    Xp=st$X;lambdaS=st$lambda;
    slack=Xp[numVar+1]
    X = Xp[1:numVar,,drop=F]
    const = A%*%X-B
    constdf=max(const[(neqtl+1):nctl]) 
    id=which.max(const[(neqtl+1):nctl])
    tolconnorm = .Machine$double.eps
    if(slack > tolconnorm){
      if(slack>1e-8){
        flag=-2
      } else{
        flag=-4
      }
      lambda[indall,,drop=F] = (lambdaS[(1:nctl),,drop=F]/Astan);
      if(neqtl){
        indxact = 1:neqtl
      }else{indxact=NULL}
      indxact = indall[indxact]
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))     
    } else {
      if(neqtl){
        indxact = 1:neqtl
      }else{indxact=NULL}
      Act = A[indxact,,drop=F]
      iact = length(indxact)
      xieq = matrix(0,nctl,1)
      xieq[indxact] = 1
      if(iact == 0){
        Q = matrix(0,numVar,numVar)
        R=NULL
        Z=1
      } else{
        s=qr(t(Act))
        Q=qr.Q(s,complete=T)
        R=qr.R(s,complete=T)
        Z = Q[,(iact+1):numVar]
      }   
    }
  } 
  
  if (iact >= numVar - 1) {
    iterflag = 1
  } 
  
  if(qpflag){
    ggf=H%*%X+f 
    ss=getDir(Z,H,ggf,numVar,f)
    Dir=ss$Dir
    searchDir=getElement(ss,"searchDir")
  }else{
    ggf=f
    if(length(Z)>1){
      Dir=-Z%*%t(Z)%*%ggf
    } else{
      Dir=-Z*Z*ggf  
    }
    searchDir='SteepDescent'
    if(norm(Dir)<1e-10 & neqtl){
      lambdact = -ginv(R)%*%(t(Q)%*%ggf)
      lambda[indall[indxact],,drop=F] = (lambdact/Astan[indxact]);
      indxact = indall[indxact];
      return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))    
    }
  }
  
  
  ind_old = 0
  while (iter < imax){
    iter = iter + 1
    GDir=A%*%Dir
    indf = which((GDir > tolerr * norm(Dir,"2"))  &  !as.matrix(xieq))
    if(length(indf)==0){
      minalpha = 1e16
      ind=NULL
    }else{
      dist = abs(const[indf,,drop=F])/GDir[indf,,drop=F];
      minalpha = min(dist);
      ind2 = which(dist == minalpha);
      ind = indf[min(ind2)]
    }
    
    flagremove = 0
    if(length(indf)>0 & is.finite(minalpha)){
      if(identical(searchDir,'Newton')){
        if(minalpha>1){
          minalpha=1
          flagremove = 1
        }
        X = X+minalpha*Dir
      }else{ 
        X = X+minalpha*Dir        
      } 
    } else { 
      if(identical(searchDir,'Newton')){
        minalpha = 1
        X = X + Dir
        flagremove = 1
        
      } else{
        if(!qpflag | identical(searchDir,'NegCurv')){
          if(norm(Dir,"2")>tolerr){
            minalpha = 1e16
            X = X + minalpha*Dir
            flag = -3                   
          } else { 
            flag = -7 
          }
          indxact = indall(indxact)
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
          
        }else{
          if(qpflag){
            ZHZ = t(Z)%*%H%*%Z; 
            Zggf = t(Z)%*%ggf;
            Dirproj = ginv(ZHZ)%*%(-Zggf)
          }
          
          if (qpflag & (norm(ZHZ%*%Dirproj+Zggf,"2") > 10*.Machine$double.eps*(norm(ZHZ,"2") + norm(Zggf,"2")))){
            if(norm(Dir,"2") > tolerr){
              minalpha = 1e16
              X = X + minalpha*Dir
              flag = -3
            } else {
              flag = -7
            }
            indxact = indall[indxact]
            return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
            
          } else {
            Dir = Z%*%Dirproj   
            if(t(ggf)%*%Dir > 0){
              Dir = -Dir
            }
            searchDir = 'singular'
            GDir=A%*%Dir
            indf = which((GDir > tolerr * norm(Dir,"2"))  &  !as.matrix(xieq))
            if(length(indf)==0){
              minalpha = 1e16
              ind=NULL
            }else{
              dist = abs(const[indf,,drop=F])/GDir[indf,,drop=F];
              minalpha = min(dist);
              ind2 = which(dist == minalpha);
              ind = indf[min(ind2)]
            }
            if (minalpha > 1){
              minalpha = 1;
              flagremove = 1
            }
            X = X + minalpha*Dir           
          }     
        }      
      } 
    } 
    
    if(qpflag){
      ggf=H%*%X+f
    }
    
    const = A%*%X-B;
    const[ieq] = abs(const[ieq,,drop=F]);
    if(neqtl<nctl){
      constdf = max(const[(neqtl+1):nctl])
    } else{constdf=0}
    
    if (flagremove){
      if (iact>0){
        rlambda = -ginv(R)%*%(t(Q)%*%ggf)
        lambdact = rlambda;
        lambdact[ieq] = abs(rlambda[ieq,,drop=F]);
        indlam = which(lambdact < 0)
        if(!length(indlam)){
          lambda[indall[indxact]] = (rlambda/Astan[indxact]);
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
          
        }
        
        minindact = which(indxact == min(indxact[indlam]))[1]
        if(!is.na(minindact)){
          Act=Act[-minindact,,drop=F] 
          xieq[indxact[minindact]] = 0;
          
          st=qprm(Q,R,minindact);
          Q=st$Q
          R=st$R
          indxact=indxact[-minindact]
          iact = length(indxact)
          iterflag = 0;
          ind = 0;
        }
      }else {
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))   
      } 
      flagremove = 0
    }
    
    if (constdf > 1e5*tolerr){
      flag = -1
    }else{ 
      flag = 1;
    } 
    
    if(!is.null(ind) ){
      xieq[ind,]=1;
      indplus = length(indxact) + 1
      Act=rbind(Act,A[ind,])
      indxact=c(indxact,ind)
      m=nrow(Act)
      n=ncol(Act)
      
      st = qrin(Q,R,indplus,t(A[ind,,drop=F]))
      Q=st$Q
      R=st$R
      iact = length(indxact)
    }
    
    if(!iterflag){
      m=nrow(Act)
      n=ncol(Act)
      Z=Q[,(m+1):n,drop=F]
      if(iact == numVar - 1){
        iterflag = 1
      }
      ind_old = 0
    } else {
      rlambda = -ginv(R)%*%(t(Q)%*%ggf)
      if (is.infinite(rlambda[1]) && rlambda[1] < 0){
        rlambda = -ginv(t(Act + sqrt(.Machine$double.eps)*matrix(rnorm(m*n),nrow=m)))%*%ggf  
      }
      lambdact = rlambda;
      lambdact[ieq]=abs(lambdact[ieq,,drop=F]);
      indlam =which(lambdact<0)
      if(length(indlam)){ 
        if (minalpha > tolerr){
          minl=min(lambdact)
          minindact=which(lambdact == min(lambdact))[1]
        }else{
          minindact = which(indxact == min(indxact[indlam,,drop=F]))[1]
        }
        Act=Act[-minindact,,drop=F]
        xieq[indxact[minindact]] = 0
        st=qprm(Q,R,minindact);
        Q=st$Q
        R=st$R
        Z = Q[,numVar,drop=F]
        ind_old = indxact[minindact]
        indxact=indxact[-minindact]
        iact = length(indxact)
        
      }else{ 
        lambda[indall[indxact]]= (rlambda/Astan[indxact,,drop=F])
        indxact = indall[indxact]
        return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
      }
    }
    
    if(qpflag){
      Zggf = t(Z)%*%ggf
      if (length(Zggf)>0 & (norm(Zggf,"2") < 1e-15)){
        Dir = matrix(0,numVar,1)
        searchDir = 'ZeroStep'
      }else{
        ss=getDir(Z,H,ggf,numVar,f)
        Dir=ss$Dir
        searchDir=getElement(ss,"searchDir")
      }  
    }else{
      if(!iterflag){
        Dir=-Z%*%(t(Z)%*%ggf)
        gradDir = norm(Dir,'2')
      }else{
        gradDir = t(Z)%*%ggf
        if(gradDir>0){
          Dir=-Z
        }else{
          Dir = Z
        }
      } 
      
      if(abs(gradDir) < 1e-10){
        if(!ind_old){
          rlambda = -ginv(R)%*%(t(Q)%*%ggf)
          indxacttmp = indxact; Qtmp = Q; Rtmp = R;
        }else{ 
          indxacttmp = indxact
          indxacttmp[(minindact+1):(iact+1)] = indxact[minindact:iact];
          indxacttmp[minindact] = ind_old;
          st = qrin(Q,R,minindact,t(A[ind_old,]))
          Qtmp=st$Q
          Rtmp=st$R
        } 
        lambdact = rlambda;
        if(neqtl){
          lambdact[1:neqtl] = abs(lambdact[1:neqtl])
        }
        indlam = which(lambdact < tolerr);
        lambda[indall[indxacttmp]] = (rlambda/Astan[indxacttmp,,drop=F]);
        if(!length(indlam)){
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        }
        
        
        indplusmax = length(indlam);
        indpluscnt = 0;
        m = length(indxacttmp)
        while ((abs(gradDir) < 1e-10) && (indpluscnt < indplusmax)){
          indpluscnt = indpluscnt + 1;
          minindact = indlam[indpluscnt]
          st=qprm(Qtmp,Rtmp,minindact);
          Q=st$Q
          R=st$R
          Z = Q[,m:numVar]
          if(m!=numVar){
            if(length(Z)>1){
              Dir=-Z%*%t(Z)%*%ggf
            } else{
              Dir=-Z*Z*ggf  
            }
            gradDir = norm(Dir,'2')
          }else{
            gradDir = t(Z)%*%ggf
            if(gradDir > 0){
              Dir = -Z
            }else{
              Dir = Z
            }
          }     
        }
        if(abs(gradDir) < 1e-10){
          indxact = indall[indxact]
          return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
        }else{ 
          indxact = indxacttmp
          indxact=indxact[-minindact]
          xieq = matrix(0,nctl,1)
          iact = length(indxact)
          Act = A(indxact,,drop=F)
        } 
        lambda = matrix(0,orignctl,1)
      }
    } 
  } 
  if(iter >= imax){
    indxact = indall[indxact]
    flag = 0
  }
  return(list(X=X,lambda=lambda,iter=iter,indxact=indxact))
  
}

getDir=function(Z,H,ggf,numVar,f){
  if(length(Z)==1){flag=1}else{flag=0}
  if(flag==1){
    ZHZ = Z*H*Z
  }else{
    ZHZ = t(Z)%*%H%*%Z
  }
  
  if(any(eigen(ZHZ)$values>0)){
    tep=chol(ZHZ)
    if(flag==1){
      Dir = - Z*(ginv(tep)%*%( ginv(t(tep))%*%(Z*ggf)))
    }else{
      Dir = - Z%*%(ginv(tep)%*%( ginv(t(tep))%*%(t(Z)%*%ggf)))
    }
    searchDir = "Newton"
    
  } else {
    ss=cholfac(ZHZ)
    sneg=ss$sneg
    if((!is.null(sneg)) & t(sneg)%*%ZHZ%*%sneg < -sqrt(.Machine$double.eps)) {#negative enough
      if(flag==1){
        Dir=Z*sneg
      }else{
        Dir=Z%*%sneg
      }
      searchDir = 'NegCurv'
    }else{
      if(flag==1){
        Dir = - Z*(Z*ggf)
      }else{
        Dir = - Z%*%(t(Z)%*%ggf)
      }
      searchDir = 'SteepDescent'
    }
  }
  return(list(Dir=Dir,searchDir=searchDir))
} 

cholfac=function(ZHZ) {
  sneg=NULL
  n=nrow(ZHZ)
  L=diag(n)
  tol=0
  for(k in 1:(n-1)){
    if(ZHZ[k,k]<=tol){
      elem = matrix(0,length(ZHZ),1)
      elem[k,1] = 1
      sneg = t(L) %*% ginv(elem)
      return(list(L=L,sneg=sneg))
    }else{
      L[k,k] = sqrt(ZHZ[k,k])
      s = (k+1):n;
      L[s,k] = ZHZ[s,k]/L[k,k]
      ZHZ[(k+1):n,(k+1):n]=ZHZ[(k+1):n,(k+1):n]-lower.tri(L[(k+1):n,k]%*%t(L[(k+1):n,k]))
    }
  }
  if(ZHZ[n,n] <= tol){
    elem = matrix(0,length(ZHZ),1)
    elem[n,1] = 1
    sneg = t(L) %*% ginv(elem)
  }else{
    L[n,n] = sqrt(ZHZ[n,n])
  }
  return(list(L=L,sneg=sneg))
} 


qprm=function(Q,R,j){ 
  R=R[,-j,drop=F]
  m=nrow(R)
  n=ncol(R)
  
  for (k in j:min(n,m-1)){
    p=k:(k+1)
    teppl=rotateplan(R[p,k])
    R[p,k]=teppl$x
    G=teppl$G
    if(k<n){
      R[p,((k+1):n)]=G%*%R[p,((k+1):n)]
    }
    Q[,p]=Q[,p]%*%t(G)  
  }
  return(list(Q=Q,R=R)) 
}

rotateplan=function(x){
  if (x[2]!=0){
    r=norm(x,"2")
    G=rbind(x,cbind(-x[2],x[1]))/r
    x=rbind(r,0)
  }else{
    G=diag(length(x))
  }
  return(list(G=G,x=x))
}

qrin=function(Q,R,j,x){
  x=as.matrix(x)
  if(length(R)>0){
    m=nrow(R)
    n=ncol(R)
  }else{
    m=0;n=0;
  }
  
  if(n==0){
    ss=qr(x)
    Q=qr.Q(ss,complete=T)
    R=qr.R(ss,complete=T)
    return(list(Q=Q,R=R))
  }
  R=cbind(R,matrix(0,nrow(R),1))
  if(j<n+1){
    R[,(j+1):(n+1)] = R[,j:n]
    R[,j] = t(Q)%*%x
  }else if(j==n+1){
    R[,j] = t(Q)%*%x
  }
  n=n+1
  
  if(m>j){
    for (k in (m-1):j){
      p=k:(k+1)
      teppl=rotateplan(R[p,j])
      R[p,j]=teppl$x
      G=teppl$G
      if(k<n){
        R[p,((k+1):n)]=G%*%R[p,((k+1):n)]
      }
      Q[,p]=Q[,p]%*%t(G) 
    }
  }
  return(list(Q=Q,R=R)) 
}






