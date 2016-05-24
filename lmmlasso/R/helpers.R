ArmijoRule <- function(xGroup,yGroup,LGroup,b,j,cut,HkOldJ,HkJ,JinNonpen,lambda,weights,nonpen,ll1,ll2,converged,control)
{
 b.new <- b
 bJ <- b[j]
 grad <- -cut + bJ*HkOldJ
 if (JinNonpen) {dk <- -grad/HkJ} else {dk <- MedianValue(grad,HkJ,lambda/weights[j],bJ)}
       
 if (dk!=0)
  { 
     # calculate delta_k
     if (JinNonpen) deltak <- dk*grad + control$gamma*dk^2*HkJ
     else deltak <- dk*grad + control$gamma*dk^2*HkJ + lambda/weights[j]*(abs(bJ+dk)-abs(bJ))

     fctOld <- ObjFunction(xGroup=xGroup,yGroup=yGroup,LGroup=LGroup,b=b,weights=weights,
                    lambda=lambda,nonpen=nonpen,ll1=ll1,ll2=ll2)
     for (l in 0:control$max.armijo)
      { 
         b.new[j] <- bJ + control$a_init*control$delta^l*dk
           
         fctNew <- ObjFunction(xGroup=xGroup,yGroup=yGroup,LGroup=LGroup,b=b.new,weights=weights,
                        lambda=lambda,nonpen=nonpen,ll1=ll1,ll2=ll2)
         addDelta <- control$a_init*control$delta^l*control$rho*deltak
         if (fctNew <= fctOld + addDelta)
           {
             b[j] <- bJ + control$a_init*control$delta^l*dk
             fct <- fctNew
             break
           }
         if (l==control$max.armijo)
           {
             if (trace>2) cat("Armijo for b_",j," not successful","\n")
             converged <- converged + 2
             fct <- fctOld
           }
      }
  } 
 return(list(b=b,fct=fct,converged=converged))
}

HessianMatrix <- function(xGroup,LGroup,activeSet,N,hessian,mat)
{
  for (i in 1:N)
    {  
     mat[i,] <- diag(t(xGroup[[i]][,activeSet,drop=FALSE])%*%LGroup[[i]]%*%xGroup[[i]][,activeSet,drop=FALSE])
    }
  hessian[activeSet] <- apply(mat,2,sum)
  return(hessian)
}

LambdaInv <- function(Z,ZId,Psi,sigma)
{
  LambdaInverse <- solve(sigma^2*ZId + quad.tform(Psi,Z)) 
  return(list(LambdaInverse))
}

MLloglik <- function(xGroup,yGroup,LGroup,b,ntot,N,activeSet)
{
 
  l1 <- l2 <- numeric(N) ; ll2b <- 0
  for (i in 1:N)
   {
    l1[i] <- -determinant(LGroup[[i]])$modulus
    l2[i] <- quad.form(LGroup[[i]],yGroup[[i]]-xGroup[[i]]%*%b)
   }

  ll <- - 1/2*(sum(l1) + sum(l2) + ntot*log(2*pi) + ll2b)
  return(loglik=ll)    
}

MLpdDiag <- function(zGroup,resGroup,zIdGroup,sigma,pars,q,thres=10^(-4),ll1,ll4,trace,CovOpt,VarInt)
{
 LPsiIter <- diag(pars,nrow=length(pars))
 
 for (s in 1:q)
  {
   # optimization with optimize
   if (CovOpt=="optimize")
     {
       optRes <- optimize(MLpdSymFct,interval=VarInt,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                              sigma=sigma,a=s,b=s,LPsi=LPsiIter)
       if (optRes$minimum<thres)
         {
           LPsiIter[s,s] <- 0
           if (trace>2) cat("thres was used","\n")
         } else LPsiIter[s,s] <- optRes$minimum
   # optimization with nlminb
     } else if (CovOpt=="nlminb")
     {  
       optRes <- nlminb(pars[s],MLpdSymFct,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                              sigma=sigma,a=s,b=s,LPsi=LPsiIter)
       if (optRes$par<thres)
         {
           LPsiIter[s,s] <- 0
           if (trace>2) cat("thres was used","\n")
         } else LPsiIter[s,s] <- optRes$par
     }

   ll <- ll1 + optRes$objective + ll4
   if (trace>3)  print(ll) 
  }
  
 pars <- diag(LPsiIter)
 
 return(list(pars=pars,fct=ll))
  
}

MLpdIdent <- function(zGroup,resGroup,zIdGroup,sigma,pars,q,thres=10^(-4),ll1,ll4,trace,CovOpt,VarInt)
{
 LPsiIter <- diag(q)

 # optimization with optimize (no starting value)
 if (CovOpt=="optimize")
  {
   optRes <- optimize(MLpdIdentFct,interval=VarInt,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                          sigma=sigma,LPsi=LPsiIter)
   if (optRes$minimum<thres)
     {
       pars <- 0
       if (trace>2) cat("thres was used","\n")
     } else pars <- optRes$minimum
   
 # optimization with nlminb (starting value provided)
  } else if (CovOpt=="nlminb")
  {
   optRes <- nlminb(pars,MLpdIdentFct,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                          sigma=sigma,LPsi=LPsiIter)
   if (optRes$par<thres)
     {
       pars <- 0
       if (trace>2) cat("thres was used","\n")
     } else pars <- optRes$par
  }

 ll <- ll1 + optRes$objective + ll4
 if (trace>3)  print(ll)
 
 return(list(pars=pars,fct=ll))
  
}

MLpdIdentFct <- function(thetak,zGroup,zIdGroup,resGroup,sigma,LPsi)
{
  Psi <- thetak^2*LPsi
  LambdaGroup <- mapply(MLpdSymLambda,Z=zGroup,ZId=zIdGroup,MoreArgs=list(sigma=sigma,Psi=Psi))

  ll2 <- mapply(MLpdSymObj,LambdaGroup,resGroup)
  1/2*sum(ll2)
}

MLpdSym <- function(zGroup,resGroup,zIdGroup,sigma,pars,q,thres=10^(-4),ll1,ll4,trace,CovOpt,CovInt,VarInt)
{
 LPsiIter <- t(triang(pars,q))

 # optimize over the variance parameters
 for (s in 1:q)
  {
   # optimization with optimize
   if (CovOpt=="optimize")
     {
       optRes <- optimize(MLpdSymFct,interval=VarInt,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                              sigma=sigma,a=s,b=s,LPsi=LPsiIter)
       if (optRes$minimum<thres)
         {
           LPsiIter[s,s] <- 0
           if (trace>2) cat("thres was used","\n")
         } else LPsiIter[s,s] <- optRes$minimum
   # optimization with nlminb
     } else if (CovOpt=="nlminb")
     {  
       optRes <- nlminb(LPsiIter[s,s],MLpdSymFct,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                              sigma=sigma,a=s,b=s,LPsi=LPsiIter)
       if (optRes$par<thres)
         {
           LPsiIter[s,s] <- 0
           if (trace>2) cat("thres was used","\n")
         } else LPsiIter[s,s] <- optRes$par
     }
   
   ll <- ll1 + optRes$objective + ll4
   if (trace>3)  print(ll) 
  }

 # optimize over the covariance parameters

 if (q>1)
  {
   for (l in 2:q)
    {
     for (r in 1:(l-1))
      {
        # optimization with optimize
        if (CovOpt=="optimize")
          {
            optRes <- optimize(MLpdSymFct,interval=CovInt,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                                   sigma=sigma,a=l,b=r,LPsi=LPsiIter)
            if (optRes$minimum<thres)
              {
                LPsiIter[r,l] <- 0
                if (trace>2) cat("thres was used","\n")
              } else LPsiIter[l,r] <- optRes$minimum
        # optimization with nlminb
          } else if (CovOpt=="nlminb")
          {
            optRes <- nlminb(LPsiIter[l,r],MLpdSymFct,zGroup=zGroup,zIdGroup=zIdGroup,resGroup=resGroup,
                                   sigma=sigma,a=l,b=r,LPsi=LPsiIter)
            if (optRes$par<thres)
              {
                LPsiIter[r,l] <- 0
                if (trace>2) cat("thres was used","\n")
              } else LPsiIter[l,r] <- optRes$par
          }     
       ll <- ll1 + optRes$objective + ll4
       if (trace>3)  print(ll) 
      }
    }
  }
   
  # return the parameters
 pars <- vecli(t(LPsiIter))
 
 return(list(pars=pars,fct=ll))
}

MLpdSymFct <- function(thetak,zGroup,zIdGroup,resGroup,sigma,a,b,LPsi)
{
  LPsi[a,b] <- thetak
  Psi <- tcrossprod(LPsi)
  LambdaGroup <- mapply(MLpdSymLambda,Z=zGroup,ZId=zIdGroup,MoreArgs=list(sigma=sigma,Psi=Psi))

  ll2 <- mapply(MLpdSymObj,LambdaGroup,resGroup)
  1/2*sum(ll2)
  
}

MLpdSymLambda <- function(Z,ZId,Psi,sigma)
{
  Lambda <- sigma^2*ZId + quad.tform(Psi,Z)
  return(list(Lambda))
}

MLpdSymObj <- function(Lambda,res)  determinant(Lambda)$modulus + MyQuadFormInv(Lambda,res)

MLsigma <- function(zGroup,zIdGroup,resGroup,q,ll1,ll4,true.sigma,Psi,trace,CovOpt,VarInt)
{
  
  ZPZtGroup <- lapply(zGroup,MLsigmaZPsiZt,Psi=Psi)

  if (CovOpt=="optimize")
   {
  	optRes <- optimize(MLsigmaFct,interval=VarInt,ZPZtGroup=ZPZtGroup,zIdGroup=zIdGroup,resGroup=resGroup)
  	sigma <- optRes$minimum
   } else if (CovOpt=="nlminb")
   {
     optRes <- nlminb(true.sigma,MLsigmaFct,ZPZtGroup=ZPZtGroup,zIdGroup=zIdGroup,resGroup=resGroup)
     sigma <- optRes$par
   }  

  ll <- ll1 + optRes$objective + ll4
  if (trace>3) print(ll)
  
  return(list(sigma=sigma,fct=ll))  
}

MLsigmaLambda <- function(ZId,ZPZt,sigma)
{
  Lambda <- sigma^2*ZId + ZPZt
  return(list(Lambda))
}

MLsigmaFct <- function(sigma,ZPZtGroup,zIdGroup,resGroup)
{
  LambdaGroup <- mapply(MLsigmaLambda,zIdGroup,ZPZtGroup,MoreArgs=list(sigma=sigma))

  ll2 <- mapply(MLpdSymObj,LambdaGroup,resGroup)
  1/2*sum(ll2)
}

MLsigmaZPsiZt <- function(Z,Psi) ZPsiZt <- quad.tform(Psi,Z)

MedianValue <- function(grad,hessian,lambda,bj)
 {
  median(c((lambda-grad)/hessian,-bj,(-lambda-grad)/hessian))
 }

MyQuadFormInv <- function(M, x) crossprod(x, solve(M, x))

ObjFunction <- function(xGroup,yGroup,LGroup,b,weights,lambda,nonpen,resGroup=NULL,ll1=ll1,ll2=NULL)
{
  ResAs <- function(x,y,b,activeSet) y-x[,activeSet,drop=FALSE]%*%b[activeSet,drop=FALSE]
  
  tResLRes <- function(L,res) 1/2*quad.form(L,res) 

  if (missing(resGroup))
    {
     activeSet <- which(b!=0)
     resGroup <- mapply(ResAs,x=xGroup,y=yGroup,MoreArgs=list(b=b,activeSet=activeSet),SIMPLIFY=FALSE)
    }

  if (missing(ll2)) ll2 <- nlogdet(LGroup=LGroup)
  
  ll <- ll1 + ll2 + sum(mapply(tResLRes,LGroup,resGroup)) + lambda*sum(abs(b[-nonpen])/weights[-nonpen])

 return(Fc=ll)      
}

ResAsSplit <- function(x,y,b,f,activeset)
{
  r <- y-x[,activeset,drop=FALSE]%*%b[activeset,drop=FALSE]
  resGroup <- split(r,f)
  return(resGroup)
}

SoftThreshold <- function(z,g)
 {
  sign(z)*max(abs(z)-g,0)
 }

ZIdentity <- function(Z)
{
 ZId <- diag(dim(Z)[[1]])
 return(list(ZId))
}

as1 <- function(xGroup,LGroup,activeSet,N)
  {
    fs <- function(x,l,a) {l%*%x[,a]}
    
    SGroup <- mapply(fs,xGroup,LGroup,MoreArgs=list(a=activeSet),SIMPLIFY=FALSE)

    return(SGroup)
  }

as2 <- function(x,y,b,j,activeSet,group,sGroup)
{
  r <- y-x[,-c(j),drop=FALSE]%*%b[-j]
  rGroup <- split(r,group)

  as3 <- function(s,r,j) crossprod(r,s[,j])

  ma <- mapply(as3,sGroup,rGroup,MoreArgs=list(j=match(j,activeSet)))

  sumMa <- sum(ma)
  
  return(sumMa)
}

covStartingValues <- function(xGroup,yGroup,zGroup,zIdGroup,b,ntot,N,lower=-10,upper=10)
{
  optimize1 <- function(x,y,b) y-x%*%b

  optimize2 <- function(gamma,zId,ZtZ)
    {
      H <- zId + exp(2*gamma)*ZtZ
      return(list(H=H))
    }

  optimize3 <- function(res,zId,ZtZ,gamma)
     {
       lambda <- optimize2(gamma,zId,ZtZ)
       logdetH <- determinant(lambda$H)$modulus
       quadH <- quad.form.inv(lambda$H,res)
       return(c(logdetH,quadH))
     }
  
  optimize4 <- function(gamma)
    {
      optH <- mapply(optimize3,resGroup,zIdGroup,ZtZ=ztzGroup,MoreArgs=list(gamma=gamma))
      H1 <- optH[1,]
      H2 <- optH[2,]

      fn <- ntot*log(sum(H2)) + sum(H1)
      fn
    }
  
  optimize5 <- function(z) tcrossprod(z)

  resGroup <- mapply(optimize1,x=xGroup,y=yGroup,MoreArgs=list(b=b),SIMPLIFY=FALSE)
  ztzGroup <- mapply(optimize5,z=zGroup,SIMPLIFY=FALSE)
  
  optRes <- optimize(f=optimize4,interval=c(lower,upper))

  gamma <- optRes$minimum

  quadH <- mapply(optimize3,resGroup,zIdGroup,ztzGroup,MoreArgs=list(gamma=gamma))[2,]

  sig <- sqrt(1/ntot*sum(quadH))
  tau <- exp(gamma)*sig
  objfct <- 1/2*(optRes$objective + ntot*(1-log(ntot)))

  return(list(tau=tau,sigma=sig,opt=objfct))
}

nlogdet <- function(LGroup)
 {
   nlogdetfun <- function(L)
    {
     -1/2*determinant(L)$modulus[1]
    }
   
   sum(mapply(nlogdetfun,LGroup))
}
