#Adapted from the function 'updateBeta' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
### Method to update coefficients.  Goes to a maximum of 10 iterations, or when
### rough convergence has been obtained.
updateBeta = function(Y, X,X.t,X.c, B, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, sqrtW,W,included,typeweights,pi.a){
  beta.new <- beta
  conv=F
  for(i in 1:10){
    eta <- as.vector(X%*%beta.new) + off
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)	
    diag(StdErr) <- sqrt(1/VarFun(mu))
   
    if(is.null(B)){  
      if(is.null(typeweights)){
        hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*%X, R.alpha.inv %*% sqrtW %*% StdErr %*%dInvLinkdEta %*% X)
        esteq <- crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% StdErr %*% as.matrix(Y - mu))          
      }else{
        if(typeweights=="GENMOD"){
          hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*%X, R.alpha.inv %*% sqrtW %*% StdErr %*%dInvLinkdEta %*% X)
          esteq <- crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% StdErr %*% as.matrix(Y - mu))          
          } else{
          hess <- crossprod(  StdErr %*% dInvLinkdEta %*%X, R.alpha.inv %*% W %*% StdErr %*%dInvLinkdEta %*% X)
          esteq <- crossprod(   StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% W %*% StdErr %*% as.matrix(Y - mu))   
        }
      }

    } else{
            
      nn<-length(Y)
      StdErr.c <- Diagonal(nn)
      dInvLinkdEta.c <- Diagonal(nn)
      eta.c <- as.vector(X.c%*%beta.new) + off  	
      diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
      mu.c <- InvLink(eta.c)	    
      diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
      
      nn<-length(Y)
      StdErr.t <- Diagonal(nn)
      dInvLinkdEta.t <- Diagonal(nn)
      eta.t <- as.vector(X.t%*%beta.new) + off
      diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
      mu.t <- InvLink(eta.t)      
      diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
      
      if(is.null(typeweights)){
        hess <- (1-pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*%X.c, R.alpha.inv  %*% StdErr.c %*%dInvLinkdEta.c %*% X.c)+
          (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*%X.t, R.alpha.inv %*% StdErr.t %*%dInvLinkdEta.t %*% X.t)
        
        esteq <- crossprod( StdErr %*%dInvLinkdEta %*%X , R.alpha.inv  %*% StdErr %*% as.matrix(Y - B[,"Bi"])) +
          (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% as.matrix(B[,"B.c"]-mu.c))+
          (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% as.matrix(B[,"B.t"]-mu.t))    
        
      }else{
        if(typeweights=="GENMOD"){
          
          hess <- (1-pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*%X.c, R.alpha.inv  %*% StdErr.c %*%dInvLinkdEta.c %*% X.c)+
            (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*%X.t, R.alpha.inv %*% StdErr.t %*%dInvLinkdEta.t %*% X.t)
          
          esteq <- crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% StdErr %*% as.matrix(Y - B[,"Bi"])) +
            (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% as.matrix(B[,"B.c"]-mu.c))+
            (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% as.matrix(B[,"B.t"]-mu.t))    
          #print(hess)  
          #print(esteq)
        }else{

          hess <- (1-pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*%X.c, R.alpha.inv  %*% StdErr.c %*%dInvLinkdEta.c %*% X.c)+
            (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*%X.t, R.alpha.inv %*% StdErr.t %*%dInvLinkdEta.t %*% X.t)
      
          esteq <- crossprod(  StdErr %*%dInvLinkdEta %*%X ,R.alpha.inv%*%  StdErr %*% W %*% as.matrix(Y - B[,"Bi"])) +
            (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% as.matrix(B[,"B.c"]-mu.c))+
            (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% as.matrix(B[,"B.t"]-mu.t))  
          #print(hess)  
          #print(esteq)
          }
      }             
    }
  
    update <- solve(hess, esteq)
    #print(update)
    if(max(abs(update)/beta.new) < 100*tol){break}
    beta.new <- beta.new + as.vector(update)
    #print(beta.new)
  }
  return(list(beta = beta.new, hess = hess, esteq=esteq))
}

