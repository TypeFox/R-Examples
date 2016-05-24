### Calculate the sandwich estimator
#Adapted from the function 'getSandwich' included in the package 'geeM',
#authored by Lee S. McDaniel and Nick Henderson
#under the GPL-2 license.
getSandwich = function(Y, X,X.t,X.c, B, beta,off, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag, sqrtW,W,included,typeweights,pi.a,print.log){

  eta <- as.vector(X%*%beta) + off
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)			
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  if(is.null(B)){  
    if(is.null(typeweights)){
      scoreDiag <- Diagonal(x= Y - mu)
      BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
      numsand <- as.matrix(crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% BlockDiag %*% StdErr %*% R.alpha.inv %*% StdErr %*% dInvLinkdEta %*% X))
    }else{
      scoreDiag <- Diagonal(x= Y - mu)
      BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
      if(typeweights=="GENMOD"){
        numsand <- as.matrix(crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% BlockDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% X))
      }else{
        numsand <- as.matrix(crossprod(  StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% W %*% StdErr %*% BlockDiag %*% StdErr %*% W %*% R.alpha.inv %*%  StdErr %*% dInvLinkdEta %*% X))
      }
    }
  }else{

    
    nn<-length(Y)
    StdErr.c <- Diagonal(nn)
    dInvLinkdEta.c <- Diagonal(nn)
    eta.c <- as.vector(X.c%*%beta) + off    
    diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
    mu.c <- InvLink(eta.c)	    
    diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
    
    nn<-length(Y)
    StdErr.t <- Diagonal(nn)
    dInvLinkdEta.t <- Diagonal(nn)
    eta.t <- as.vector(X.t%*%beta) + off
    diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
    mu.t <- InvLink(eta.t)      
    diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
    
    scoreDiag <- Diagonal(x= Y - B[,"Bi"])
    scoreDiag.t <-Diagonal(x= B[,"B.t"]-mu.t)
    scoreDiag.c <-Diagonal(x= B[,"B.c"]-mu.c)
    
    if(is.null(typeweights)){
      aa <- crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv  %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag %*% StdErr  %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
      bb <- (pi.a**2)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
      cc <- ((1-pi.a)**2)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.c %*% StdErr.c %*% R.alpha.inv  %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
      ab <- pi.a*crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv  %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
      ac <- (1-pi.a)*crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
      ba <- (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag %*% StdErr %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
      bc <- (pi.a*(1-pi.a))*crossprod(  StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
      ca <- ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag %*% StdErr %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
      cb <- ((1-pi.a)*pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.t %*% StdErr.t %*% R.alpha.inv %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
      numsand <- as.matrix((aa+bb+cc+ab+ac+ba+bc+ca+cb))
    }else{ 
     if(typeweights=="GENMOD"){
       aa <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% X)
       bb <- (pi.a**2)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
       cc <- ((1-pi.a)**2)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.c %*% StdErr.c %*% R.alpha.inv  %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
       ab <- pi.a*crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
       ac <- (1-pi.a)*crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
       ba <- (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% X)
       bc <- (pi.a*(1-pi.a))*crossprod(  StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
       ca <- ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag %*% StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% X)
       cb <- ((1-pi.a)*pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.t %*% StdErr.t %*% R.alpha.inv %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
       numsand <- as.matrix((aa+bb+cc+ab+ac+ba+bc+ca+cb))
      } else{
        aa <- crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv  %*% StdErr %*% W %*% scoreDiag %*% BlockDiag %*% scoreDiag %*% W %*% StdErr  %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
        bb <- (pi.a**2)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
        cc <- ((1-pi.a)**2)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.c %*% StdErr.c %*% R.alpha.inv  %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
        ab <- pi.a*crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv  %*% StdErr %*% W %*% scoreDiag %*% BlockDiag %*% scoreDiag.t %*% StdErr.t  %*% R.alpha.inv  %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
        ac <- (1-pi.a)*crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% W %*% scoreDiag %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
        ba <- (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag %*% W %*% StdErr  %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
        bc <- (pi.a*(1-pi.a))*crossprod(  StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t %*% BlockDiag %*% scoreDiag.c %*% StdErr.c  %*% R.alpha.inv %*% StdErr.c %*% dInvLinkdEta.c %*% X.c)
        ca <- ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag %*% W %*% StdErr %*% R.alpha.inv  %*% StdErr %*% dInvLinkdEta %*% X)
        cb <- ((1-pi.a)*pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c %*% BlockDiag %*% scoreDiag.t %*% StdErr.t %*% R.alpha.inv %*% StdErr.t %*% dInvLinkdEta.t %*% X.t)
        numsand <- as.matrix((aa+bb+cc+ab+ac+ba+bc+ca+cb))
        #print("old")
        #print(numsand)
        #SEE<-crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% W %*% scoreDiag)+
        #  (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
        #  ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c) 
        #numsand<-SEE%*%BlockDiag%*%t(SEE)
        #print("new")
        #print(numsand)
      }       
    }
  }
  hessMat<-as.matrix(hessMat)
  sandvar <- t(solve(hessMat, numsand))
  sandvar <- t(solve(t(hessMat), sandvar))
  

  if(print.log){

  print("Sandwich variance")
  print(sandvar)
  }
  
  return(list(sandvar = sandvar, numsand = numsand, hessMat=hessMat))
}