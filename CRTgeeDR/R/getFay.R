### Calculate the Fay adjusted estimator of the  variance
#Adapted from the function 'gee.var.fg' included in the package 'geesmv',
#authored by Ming Wang
#under the GPL-2 license.
getFay = function(formula,id,family,data,corstr,b,beta,alpha,scale,Y,X,hessMAT,
                  X.t,X.c, B, off, 
                  R.alpha.inv, phi,
                  InvLinkDeriv, InvLink, VarFun,
                  StdErr, dInvLinkdEta,
                  BlockDiag, sqrtW,W,included,typeweights,pi.a,
                  nameTRT,propensity.score,om.t,om.c,
                  nameY,nameMISS,print.log){

  
  data<-as.data.frame(cbind(id,Y))
  names(data)<-c("id","response")
  mat<-as.data.frame(X)
  mat.c<-as.data.frame(X.c)
  mat.t<-as.data.frame(X.t)
  
  
  ### Fit the GEE model to get the estimate of parameters \hat{\beta};
  #library(stats)
  #gee.fit <- gee(formula,data=data,id=id,family=family,corstr=corstr)
  beta_est <- beta
  alpha <- alpha
  len <- length(beta_est)
  len_vec <- len^2
  
  ### Estimate the robust variance for \hat{\beta}
  #data$id <- id
  cluster<-cluster.size(data$id)
  ncluster<-max(cluster$n)
  size<-cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  mat.t$subj <- rep(unique(data$id), cluster$n)
  mat.c$subj <- rep(unique(data$id), cluster$n)
  
  if(is.character(corstr)){
    var <- switch(corstr,
                  "independence"=cormax.ind(ncluster),
                  "exchangeable"=cormax.exch(ncluster, alpha),
                  "AR-M"=cormax.ar1(ncluster, alpha))
  }else{
    print(corstr)
    stop("'working correlation structure' not recognized")
  }   
  if(is.character(family)){
    family <- switch(family,
                     "gaussian"="gaussian",
                     "binomial"="binomial",
                     "poisson"="poisson")
  }else{ 
    if(is.function(family)){
      family <- family()[[1]]
    }else{
      print(family)
      stop("'family' not recognized")
    }    
  }

  cov.beta<-unstr<-matrix(0,nrow=len,ncol=len)
  step11<-hessMAT
  
  step12<-matrix(0,nrow=len,ncol=len)
  step13<-matrix(0,nrow=len_vec,ncol=1)
  step14<-matrix(0,nrow=len_vec,ncol=len_vec)
  p<-matrix(0,nrow=len_vec,ncol=size)
  for (i in 1:size){
    y<-as.matrix(data$response[data$id==unique(data$id)[i]])
    covariate<-as.matrix(subset(mat[,-length(mat[1,])], mat$subj==unique(data$id)[i]))
    covariate.t<-as.matrix(subset(mat.t[,-length(mat.t[1,])], mat.t$subj==unique(data$id)[i]))
    covariate.c<-as.matrix(subset(mat.c[,-length(mat.c[1,])], mat.c$subj==unique(data$id)[i]))
    ncluster=cluster$n[i]
    var1=var[1:ncluster,1:ncluster]
    hessi<-hessianCommunity(Y=y, X=covariate,X.t=covariate.t,X.c=covariate.c,
                            B=B[which(data$id==unique(data$id)[i]),], beta=beta, off=off[which(data$id==unique(data$id)[i])], InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun,
                            R.alpha.inv=R.alpha.inv[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])], StdErr=StdErr[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])],
                            dInvLinkdEta=dInvLinkdEta[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])], 
                            sqrtW[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])],W[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])],
                            included[which(data$id==unique(data$id)[i]),which(data$id==unique(data$id)[i])],typeweights=typeweights,pi.a=pi.a)
    xx<-hessi$hess
    Qi <- xx%*%solve(step11)
    Ai<-diag((1-pmin(b,diag(Qi)))^(-0.5))
    xy<-Ai%*%hessi$esteq
    step12<-step12+xy%*%t(xy)
    step13<-step13+vec(xy%*%t(xy))
    p[,i]<-vec(xy%*%t(xy))  
  }
  for (i in 1:size){
    dif<-(p[,i]-step13/size)%*%t(p[,i]-step13/size)
    step14<-step14+dif
  }
  cov.beta<-solve(step11)%*%(step12)%*%solve(step11)
  cov.var<-size/(size-1)*kronecker(solve(step11), solve(step11))%*%step14%*%kronecker(solve(step11), solve(step11))
  return(list(sandadjfay=cov.beta, cov.var=cov.var))
}

hessianCommunity = function(Y, X,X.t,X.c, B, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, sqrtW,W,included,typeweights,pi.a){
    beta.new<-beta
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
          
        }else{
          hess <- (1-pi.a)*crossprod( StdErr.c %*% dInvLinkdEta.c %*%X.c, R.alpha.inv  %*% StdErr.c %*%dInvLinkdEta.c %*% X.c)+
            (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*%X.t, R.alpha.inv %*% StdErr.t %*%dInvLinkdEta.t %*% X.t)
          esteq <- crossprod(  StdErr %*%dInvLinkdEta %*%X ,R.alpha.inv%*% W  %*% StdErr %*% as.matrix(Y - B[,"Bi"])) +
            (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% as.matrix(B[,"B.c"]-mu.c))+
            (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% as.matrix(B[,"B.t"]-mu.t))  
        }
      }             
    }
    
  return(list(beta = beta.new, hess = hess, esteq=esteq))
}


