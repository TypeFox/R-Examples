

################################
#####  Cox model: local    #####
#####  Net, Enet & Lasso   #####
#####  Coordinate Descent  #####
################################


##########################################################
#####  Cox with Regularization (coordinate descent)  #####
##########################################################
###  lambda1*||_1+lambda2/2*||_2 <-> lambda*(alpha_1+(1-alpha)/2_2)

loCoxnet=function(x, y, w, w0=NULL, h=NULL, hnext=NULL, Omega=NULL, penalty=c("Lasso","Enet", "Net"), alpha=1, lambda=NULL, nlambda=50, rlambda=NULL, nfolds=1, foldid=NULL, adaptive=c(FALSE,TRUE), aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-6, thresh2=1e-10, maxit=1e+5) {
  
  #fcall=match.call()
  penalty=match.arg(penalty)
  if (penalty=="Lasso") {
    penalty="Enet"
    alpha=1.0
  }
  
  if (penalty=="Net" & is.null(Omega)) {
    penalty="Enet"
    cat("Enet was performed as no input of Omega")
  }
  
  if (is.null(w0)) {
    w0=seq(min(w), max(w), length=10)
  }
  
  if (is.null(hnext)) {
    hnext=h+0.01
  }
  hnext=c(h, h+hnext); hnext=unique(hnext)
  
  fit=switch(penalty,
             "Enet"=locoxEnet(x,y,w,w0,h,hnext,alpha,lambda,nlambda,rlambda,nfolds,foldid,adaptive[1],aini,isd,keep.beta,thresh,thresh2,maxit),
             "Net"=locoxNet(x,y,w,w0,h,hnext,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,adaptive,aini,isd,keep.beta,thresh,thresh2,maxit))
  
  #fit$call=fcall
  class(fit)=c("loCoxnet", "Coxnet")
  return(fit)
}



##########################
#####  Enet (Lasso)  #####
##########################

locoxEnet=function(x, y, w, w0, h, hnext=h+0.001, alpha=1, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, adaptive=FALSE, aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-5, thresh2=1e-10, maxit=1e+5){
  
  penalty=ifelse(alpha==1,"Lasso","Enet")
  
  N0=nrow(x);p=ncol(x);nw0=length(w0)
  
  ### scaleC and standardized
  tem=scaleC(x)
  xscale=tem$sd; x=tem$x
  rm(tem)
  
  
  ###  Full data  ###
  prep0=list()
  for (iw0 in 1:nw0) {
    prep0[[iw0]]=locoxprep(x,y,w,w0[iw0],h)
  }
  
  ### Adaptive based on Ridge (L2)  
  if (adaptive) {
    if (is.null(aini)) 
      aini=locoxini(x,y,w,w0,h)
    wbeta=aini$wbeta
    rm(aini)
  } else {
    wbeta=matrix(1,ncol=nw0,nrow=p)
  }
  
  
  ### Lambda path
  if (is.null(lambda) ){
    tem_lambda_max=numeric(nw0)
    for (iw0 in 1:nw0) {
      tem_lambda_max[iw0]=max_loclambdaC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,alpha,wbeta[,iw0],N0)
    }
    lambda_max=max(tem_lambda_max);temm=min(tem_lambda_max)
    lambda_min=ifelse(is.null(rlambda),ifelse(N0>=p,temm*0.0001,temm*0.05),temm*rlambda)
    lambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    nlambda=length(lambda)
  }
  
  
  #####  Run  #####
  out=list();thresh2i=thresh2;ithresh2=0;
  repeat{
    nlambdai=nlambda;lambdai=lambda;ithresh2=ithresh2+1
    for(iw0 in 1:nw0){
      out[[iw0]]=locoxenetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,thresh2i,maxit)
      #out[[iw0]]$Beta=out[[iw0]]$Beta/xscale
      nlambdai=out[[iw0]]$nlambda;lambdai=lambda[1:nlambdai]
      if(nlambdai==0) break
      #out[[iw0]]$nlambda=NULL
    }
    thresh2i=thresh2i/10
    if(all(sapply(out,function(x){x$nlambda>0})))break
    if(ithresh2==10)stop("Need larger lambda!")
  }
  
  ### re-fit using thresh2=0
  out=list(); flag=numeric(nw0)
  for(iw0 in 1:nw0){
    out[[iw0]]=locoxenetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,0,maxit)
    if (!isd) out[[iw0]]$Beta=out[[iw0]]$Beta/xscale
    flag[iw0]=out[[iw0]]$flag
  }
  
  if (nfolds==1 & is.null(foldid)) {
    if (nw0==1) {
      Beta=out[[1]]$Beta
      fit=data.frame(lambda=lambdai, nzero=apply(Beta!=0, 2, sum))
      return(list(Beta=Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=flag))
    } else {
      Beta=list()
      for(il in 1:nlambdai){Beta[[il]]=sapply(out,function(x){x$Beta[,il]})}
      fit=data.frame(lambda=lambdai, nzero=sapply(Beta, function(x){sum(apply(x!=0,1,any))}))
      return(list(Beta=Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=flag))
    }
  } else {
    
    ###  Split data  ###
    if(is.null(foldid)){
      foldid=coxsplitw(w, nfolds)
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid);N0i=numeric(nfolds)
    for(i in 1:nfolds){N0i[i]=sum(tb[-i])}
    prepk=list()
    for(iw0 in 1:nw0){
      prepki=list()
      for(i in 1:nfolds){
        temid=which(foldid!=i)
        prepki[[i]]=locoxprep(x[temid,],y[temid,],w[temid],w0[iw0],h)
      }
      prepk[[iw0]]=prepki
      if(any(sapply(prepki,function(x){is.null(x)}))){stop("Larger bandwidth")}
    }
    
    ###  Cross-validation PL  ###
    cvPLs=list()
    for(i in 1:nfolds){
      cvPLi=matrix(0,nrow=nlambdai,ncol=nw0)
      for(iw0 in 1:nw0){
        outi=cvlocoxenetC(prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$Kh1,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
        cvPLi[1:nlambdai,iw0]=outi$lf-outi$ll
        if(outi$nlambda==0){stop("Need larger lambda!")}
      }
      cvPLs[[i]]=cvPLi
    }
    
    cvm=numeric(nlambdai);cvse=numeric(nlambdai)
    for(i in 1:nlambdai){
      temcv=matrix(sapply(cvPLs,function(x){x[i,]}),nrow=nw0)
      temcv=apply(temcv,2,mean)
      cvm[i]=sum(temcv)/N0;cvse[i]=sd(temcv)*sqrt(nfolds/N0^2)
    }
    indexm=which.max(cvm)
    
    if(nlambdai<nlambda & indexm==nlambdai){
      nlambdai=nlambdai+1
      repeat{
        for(iw0 in 1:nw0){
          outsi=locoxenetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambda[nlambdai],1,wbeta[,iw0],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,0,maxit)
          outsi$Beta=outsi$Beta/xscale
          if(outsi$nlambda==0){stop("Need larger lambda!")}
          
          out[[iw0]]$Beta=cBind(out[[iw0]]$Beta,outsi$Beta)
          out[[iw0]]$ll=c(out[[iw0]]$ll,outsi$ll)
          out[[iw0]]$nzero=c(out[[iw0]]$nzero,outsi$nzero)
        }
        
        cvPLsi=NULL
        for(i in 1:nfolds){
          cvPLi=matrix(0,nrow=1,ncol=nw0)
          for(iw0 in 1:nw0){
            outi=cvlocoxenetC(prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$tevent,alpha,lambda[nlambdai],1,wbeta[,iw0],prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$Kh1,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
            cvPLi[1,iw0]=outi$lf-outi$ll
            if(outi$nlambda==0){stop("Need larger lambda!")}
          }
          cvPLsi=rbind(cvPLsi,cvPLi)
        }
        temcv=apply(cvPLsi,1,mean)
        cvm=c(cvm,sum(temcv)/N0);cvse=c(cvse,sd(temcv)*sqrt(nfolds/N0^2))
        
        indexm=which.max(cvm)
        if(nlambdai==nlambda | indexm<nlambdai)break
        nlambdai=nlambdai+1
      }
    }
    
    indexi=which.max(cvm)
    indexj=which(cvm>(cvm-cvse)[indexi])[1]
    temi=rep("",nlambdai);temi[indexi]="max"; # temi[indexj]="se"
    temCV=data.frame(lambda=lambda[1:nlambdai],cvm=cvm,cvse=cvse,index=temi,stringsAsFactors=FALSE)
    
    ###  Cross-validation PL for bandwidth based on max  ###
    temH=NULL;lambdao=lambda[indexi]
    if(!is.null(hnext)){
      nh=length(hnext);cvm=numeric(nh);cvse=numeric(nh)
      for(ih in 1:nh){
        preph=list()
        for(iw0 in 1:nw0){
          prephi=list()
          for(i in 1:nfolds){
            temid=which(foldid!=i)
            prephi[[i]]=locoxprep(x[temid,],y[temid,],w[temid],w0[iw0],hnext[ih])
          }
          preph[[iw0]]=prephi
        }
        
        cvPLh=matrix(0,nrow=nfolds,ncol=nw0)
        for(i in 1:nfolds){
          cvPLi=numeric(nw0)
          for(iw0 in 1:nw0){
            outj=locoxenetC(preph[[iw0]][[i]]$x,preph[[iw0]][[i]]$tevent,alpha,lambdao,1,wbeta[,iw0],preph[[iw0]][[i]]$Kh,preph[[iw0]][[i]]$Kh1,preph[[iw0]][[i]]$N,preph[[iw0]][[i]]$nevent,preph[[iw0]][[i]]$nevent1,preph[[iw0]][[i]]$loc1,preph[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit)
            betai=matrix(outj$Beta,ncol=1)
            lfi=loclbetaC(betai,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
            lli=loclbetaC(betai,prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n)
            cvPLi[iw0]=lfi-lli
          }
          cvPLh[i,]=cvPLi
        }
        temcv=apply(cvPLh,1,mean)
        cvm[ih]=sum(temcv)/N0;cvse[ih]=sd(temcv)*sqrt(nfolds/N0^2)
      }
      temH=data.frame(h=hnext[1:ih],cvm=cvm[1:ih],cvse=cvse[1:ih])
    }
    
    
    if (nw0==1) {
      
      Beta=out[[1]]$Beta
      temCV=data.frame(lambda=temCV$lambda, cvm=temCV$cvm, cvse=temCV$cvse, nzero=apply(Beta!=0, 2, sum), index=temCV$index)
      return(list(Beta=Beta, fit=temCV, lambda.max=lambdai[indexi], cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      
    } else {
      
      Beta=list()
      for(il in 1:nlambdai){Beta[[il]]=sapply(out,function(x){x$Beta[,il]})}
      temCV=data.frame(lambda=temCV$lambda, cvm=temCV$cvm, cvse=temCV$cvse, nzero=sapply(Beta, function(x){sum(apply(x!=0,1,any))}), index=temCV$index)
      
      if (!keep.beta) {
        return(list(Beta=Beta[[indexi]],fit=temCV,lambda.max=lambdai[indexi],cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      } else {
        return(list(Beta=Beta,fit=temCV,lambda.max=lambdai[indexi],cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      } 
    }
    
  }
}



#############################
#####  Network (local)  #####
#############################

locoxNet=function(x, y, w, w0, h, hnext=h+0.001, Omega=NULL, alpha=1, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, adaptive=c(FALSE,TRUE), aini=NULL, isd=FALSE, keep.beta=FALSE, thresh=1e-5, thresh2=1e-10, maxit=1e+5){
  
  penalty=ifelse(alpha==1,"Lasso","Net")
  
  N0=nrow(x);p=ncol(x);nw0=length(w0)
  
  ### centered and standardized
  tem=scaleC(x)
  xscale=tem$sd; x=tem$x
  rm(tem)
  
  ###  Full data  ###
  prep0=list()
  for(iw0 in 1:nw0){
    prep0[[iw0]]=locoxprep(x,y,w,w0[iw0],h)
  } 
  
  ### Adaptive based on Ridge (L2)
  if (any(adaptive)>0) {
    if (is.null(aini)) 
      aini=locoxini(x,y,w,w0,h)
    if (adaptive[1] & !adaptive[2]) {
      wbeta=aini$wbeta
      sgn=matrix(1,ncol=nw0,nrow=p)
    } else if (!adaptive[1] & adaptive[2]) {
      wbeta=matrix(1,ncol=nw0,nrow=p)
      sgn=aini$sgn
    } else {
      wbeta=aini$wbeta
      sgn=aini$sgn
    }
    rm(aini)
  } else {
    wbeta=matrix(1,ncol=nw0,nrow=p)
    sgn=matrix(1,ncol=nw0,nrow=p)
  }
  
  
  ### Lambda path
  if(is.null(lambda)){
    tem_lambda_max=numeric(nw0)
    for(iw0 in 1:nw0){
      tem_lambda_max[iw0]=max_loclambdaC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,alpha,wbeta[,iw0],N0)
    }
    lambda_max=max(tem_lambda_max);temm=min(tem_lambda_max)
    lambda_min=ifelse(is.null(rlambda),ifelse(N0>=p,temm*0.0001,temm*0.05),temm*rlambda)
    lambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
  }else{
    nlambda=length(lambda)
  }
  
  ### Correlation/Adjacent matrix
  diag(Omega)=0; Omegas=list()
  Adjm=Omega
  if(any(Omega!=0)){
    nadj=apply(Adjm,2,function(x){sum(x>0)}) # NO of adj
    loc=matrix(0,ncol=p,nrow=max(nadj)) # position of adj
    ndegree=rowSums(Adjm) # d=\sum w or d=#{Adjm>0}=nadj
    
    for(iw0 in 1:nw0){
      Omegas[[iw0]]=Omega
      for(i in 1:p){
        if(nadj[i]>0){
          loc[1:nadj[i],i]=which(Adjm[,i]>0)
          Omegas[[iw0]][loc[,i],i]=Omega[loc[,i],i]/sqrt(ndegree[i]*ndegree[loc[,i]])*sgn[i,iw0]*sgn[loc[,i],iw0]
        }
      }
    }
  }else{loc=matrix(0,ncol=p,nrow=1);nadj=numeric(p);ndegree=rep(0,p)}
  
  si=numeric(p)
  if(any(ndegree>0)){indexi=which(ndegree>0);si[indexi]=1/sqrt(ndegree[indexi])}
  Ls=list();L=-Omega;diag(L)=ndegree
  for(iw0 in 1:nw0){sj=si*sgn[,iw0];Si=diag(sj);Ls[[iw0]]=Si%*%L%*%Si}
  
  
  #####  Run  #####
  out=list();thresh2i=thresh2;ithresh2=0;
  repeat{
    nlambdai=nlambda;lambdai=lambda;ithresh2=ithresh2+1
    for(iw0 in 1:nw0){
      out[[iw0]]=locoxnetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,thresh2i,maxit)
      #out[[iw0]]$Beta=out[[iw0]]$Beta/xscale
      nlambdai=out[[iw0]]$nlambda;lambdai=lambda[1:nlambdai]
      if(nlambdai==0)break
    }
    thresh2i=thresh2i/10
    if(all(sapply(out,function(x){x$nlambda>0})))break
    if(ithresh2==10)stop("Need larger lambda!")
  }
  
  ### re-fit using thresh2=0
  out=list(); flag=numeric(nw0)
  for(iw0 in 1:nw0){
    out[[iw0]]=locoxnetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,0,maxit)
    if (!isd) out[[iw0]]$Beta=out[[iw0]]$Beta/xscale
    flag[iw0]=out[[iw0]]$flag
  }
  
  if (nfolds==1 & is.null(foldid)) {
    if (nw0==1) {
      Beta=out[[1]]$Beta
      fit=data.frame(lambda=lambdai, nzero=apply(Beta!=0, 2, sum))
      return(list(Beta=Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=flag))
    } else {
      Beta=list()
      for(il in 1:nlambdai){Beta[[il]]=sapply(out,function(x){x$Beta[,il]})}
      fit=data.frame(lambda=lambdai, nzero=sapply(Beta, function(x){sum(apply(x!=0,1,any))}))
      return(list(Beta=Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=flag))
    }
  } else {
    
    ###  Split data  ###
    if (is.null(foldid)) {
      foldid=coxsplitw(w, nfolds)
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid);N0i=numeric(nfolds)
    for(i in 1:nfolds){N0i[i]=sum(tb[-i])}
    prepk=list()
    for(iw0 in 1:nw0){
      prepki=list()
      for(i in 1:nfolds){
        temid=which(foldid!=i)
        prepki[[i]]=locoxprep(x[temid,],y[temid,],w[temid],w0[iw0],h)
      }
      prepk[[iw0]]=prepki
      if(any(sapply(prepki,function(x){is.null(x)}))){stop("Larger bandwidth")}
    }
    
    ###  Cross-validation PL  ###
    cvPLs=list()
    for(i in 1:nfolds){
      cvPLi=matrix(0,nrow=nlambdai,ncol=nw0)
      for(iw0 in 1:nw0){
        outi=cvlocoxnetC(prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$tevent,alpha,lambdai,nlambdai,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$Kh1,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
        cvPLi[1:nlambdai,iw0]=outi$lf-outi$ll
        if(outi$nlambda==0){stop("Need larger lambda!")}
      }
      cvPLs[[i]]=cvPLi
    }
    cvm=numeric(nlambdai);cvse=numeric(nlambdai)
    for(i in 1:nlambdai){
      temcv=matrix(sapply(cvPLs,function(x){x[i,]}),nrow=nw0)
      temcv=apply(temcv,2,mean)
      cvm[i]=sum(temcv)/N0;cvse[i]=sd(temcv)*sqrt(nfolds/N0^2)
    }
    indexm=which.max(cvm)
    
    if(nlambdai<nlambda & indexm==nlambdai){
      nlambdai=nlambdai+1
      repeat{
        
        for(iw0 in 1:nw0){
          outsi=locoxnetC(prep0[[iw0]]$x,prep0[[iw0]]$tevent,alpha,lambda[nlambdai],1,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],prep0[[iw0]]$Kh,prep0[[iw0]]$Kh1,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n,p,N0,thresh,0,maxit)
          outsi$Beta=outsi$Beta/xscale
          if(outsi$nlambda==0){stop("Need larger lambda!")}
          
          out[[iw0]]$Beta=cBind(out[[iw0]]$Beta,outsi$Beta)
          out[[iw0]]$ll=c(out[[iw0]]$ll,outsi$ll)
          out[[iw0]]$nzero=c(out[[iw0]]$nzero,outsi$nzero)
        }
        
        cvPLsi=NULL
        for(i in 1:nfolds){
          cvPLi=matrix(0,nrow=1,ncol=nw0)
          for(iw0 in 1:nw0){
            outi=cvlocoxnetC(prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$tevent,alpha,lambda[nlambdai],1,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$Kh1,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
            cvPLi[1,iw0]=outi$lf-outi$ll
            if(outi$nlambda==0){stop("Need larger lambda!")}
          }
          cvPLsi=rbind(cvPLsi,cvPLi)
        }
        temcv=apply(cvPLsi,1,mean)
        cvm=c(cvm,sum(temcv)/N0);cvse=c(cvse,sd(temcv)*sqrt(nfolds/N0^2))
        
        indexm=which.max(cvm)
        if(nlambdai==nlambda | indexm<nlambdai)break
        nlambdai=nlambdai+1
      }
    }
    
    indexi=which.max(cvm)
    indexj=which(cvm>(cvm-cvse)[indexi])[1]
    temi=rep("",nlambdai);temi[indexi]="max"; # temi[indexj]="se"
    temCV=data.frame(lambda=lambda[1:nlambdai],cvm=cvm,cvse=cvse,index=temi,stringsAsFactors=FALSE)
    
    ###  Cross-validation PL for bandwidth based on max  ###
    temH=NULL;lambdao=lambda[indexi]
    if(!is.null(hnext)){
      nh=length(hnext);cvm=numeric(nh);cvse=numeric(nh)
      for(ih in 1:nh){
        preph=list()
        for(iw0 in 1:nw0){
          prephi=list()
          for(i in 1:nfolds){
            temid=which(foldid!=i)
            prephi[[i]]=locoxprep(x[temid,],y[temid,],w[temid],w0[iw0],hnext[ih])
          }
          preph[[iw0]]=prephi
        }
        
        cvPLh=matrix(0,nrow=nfolds,ncol=nw0)
        for(i in 1:nfolds){
          cvPLi=numeric(nw0)
          for(iw0 in 1:nw0){
            outj=locoxnetC(preph[[iw0]][[i]]$x,preph[[iw0]][[i]]$tevent,alpha,lambdao,1,wbeta[,iw0],Ls[[iw0]],Omegas[[iw0]],preph[[iw0]][[i]]$Kh,preph[[iw0]][[i]]$Kh1,preph[[iw0]][[i]]$N,preph[[iw0]][[i]]$nevent,preph[[iw0]][[i]]$nevent1,preph[[iw0]][[i]]$loc1,preph[[iw0]][[i]]$n,p,N0i[i],thresh,0,maxit)
            betai=matrix(outj$Beta,ncol=1)
            lfi=loclbetaC(betai,prep0[[iw0]]$x,prep0[[iw0]]$Kh,prep0[[iw0]]$N,prep0[[iw0]]$nevent,prep0[[iw0]]$nevent1,prep0[[iw0]]$loc1,prep0[[iw0]]$n)
            lli=loclbetaC(betai,prepk[[iw0]][[i]]$x,prepk[[iw0]][[i]]$Kh,prepk[[iw0]][[i]]$N,prepk[[iw0]][[i]]$nevent,prepk[[iw0]][[i]]$nevent1,prepk[[iw0]][[i]]$loc1,prepk[[iw0]][[i]]$n)
            cvPLi[iw0]=lfi-lli
          }
          cvPLh[i,]=cvPLi
        }
        temcv=apply(cvPLh,1,mean)
        cvm[ih]=sum(temcv)/N0;cvse[ih]=sd(temcv)*sqrt(nfolds/N0^2)
      }
      temH=data.frame(h=hnext[1:ih],cvm=cvm[1:ih],cvse=cvse[1:ih])
    }
    
    if (nw0==1) {
      
      Beta=out[[1]]$Beta
      temCV=data.frame(lambda=temCV$lambda, cvm=temCV$cvm, cvse=temCV$cvse, nzero=apply(Beta!=0, 2, sum), index=temCV$index)
      return(list(Beta=Beta, fit=temCV, lambda.max=lambdai[indexi], cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      
    } else {
      
      Beta=list()
      for(il in 1:nlambdai){Beta[[il]]=sapply(out,function(x){x$Beta[,il]})}
      temCV=data.frame(lambda=temCV$lambda, cvm=temCV$cvm, cvse=temCV$cvse, nzero=sapply(Beta, function(x){sum(apply(x!=0,1,any))}), index=temCV$index)
      
      if (!keep.beta) {
        return(list(Beta=Beta[[indexi]],fit=temCV,lambda.max=lambdai[indexi],cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      } else {
        return(list(Beta=Beta,fit=temCV,lambda.max=lambdai[indexi],cvh=temH, penalty=penalty, adaptive=adaptive, flag=flag))
      } 
    }
    
  }
}




