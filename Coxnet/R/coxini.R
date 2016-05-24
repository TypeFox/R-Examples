

#################################################
#####  Initial values for adaptive methods  #####
#################################################

coxini=function(x, y){
  
  N0=nrow(x);p=ncol(x)
  nalambda=10
  
  prep0=coxprep(x, y); wbeta=rep(1, p)
  ### Lambda path
  lambda_max=max_lambdaC(prep0$x, prep0$tevent, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, 1, wbeta, N0)
  lambda_min=ifelse(N0>=p, lambda_max*0.0001, lambda_max*0.01)
  alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  
  repeat {
    outi=coxEnet(x, y, alpha=0.0, lambda=alambda, isd=FALSE, ifast=TRUE, keep.beta=TRUE)
    if(!is.null(outi))break
    alambda=alambda*2.0
  }
  
  indexi=ncol(outi$Beta)
  beta0=outi$Beta[, indexi]
  wbeta=1/abs(beta0);sgn=sign(beta0[1:p])
  return(list(wbeta=wbeta, sgn=sgn, lambda=alambda[indexi]))
}



###############
###  local  ###
###############

locoxini=function(x, y, w, w0, h){
  # alambda=NULL; nalambda=10; rlambda=NULL; isd=TRUE; thresh3=1e-15
  N0=nrow(x); p=ncol(x); nw0=length(w0);
  nalambda=10; thresh3=1e-15
    
  tem_lambda_max=numeric(nw0)
  for(iw0 in 1:nw0){
    prep0=locoxprep(x, y, w, w0[iw0], h)
    tem_lambda_max[iw0]=max_loclambdaC(prep0$x, prep0$tevent, prep0$Kh, prep0$Kh1, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, 1.0, rep(1.0, p), N0)
  }
  
  lambda_max=max(tem_lambda_max); temm=min(tem_lambda_max)
  lambda_min=ifelse(N0>=p, temm*0.0001, temm*0.05)
  alambda=lambda_max*(lambda_min/lambda_max)^(c(0:(nalambda-1))/(nalambda-1))
  
  outi=list(); alambdaj=alambda
  repeat {
    nalambdai=nalambda; alambdai=alambdaj
    for(iw0 in 1:nw0){
      outi[[iw0]]=locoxEnet(x, y, w, w0[iw0], h, alpha=0.0, lambda=alambdai, isd=FALSE, thresh2=thresh3)
      nalambdai=length(outi[[iw0]]$fit$lambda);alambdai=alambdaj[1:nalambdai]
      if(nalambdai==0)break
    }
    #thresh3i=thresh3i/10
    alambdaj=alambdaj*2.0
    if(all(sapply(outi, function(x){length(x$fit$lambda)>0})))break
    #if(ithresh3==10)stop("Need larger lambda!")
  }
  
  ### re-fit using thresh2=0
  outi=list()
  for(iw0 in 1:nw0){
    outi[[iw0]]=locoxEnet(x, y, w, w0[iw0], h, alpha=0.0, lambda=alambdai, keep.beta=TRUE, isd=FALSE, thresh2=0)
  }
  
  beta0=sapply(outi, function(x){x$Beta[, nalambdai]})[1:p, ]
  wbeta=matrix(1/abs(beta0), ncol=nw0);sgn=matrix(sign(beta0), ncol=nw0)
  
  return(list(wbeta=wbeta, sgn=sgn))
}




