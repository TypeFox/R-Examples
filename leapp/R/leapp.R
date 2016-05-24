leapp <-
function(data,pred.prim,pred.covar=NULL,
         O = NULL, num.fac = "buja", method = "hard", sparse = TRUE, centered = FALSE, verbose = FALSE, perm.num = 50, TOL = 1e-4, length.out = 50){

  # Data transformation 
  N = nrow(data)
  n = ncol(data)
  pred.prim = matrix(pred.prim, n,1)

  if(sum(pred.prim) !=0)
    pred.prim = pred.prim - mean(pred.prim)

  if (!centered)
    data = t(apply(data,1,scale, center = TRUE, scale = FALSE))
  
  if (num.fac == "buja")
    num.fac = num.sv(data,cbind(pred.prim,pred.covar),B = perm.num)

    
  if (is.null(O))
     O = qr.Q( qr(pred.prim), complete = TRUE)
  
  if ((O%*%pred.prim)[1]<0)
    O = -O
  
  if (!is.null(pred.covar)){
    s = ncol(pred.covar)
    pred.covar = O%*%pred.covar
    pred.covar.v1 = pred.covar[1,]
    pred.covar.rest = pred.covar[-1,]
  }else{
    s = 0
    pred.covar.rest = NULL
  }

  data = data%*%t(O)
  Y = data[,1]
  result.svd = AlternateSVD(data[,-1],pred=pred.covar.rest,r=num.fac)  
  sigma.all = result.svd$sigma
  X = result.svd$uest
  if(is.null(pred.covar)){
     Y = Y/sigma.all
  }else{
     Y = (Y - result.svd$coef%*%t(pred.covar.v1))/sigma.all
  }
  if (!is.null(X)){
    X = X/sigma.all
    H = X%*%solve(t(X)%*%X)%*%t(X)
  }else{
    H = NULL
  }
  
  
  if (!sparse){
    # ridge type regression
    result = ridge(X,Y,H,sigma=(n - s - num.fac - 1)/(n - s - num.fac - 3))
    
  }else{
  # IPOD algorithm
    result = IPOD(X,Y,H, method, TOL = TOL, length.out = length.out )
  }
  RetList = list(p = result$p, vest = result.svd$vest,uest = result.svd$uest, gamma = as.vector(sigma.all)*result$gamma, sigma = sigma.all )
  if(verbose){
  print(paste("number of factors:", num.fac))
  if(!is.null(X)){
  png("outlierplot.png")
  plot(X[,1], Y)
  dev.off()
  }
  }
  return(RetList)
}


ridge <-
  function(X,Y,H,sigma){
    N = length(Y)
    if(!is.null(X)){
      res = Y - H%*%Y
    }else{
      res = Y
    }
    sigmahat = sd(res)
    if(sigmahat > 1){
        lambda = 1/(sigmahat -1)
        gamma = 1/(1+lambda)*res
    }else{
        gamma = rep(0,N)
    }
    p = 2*pnorm(-abs(res))
    RetList = list(p = p, gamma = gamma)
    return(RetList)
  }
      

IPOD <-
  function(X,Y,H,
          method ="hard", TOL = 1e-4, length.out = 50){
  if(is.null(X)){
    r = 0
  }else{
    r = ncol(X)
  }
  N = length(Y)
  
  ress = NULL
  gammas = NULL
  if (is.null(X)){
       betaInit = NULL
       lambdas = seq(round(norm(matrix(Y,N,1),'I')/1+1),0, by = -0.1)
   #     lambdas = seq(0.1,2,length.out = length.out)*sqrt(2*log(N))
  }else{
       betaInit = rlm(Y~X-1)$coefficients
       tmp = t(diag(ncol(H))-H)%*%Y/sqrt(1-diag(H))
  
      # lambdas = seq(0.1,2,length.out = length.out)*sqrt(2*log(N))    
       lambdas = seq( round(norm(tmp,'I')/1+1) ,0, by = -0.1) 
  }
  for (sigma in lambdas){
      sigma = sigma/sqrt(2*log(N))
      result = IPODFUN(X,Y,H,sigma,betaInit, method=method, TOL = TOL)
      gammas = cbind(gammas,result$gamma)
      ress = cbind(ress, result$ress)
  }

 
  DF = colSums(abs(gammas)>1e-5)
  
  if(!is.null(X)){
    
     Q = qr.Q(qr(X), complete = TRUE)
     X.unscaled = t(Q[,(r+1):N])
     Y.eff = X.unscaled%*%Y
 
    sigmaSqEsts = colSums((Y.eff%*%matrix(rep(1, ncol(gammas)), 1, ncol(gammas)) - X.unscaled%*%gammas)^2)/(length(Y.eff)-DF)
  }else{
    sigmaSqEsts = colSums((Y%*%matrix(rep(1, ncol(gammas)), 1, ncol(gammas)) - gammas)^2)/(length(Y)-DF)
  }

  
   sigmaSqEsts[sigmaSqEsts < 0] = 0;
  
  if(!all(DF <= N-r)){
    gammas = gammas[, DF<=N-r]
    sigmaSqEsts = sigmaSqEsts[DF <= N-r]
    lambdas = lambdas[DF <= N-r]
    ress = ress[,DF<=N - r]
    DF = DF[DF <= N-r]
  }

  # remove redundancy 
  redundantInds = which(sum(abs(gammas[,-1]-gammas[,1:(ncol(gammas)-1)]))<1e-3)+1
  if (length(redundantInds)>0){
  gammas=gammas[,-redundantInds]
  lambdas=lambdas[-redundantInds]
  sigmaSqEsts= sigmaSqEsts[-redundantInds]
  ress = ress[,-redundantInds]
  DF=DF[-redundantInds]

}
  redundantInds = which(abs(sigmaSqEsts[-1]-sigmaSqEsts[1:(length(sigmaSqEsts)-1)])<1e-3)+1
 if (length(redundantInds)>0){
  gammas=gammas[,-redundantInds]
  lambdas=lambdas[-redundantInds]
  sigmaSqEsts= sigmaSqEsts[-redundantInds]
  ress = ress[,-redundantInds]
  DF=DF[-redundantInds]
 
  }
  
  if(!all(DF <= N/2)){

    gammas = gammas[,DF<=N/2]
    sigmaSqEsts = sigmaSqEsts[DF <= N/2]
    lambdas = lambdas[DF <= N/2]
    ress = ress[,DF<=N/2]
    DF = DF[DF<=N/2]
   
  }
  riskEst = ((N-r)*log(sigmaSqEsts*(N-r-DF)/(N-r))+ (log(N-r)+1)*(DF+1))/(N-r)
  optSet = which(riskEst==min(riskEst))
  gammaOptInd = optSet[which(DF[optSet]==min(DF[optSet]))[1]]
  gammaOpt = gammas[,gammaOptInd]
  resOpt = ress[, gammaOptInd]
  tau = mad(ress[gammas[,gammaOptInd]==0,gammaOptInd])  # C function
  resOpt.scale = resOpt/tau
  p = 2*pnorm(-abs(resOpt.scale))

  return(list( p = p, resOpt.scale = resOpt.scale, gamma = gammaOpt))
}


# utility functions  




FindFpr <-
  function(pvalue,ind,topk){
    N = length(pvalue)
    k = length(topk)
    fpr = rep(0,k) 
    for (i in seq(k)){
  
      indest = order(pvalue)[1:topk[i]]
      fpr[i] = (length(indest)-length(intersect(indest,ind)))/(N-length(ind))
    }
    fpr
}


FindTpr <-
  function(pvalue,ind,topk){
    N = length(pvalue)
    k = length(topk)
    tpr = rep(0,k)
  
    for (i in seq(k)){
   
      indest=order(pvalue)[1:topk[i]]
      tpr[i]=length(intersect(indest,ind))/length(ind)
    }

    tpr
}



FindPrec <-
  function(pvalue,ind,topk){
    k=length(topk)
    prec=rep(0,k) 
    for (i in seq(k)){
      indest=order(pvalue)[1:topk[i]]
      prec[i]=length(intersect(indest,ind))/length(indest)
    }
    prec
}


FindRec <-
  function(pvalue,ind,topk){
  k=length(topk)
  rec=rep(0,k) 
  for (i in seq(k)){ 
    indest=order(pvalue)[1:topk[i]]
    rec[i]=length(intersect(indest,ind))/length(ind)
  }
  rec
}

## wilcoxin statistics
FindAUC <-
  function(pvalue,ind){
    auc=0
    N=length(pvalue)
    pp=length(ind)
    for (i in seq(pp)){
        auc=auc+sum(pvalue[ind[i]]<pvalue[-ind])
    }
    return(auc/pp/(N-pp))
}


    

IPODFUN <-
  function(X,Y,H,sigma, betaInit,
           method = "hard", TOL = 1e-4){
    N = length(Y)
    gamma = matrix(rep(0,N),N,1)
    theta = sigma*sqrt(2*log(N))

    if (is.null(betaInit)){
      r = Y
      if (method == "hard"){
         gamma[abs(r)>theta] = r[abs(r)>theta]
      }else if(method == "soft"){
         gamma[r > theta] = r[r > theta] - theta[r > theta]
         gamma[r < -theta] = r[r < -theta] + theta[r < -theta]
      }
    }else{
      gamma.old = gamma
      r0 = Y - H%*%Y
      yEff = Y
      niter =0
      theta = sigma*sqrt(2*log(N))*sqrt(1-diag(H))

      while(1){
        niter = niter +1
        beta0 = ginv(t(X)%*%X)%*%t(X)%*%yEff
        if (niter == 1){         
           r = Y - X%*%betaInit
        }else{
           r = H%*%gamma.old +r0
        }
        if (method == "hard"){
           gamma = matrix(rep(0,N),N,1)
           gamma[abs(r)>theta] = r[abs(r)>theta]
        }else if(method == "soft"){
           gamma = matrix(rep(0,N),N,1)
           gamma[r > theta] = r[r > theta] - theta[r > theta]
           gamma[r < -theta] = r[r < -theta] + theta[r < -theta]
        }
        if(max(abs(gamma - gamma.old)) < TOL){
           break
        }else{
           gamma.old = gamma
        }
   
        yEff = Y - gamma
      }
    }

      RetList = list(gamma = gamma, ress = r)
}
    

       
      
  
AlternateSVD <-
  function(x, r,
           pred = NULL, max.iter = 10,TOL = 1e-4 ){
    N = nrow(x)
    n = ncol(x)

    if(!is.null(pred)){
      R = x%*%(diag(n) - pred%*%solve(t(pred)%*%pred)%*%t(pred))
      coef = x%*%pred%*%solve(t(pred)%*%pred)
      x = R
    }else{
      coef = NULL
    }

    if(r == 0){
      sigma = apply(x, 1, sd) 
      return(list(sigma = sigma, u=coef, v = NULL))
    }

    sigma = rep(1,N)
    iter = 0
    sigma.old = sigma + 1
    while(1){
      if (iter > max.iter | sum(abs(sigma - sigma.old))/sum(abs(sigma.old)) < TOL)
         break
      data = as.vector(1/sigma)*x
      result.svd = fast.svd(data,tol = 0)
      
      u = result.svd$u[,1:r]%*%diag(result.svd$d[1:r],r,r)
      v = result.svd$v[,1:r]
      sigma.old = sigma

      sigma = apply(x - as.vector(sigma)*u%*%t(v),1,sd)
      iter = iter +1 
    }
  
    uest = as.vector(sigma)*u
    uest = cbind(uest)  
 
  return(list(sigma = sigma ,coef = coef,  uest = uest, vest = v))
}

ROCplot <-
  function(fpr,tpr, main, name.method,
           xlim = c(0,0.2),ylim = c(0.4,1), save = TRUE,name.file = NULL){
    A =ncol(fpr)
    if(save){
    png(paste(name.file,".png",sep=""))
  }
    matplot(fpr, tpr,
            type = "b", col = seq(A),
            main = main,
            xlab = "False positive rate", ylab = "True positive rate",
            xlim = xlim, ylim = ylim)
    legend("bottomright", legend = name.method, col = seq(A), lty = 2)
    if(save){
    dev.off()
  }
  }
    
    
Pvalue <-
    function(dat,mod,mod0){

  
    n <- dim(dat)[2]
    N <- dim(dat)[1]
    df1 <- dim(mod)[2]
    if (is.null(mod0)){
      df0=0
    }else{
    df0 <- dim(mod0)[2]
  }
    p <- rep(0, N)
    Id <- diag(n)
    resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
        t(mod))
    if (is.null(mod0)){
      resid0 = dat
    }else{
    resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% 
        t(mod0))
  }
    rss1 <- resid^2 %*% rep(1, n)
    rss0 <- resid0^2 %*% rep(1, n)
    fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
    p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
    coef = dat%*%mod%*%solve(t(mod)%*%mod)
    sdd = sqrt(rowSums(resid^2)/(n-df1))
    
    tstat = coef[,1]/sdd/sqrt(solve(t(mod)%*%mod)[1,1])
    
     
    return(list(tstat = tstat, fstats = fstats,p=p, coef = coef))
  }



