MIPCA<-function(X, ncp = 2, scale = TRUE, method = c("Regularized","EM"), threshold = 1e-04, nboot = 100, method.mi="Boot",
   Lstart=1000, L=100, verbose=FALSE) 
{
  
  Estim_param<-function(X,scale=F,S,method="Regularized",threshold=10^(-6)){
    # PURPOSE: estimate the parameters phi and sigma of a PCA model
    # from an incomplete or a complete data set 
    if(sum(is.na(X))>0){
      # with missing values
      missing <- which(is.na(X))
      impute.data <- imputePCA(X, scale = scale, ncp = S, method = method, 
                               threshold = threshold)$completeObs
      reference <- PCA(impute.data, scale.unit = scale, graph = FALSE, 
                       ncp = S)
      rec <- reconst(reference, S)
      rec.pca <- as.matrix(X)
      rec.pca[missing] <- rec[missing]
      sd_resid <- apply(impute.data, 2, sd)
      resid <- rec.pca - rec
      if (scale){resid <- sweep(resid, 2, sd_resid, FUN = "/")}
      sigma <- sqrt((sum((resid[-missing])^2))/(nrow(X) * ncol(X) - (length(missing) + ncol(X) + S * (nrow(X) - 1 + ncol(X) - S))))
    }else{
      # without missing values
      reference <- PCA(X, scale.unit = scale, graph = FALSE, ncp = S)
      rec <- reconst(reference,ncp=S)
      rec.pca <- as.matrix(X)
      sd_resid <- apply(X, 2, sd)
      resid <- rec.pca - rec
      if (scale){resid <- sweep(resid, 2, sd_resid, FUN = "/")}
      sigma <- sqrt((sum((resid)^2))/(nrow(X) * ncol(X) - (ncol(X) + S * (nrow(X) - 1 + ncol(X) - S))))
    }
    phi<-diag(c((reference$eig[1:S,1]-sigma^2)/reference$eig[1:S,1],rep(0,min(nrow(X)-1,ncol(X))-S)))
    estimation<-list(sigma=sigma,phi=phi)
    return(estimation)
  }
  
  BayesMIPCA<-function(X.na,nboot=100,ncp,L=100,verbose=T,Lstart=1000,scale=F){
    if(verbose){cat("Multiple Imputation using Bayesian PCA using",nboot,"imputed arrays","\n")}
    if(scale){
      sd.init<-apply(X.na,2,FUN=sd,na.rm=T)
      X.na<-sweep(X.na,2,sd.init,FUN="/")
    }
    res.MI<-res.Over<-list()
    
    ################
    # Initialization
    ################
    
    Mean<-colMeans(X.na,na.rm=TRUE)
    X<-sweep(X.na,2,Mean,FUN="-")
    X.tild<-imputePCA(X,ncp=ncp,method="Regularized",scale=F)$fitt
    sigma<-Estim_param(X=X,S=ncp)[["sigma"]]
    
    ###############
    # Imputation
    ###############
    for (l in 1:(Lstart+nboot*L)){
      #Step I
      X[is.na(X.na)]<-X.tild[is.na(X.na)]+rnorm(sum(is.na(X.na)),mean=0,sd=sigma)#imputation aleatoire
      X<-sweep(X,2,Mean,FUN="+")#decentrage
      
      if (l %in% c(Lstart+(1:nboot)*L)){
        res.MI[[(l-Lstart)/L]]<-X
        res.Over[[(l-Lstart)/L]]<-sweep(X.tild+rnorm(nrow(X.na)*ncol(X.na),mean=0,sd=sigma),2,Mean,FUN="+")
        if(verbose){cat(paste((l-Lstart)/L,"...",sep=""))}
      }
      
      #Step P
      Mean<-colMeans(X)
      X<-sweep(X,2,Mean,FUN="-")
      res.pca<-PCA(X,ncp=ncol(X),scale.unit=F,graph=F)
      param<-Estim_param(X=X,S=ncp)
      phi<-param[["phi"]]
      sigma<-param[["sigma"]]
      Moy<-res.pca[["svd"]][["U"]]%*%diag(res.pca[["svd"]][["vs"]])%*%phi%*%t(res.pca[["svd"]][["V"]])
      Var<-diag(rep(sigma^2/min(nrow(X)-1,ncol(X))*sum(phi),ncol(X)))
      X.tild<-Moy+rmvnorm(nrow(X.na),sigma=Var)
    }
    class(res.MI)<-"BayesMIPCA"
    if(scale){
      res.MI<-lapply(res.MI,FUN=function(res){sweep(res,2,sd.init,FUN="*")})
      res.Over<-lapply(res.Over,FUN=function(res){sweep(res,2,sd.init,FUN="*")})
    }
    call<-list(X=X.na,nboot=nboot,ncp=ncp,L=L,verbose=verbose,Lstart=Lstart,scale=scale)
    return(list(res.MI=lapply(res.MI,as.data.frame),res.Over=res.Over,call=call))
  }
  
  res.Over<-list()
  method <- match.arg(method, c("Regularized", "regularized", 
                                "EM", "em"), several.ok = T)[1]
  method <- tolower(method)
  missing <- which(is.na(X))
  inputeData <- imputePCA(X, scale = scale, ncp = ncp, method = method, 
                          threshold = threshold)$completeObs
  if(method.mi%in%c("Bayes","bayes")){
    #BayesPCA
    res<-BayesMIPCA(X.na=X,nboot=nboot,ncp=ncp,L=L,verbose=verbose,Lstart=Lstart,scale=scale)
    res.MI<-res$res.MI
    res.Over<-res$res.Over 
  }else if(method.mi%in%c("Boot","boot")){
    #Bootstrap 
    if(verbose){cat("Multiple Imputation using bootstrap",method,"PCA using",nboot,"imputed arrays","\n")}
    reference <- PCA(inputeData, scale.unit = scale, graph = FALSE, 
                     ncp = ncp)
    rec <- reconst(reference, ncp)
    rec.pca <- as.matrix(X)
    rec.pca[missing] <- rec[missing]
    resid <- rec.pca - rec
    sdResid <- apply(inputeData, 2, sd)
    if (scale) 
      resid <- t(t(resid)/sdResid)
    sigma <- sqrt((sum((resid[-missing])^2))/(nrow(X) * ncol(X) - 
                                                (length(missing) + ncol(X) + ncp * (nrow(X) - 1 + ncol(X) - 
                                                                                      ncp))))
    rownames(rec.pca) <- rownames(X)
    res.MI <- list()
    for (i in 1:nboot) {
      if(verbose){cat(paste(i,"...",sep=""))}
      resid.star <- matrix(rnorm(nrow(X) * ncol(X), 0, sigma), 
                           ncol = ncol(X))
      if (scale) 
        resid.star <- t(t(resid.star) * sdResid)
      resid.star[missing] <- NA
      Xstar <- rec + resid.star - matrix(mean(resid.star, na.rm = TRUE), 
                                         ncol = ncol(resid.star), nrow = nrow(resid.star))
      acpboot <- PCA(imputePCA(Xstar, scale = scale, ncp = ncp, 
                               method = method, threshold = threshold)$completeObs, 
                     scale.unit = scale, ncp = ncp, graph = FALSE)
      residstar2 <- matrix(rnorm(nrow(X) * ncol(X), 0, sigma), 
                           ncol = ncol(X))
      if (scale) 
        residstar2 <- t(t(residstar2) * sdResid)
      rec.pca[missing] <- (reconst(acpboot, ncp) + residstar2)[missing]
      res.Over[[i]]<-reconst(acpboot, ncp) + residstar2
      res.MI[[i]] <- rec.pca
      dimnames(res.MI[[i]]) = list(rownames(X), colnames(X))
      res.MI[[i]]<-as.data.frame(res.MI[[i]])
    }
    
  }else{
    stop("method.mi is misspecified")
  }
  if(verbose){cat("\ndone!")}
  result = list(res.MI = res.MI, res.imputePCA = inputeData,
                call = list(X = X, ncp = ncp, missing = missing, nboot = nboot, 
                            scale = scale,L=L,verbose=verbose,Lstart=Lstart,res.Over=res.Over))
  class(result) <- c("MIPCA", "list")
  return(result)
}