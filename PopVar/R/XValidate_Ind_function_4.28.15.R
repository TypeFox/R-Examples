## G.CV must have individuals as rows and markers as columns... so 100 individuals with 500 markers would be a 100x500 matrix
## y.CV must be a numeric vector representing a continuous variable

XValidate.Ind <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, nFold.CV=NULL, nFold.CV.reps=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
  gc(verbose = F) ## Close unused connections
  #con.path <- getwd() ## BGLR will write temp files to the wd
  
  non.BGLR <- models.CV[models.CV %in% c("rrBLUP")]
  BGLR <- models.CV[models.CV %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")]
  if(nFold.CV >= length(y.CV)) stop("nFold too large given the TP size")
  
  for(j in 1:nFold.CV.reps){
    for(i in 1:nFold.CV){
      
      if(i==1){
        rrBLUP.cv <- c()
        BayesA.cv <- c()
        BayesB.cv <- c()
        BayesC.cv <- c()
        BL.cv <- c()
        BRR.cv <- c()
        
        for.VP.mat <- 1:length(y.CV)
        while(length(for.VP.mat) %% nFold.CV != 0) for.VP.mat <- c(for.VP.mat, NA)
        VP.mat <- matrix(for.VP.mat[order(sample(1:length(for.VP.mat), replace = F))], nrow = nFold.CV, byrow = T)
      }
      
      VP.sample <- as.numeric(VP.mat[i, ]); VP.sample <- VP.sample[which(!is.na(VP.sample))]
      TP.sample <- setdiff(1:length(y.CV), VP.sample)
      
      TP.G <- G.CV[TP.sample, ]
      TP.y <- y.CV[TP.sample]
      
      VP.G <- G.CV[VP.sample, ]
      VP.y <- y.CV[VP.sample]
      
      ### non-BGLR models below
      if("rrBLUP" %in% non.BGLR) {RR.pred <- rrBLUP::kinship.BLUP(y=TP.y, G.train=TP.G, G.pred=VP.G, K.method="RR"); rrBLUP.cv[i] <- cor(VP.y, RR.pred$g.pred, use="pairwise.complete.obs")}
      
      ### models from BGLR package
      ## Bayesian A
      if("BayesA" %in% BGLR){
        tryCatch({
          BayesA.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesA")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesA.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesA.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesA.fit)){
          mkr.effs <- as.numeric(BayesA.fit$ETA[[1]]$b); BayesA.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesA.cv[i] <- NA}
      }
      
      
      ## Bayesian B 
      if("BayesB" %in% BGLR){
        tryCatch({
          BayesB.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesB")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesB.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesB.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesB.fit)){
          mkr.effs <- as.numeric(BayesB.fit$ETA[[1]]$b)
          BayesB.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesB.cv[i] <- NA}
      }
      
      
      ## Bayesian C
      if("BayesC" %in% BGLR){
        tryCatch({
          BayesC.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesC")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesC.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesC.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesC.fit)){
          mkr.effs <- as.numeric(BayesC.fit$ETA[[1]]$b)
          BayesC.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesC.cv[i] <- NA}
      }
      
      
      ## Bayesian LASSO
      if("BL" %in% BGLR){
        tryCatch({
          BL.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BL")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BL.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BL.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BL.fit)){
          mkr.effs <- as.numeric(BL.fit$ETA[[1]]$b)
          BL.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BL.cv[i] <- NA}        
      }
          
      
      ### Bayesian Ridge Reg.
      if("BRR" %in% BGLR){
        tryCatch({
          BRR.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BRR")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BRR.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BRR.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BRR.fit)){
          mkr.effs <- as.numeric(BRR.fit$ETA[[1]]$b)
          BRR.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BRR.cv[i] <- NA}
      }
      
    }; gc(verbose=F) ## End of i-fold iteration
    
    if(j==1) cv.list <- paste(models.CV, "cv", sep=".")
    
    for(l in 1:length(cv.list)){ ## Tried using sapply() across cv.list using the "get" function, but did not work within PopVar function
      toget <- cv.list[l]
      if(l == 1) cvs <- get(toget)
      if(l > 1) cvs <- cbind(cvs, get(toget))
    }
    
    if(j == 1) cvs.all <- cvs
    if(j > 1) cvs.all <- rbind(cvs.all, cvs)

  }
  
  if(length(models.CV) == 1){
    bad.models <- F
    if(length(which(is.na(cvs.all))) > 0.025*(nFold.CV*nFold.CV.reps)) bad.models <- T
  }
  
  if(length(models.CV) > 1){
    bad.models <- apply(cvs.all, 2, function(X){if(length(which(is.na(X))) > 0.025*(nFold.CV*nFold.CV.reps)){ ## if more than 2.5% of the iterations resulted in an error 
      return(T)
    }else{return(F)}
    })
  }
  
  if(length(models.CV) > 1) {CV.results <- data.frame(Model=models.CV, r_avg=apply(cvs.all, 2, mean, na.rm=T), r_sd=apply(cvs.all, 2, sd, na.rm=T)) ; rownames(CV.results) <- NULL}
  if(length(models.CV) == 1) {CV.results <- data.frame(Model=models.CV, r_avg=mean(cvs.all), r_sd=sd(cvs.all)) ; rownames(CV.results) <- NULL}
    
  if(length(which(bad.models)) > 0){
    if(length(models.CV) == 1 | (length(which(bad.models)) == length(models.CV))) stop("All model(s) tested was/were removed due to excessive negative values of nu being returned by BGLR::BGLR")
    CV.results <- CV.results[-which(bad.models), ]
    warning(paste("Model(s)", models.CV[which(bad.models)], "was/were removed due to excessive negative values of nu being returned by BGLR::BGLR."))
  }
  
  #CV.lists <- as.data.frame(t(rbind(as.character(models.CV), matrix(c(rrBLUP.cv, BayesA.cv, BayesB.cv, BayesC.cv, BL.cv, BRR.cv), ncol=length(models.CV)))))
  return(list(CV.summary = CV.results))#, iter.CV = CV.lists))
  
} ## End of XValidate function

