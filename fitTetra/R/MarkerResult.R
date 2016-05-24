MarkerResult <-
function(marker,markername,ratio,model,mutype,sdtype="sd.const",ptype="p.free",
              ng=5,clus=TRUE,mrkresult,nsamp=0,mustart=NA,sdstart=rep(0.075,ng),
              maxiter=40,maxn.bin=200,nbin=200, plothist=FALSE,freq=TRUE) {

  modelresult<-mrkresult$stats
  probslist<-mrkresult$probs
  modelname<-getModelName(mutype,ptype)
  result <- tryCatch({
    cdsr <- CodomMarker(y=ratio,ng=ng, ptype=ptype, mutype=mutype, sdtype=sdtype,
        clus=clus, mu.start=mustart, sd.start=sdstart, maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin,
        plothist=plothist, maintitle=paste(marker,markername,modelname))
  }, error = function(x) {x} )

  modelresult$marker[model] <- marker
  modelresult$markername[model] <- markername
  modelresult$model[model] <- modelname
  modelresult$nsamp[model] <- nsamp #including the NA but excluding the select=F samples
  modelresult$nsel[model] <- length(ratio)
  modelresult$message[model] <- result$message

  if (inherits(result,"error")) {
    modelresult$npar[model] <- NA
    modelresult$iter[model] <- NA
    modelresult$dip[model] <- NA
    modelresult$LL[model] <- NA
    modelresult$AIC[model] <- NA
    modelresult$BIC[model] <- NA
  } else { #success
    modelresult$npar[model] <- result$npar
    modelresult$iter[model] <- result$iter
    modelresult$LL[model] <- result$loglik
    modelresult$AIC[model] <- result$AIC # Akaike's Information Criterion, better than BIC for small models
    modelresult$BIC[model] <- result$BIC # Bayesian Information Criterion, optimal for complex models
                                         # more severe penalty than AIC if nsel>=8
                                         # we use BIC for model selection
    if (ng==5) {
      #the following output will not fit or is not meaningful if ng!=5
      #the following statistics may reflect the quality of the fit but are not used any more:
      mudiff <- result$psi$mu[2:ng] - result$psi$mu[1:(ng-1)]  # the differences between successive mu's
      sigavs <- (result$psi$sigma[2:ng] + result$psi$sigma[1:(ng-1)]) / 2 #the average of successive sigma's
      separation <- mudiff/sigavs  # a measure of the (ng-1) separations between the peaks
      modelresult$minsepar[model] <- min(separation)
      #the following calculates the fraction of genotyped samples given various P thresholds:
      #note that these fractions relate to the total number of ratio values pased to CodomMarker;
      #this includes only the non-NA samples for which select is TRUE (in fitTetra)
      maxPna <- apply(result$post,1,max) #find the max P of each row of p-values (= of each sample)
      maxP <- maxPna[!is.na(maxPna)] # omit missing values
      modelresult$meanP[model] <- mean(maxP) #the average maxP over all samples
      modelresult$P80[model] <- length(maxP[maxP>0.80])/length(maxP) #fraction samples with maxP>0.80
      modelresult$P90[model] <- length(maxP[maxP>0.90])/length(maxP) #idem
      modelresult$P95[model] <- length(maxP[maxP>0.95])/length(maxP) #idem
      modelresult$P975[model] <- length(maxP[maxP>0.975])/length(maxP) #idem
      modelresult$P99[model] <- length(maxP[maxP>0.99])/length(maxP) #idem; this criterion is used in fitTetra
      #the following give the parameters of the fitted model and are used by fitTetra for plotting
      nfirst <- which(names(modelresult)=="mutrans0")
      modelresult[model,nfirst:(nfirst+ng-1)] <- as.numeric(result$psi$mu)
      modelresult[model,(nfirst+ng):(nfirst+2*ng-1)] <- as.numeric(result$psi$sigma)
      modelresult[model,(nfirst+2*ng):(nfirst+3*ng-1)] <- as.numeric(result$psi$p)
      modelresult[model,(nfirst+3*ng):(nfirst+4*ng-1)] <- as.numeric(result$back$mu)
      modelresult[model,(nfirst+4*ng):(nfirst+5*ng-1)] <- as.numeric(result$back$sigma)
      #after the general model data we process the data per sample:
      n <- length(ratio)
      probs <- data.frame(result$post)
      rownames(probs) <- names(ratio) #note: does not include samples rejected because select=F
      names(probs) <- c("P0","P1","P2","P3","P4")
      #list for each sample the maximum p value and the maxgeno corresponding to it:
      probs$maxgeno<-max.col(result$post)-1  # for every row (sample) the max peak (0..4 instead of 1..5)
      probs$maxP<-maxPna # for every row the P at the max
      probslist[[model]]<-probs
      #check if there is a dip (smaller peak surrounded by larger peaks):
       gsamp<-hist(probs$maxgeno,breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5),plot=F)$counts #based on maxgeno: all samples in peaks
      top1 <- rep(NA,3); top2 <- top1
      for (i in 2:4) {
        top1[i-1] <- max(gsamp[1:(i-1)])
        top2[i-1] <- max(gsamp[(i+1):5])
      }
      mintop <- pmin(top1,top2)  #pmin is "parallel min"
      dips <- (mintop>(0.01*sum(gsamp))) & (gsamp[2:4]<0.9*mintop) # &, not && : vector-and
      dip <- max(dips) # 1 (True) if at least one valley
      if (!dip && ptype=="p.free") {
        #...or a dip in the P's?
        for (i in 2:4) {
          top1[i-1] <- max(result$psi$p[1:(i-1)])
          top2[i-1] <- max(result$psi$p[(i+1):5])
        }
        mintop <- pmin(top1,top2)
        dips <- (mintop>0.01) & (result$psi$p[2:4]<0.9*mintop)
        dip <- max(dips) # 1 (True) if at least one valley
      }
      modelresult$dip[model] <- dip
      # calculate the actual mean, sd and P of the samples in each peak
      asr <- asin(sqrt(ratio))
      act <- tapply(asr,INDEX=probs$maxgeno,FUN=mean) #mean per peak
      nam <- as.numeric(names(act))
      nfirst <- which(names(modelresult)=="muact0")
      mutrans0 <- which(names(modelresult)=="mutrans0")
      for (n in 0:4) {
        w<-which(nam==n)
        if (length(w)==1)
          modelresult[[nfirst+n]][model] <- act[w]
        else modelresult[[nfirst+n]][model] <- modelresult[[mutrans0+n]][model] #if no samples in peak use model mean
      }
      act <- tapply(asr,INDEX=probs$maxgeno,FUN=sd) #sd per peak
      nam <- as.numeric(names(act))
      nfirst <- which(names(modelresult)=="sdact0")
      for (n in 0:4) {
        w<-which(nam==n)
        if (length(w)==1)
          modelresult[[nfirst+n]][model] <- act[w]
      }
      act <- tapply(asr,INDEX=probs$maxgeno,FUN=length.noNA) #counts per peak
      totact <- sum(act) 
      if (totact>0) {
        nam <- as.numeric(names(act))
        nfirst <- which(names(modelresult)=="Pact0")
        for (n in 0:4) {
          w<-which(nam==n)
          if (length(w)==1)
            modelresult[[nfirst+n]][model] <- act[w]/totact
          else modelresult[[nfirst+n]][model] <- 0 #no samples in peak
        }
      }  
    } #ng==5
  } #CodomMarker succesful
  list(stats=modelresult,probs=probslist)
}
