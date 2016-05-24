
################################AUC.P-value.permutation
compute.permutation.table <- function(dattable){
  datr <- dattable
  for (j in 1:(ncol(dattable)-1)){
    datr[,j] <- sample(datr[,j],nrow(datr),replace=F)
  }
  return(datr)
}


compute.auc.permutation <- function(aucs,dattable,repetitions=1000){

  minp <- 1/repetitions
  aucmax <- rep(0,repetitions)
  for (i in 1:repetitions){
    datr <- compute.permutation.table(dattable)
    aucvalsrs <- compute.aucs(datr)
    aucmax[i] <- max(aucvalsrs[,2])
    #aucmax[i] <- max(aucvalsrs)
  }

  p.values <- rep(1,length(aucs))
  for (i in 1:length(aucs)){
    p.values[i] <- max(c(length(subset(aucmax,aucmax>aucs[i]))/repetitions,minp))
  }
  return(p.values)
}

################################AUC.P-value.random


multi.random.aucs.lab <- function(labs,repetitions){

  n<-length(labs)

  random.aucs.n <- function(scores){
    aucv<-multiclass.roc(labs,scores)$auc
    return(max(c(aucv,1-aucv)))
  }

  x <- matrix(runif(n*repetitions),nrow=n)
  aucs <- sapply(data.frame(x),random.aucs.n)

  ord <- order(aucs,decreasing=T)
  return(aucs[ord])
}

bh.correction <- function(p.values){
  psort <- sort(p.values,index.return=T)
  p.values.bh <- p.values
  previous.p.value <- 0
  for (i in 1:length(p.values)){
    p.values.bh[psort$ix[i]] <- p.values[psort$ix[i]]*(length(p.values)+1-i)
    if (p.values.bh[psort$ix[i]]>1){p.values.bh[psort$ix[i]]<- 1}
    if (p.values.bh[psort$ix[i]]<previous.p.value){p.values.bh[psort$ix[i]]<- previous.p.value}
    previous.p.value <- p.values.bh[psort$ix[i]]
  }
  return(p.values.bh)
}

compute.auc.random <- function(aucs,dattable,repetitions=10000,correction="none"){###repetitions correction
    ####n und prevalence<-get.n.and.prevalence:n.prev[1]
    ####aucs<-compute.aucs : auc.val[,2]
  labs=dattable[,ncol(dattable)]

  aucs.rand <-multi.random.aucs.lab(labs,repetitions)

  minpv <- 1/repetitions
  pvalues.raw <- rep(1,length(aucs))
  for (i in 1:length(aucs)){
    pvalues.raw[i] <- max(c(length(subset(aucs.rand,aucs.rand>aucs[i]))/repetitions,minpv))
  }
  if (correction=="bonferroniholm"){return(bh.correction(pvalues.raw))}
  if (correction=="bonferroni"){return(pvalues.raw/length(pvalues.raw))}
  return(pvalues.raw)#####nur show the p-values
}
