## Rescaled multistage bootstrap
## Preston http://www.statcan.gc.ca/pub/12-001-x/2009002/article/11044-eng.pdf
##

mrbweights<-function(clusters,stratas,fpcs, replicates=50, multicore=getOption("survey.multicore")){
  nstages<-NCOL(clusters)
  if (is.null(fpcs$popsize)){
    warning("Design is sampled with replacement: only first stage used")
    fpcs$popsize<-matrix(Inf, ncol=1,nrow=NROW(clusters))
    nstages<-1
  }
  
  if (multicore & !require("parallel", quietly=TRUE))
    multicore<-FALSE
  do.it<-if(multicore) mclapply else lapply
  
  weightlist<-do.it(1:replicates, function(k){
    weights<-matrix(1,nrow=NROW(clusters),ncol=nstages)
    kept<-rep(TRUE, NROW(clusters)) 
    cumffs<-rep(1,NROW(clusters))
    for(i in 1:nstages){
      ustrata<-unique(stratas[,i])
      nstrata<-length(ustrata)
      for(j in 1:nstrata){
        thisstratum<-stratas[,i]==ustrata[j]
        su <- unique(clusters[thisstratum & kept,i] )
        n <-length(su)
        nstar<-floor(n/2)
        cumff<-cumffs[thisstratum][1]
        if (nstar==0) {
          wstar<-0
          keep<- rep(FALSE,sum(thisstratum))
        } else {
          fpc<- fpcs$sampsize[thisstratum,i][1]/fpcs$popsize[thisstratum,i][1]
          lambda<-sqrt(cumff*nstar*(1-fpc)/(n-nstar))
          keep<-clusters[thisstratum,i] %in% sample(su,nstar)
          wstar<-(-lambda+lambda*(n/nstar)*keep)
        }
        weights[thisstratum, i]<-wstar*weights[thisstratum, i]
        if (nstar>0 & i<nstages){
          weights[thisstratum & kept,(i+1):nstages]<-weights[thisstratum & kept,(i+1):nstages]*sqrt(n/nstar)*keep
        }
        kept[thisstratum] <- kept[thisstratum] & keep
        cumffs[thisstratum]<-cumffs[thisstratum] * fpc   	
      }		
      
    }
    rowSums(weights)
  })	
  list(repweights=(1+do.call(cbind,weightlist)), scale=1, rscales=rep(1/(replicates-1),replicates))
}	

