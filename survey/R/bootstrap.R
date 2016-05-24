
subbootstratum<-function(psu,replicates){
 upsu<-sample(unique(psu))
 n<-length(upsu)
 replicate(replicates,
           table(factor(sample(upsu, length(upsu)-1,replace=TRUE),
                        levels=unique(psu))))*n/(n-1)
}

bootstratum<-function(psu, popsize, replicates){
  upsu<-sample(unique(psu))
  if (is.null(popsize)){
    replicate(replicates,
              table(factor(sample(upsu,length(upsu),replace=TRUE),
                                       levels=unique(psu))))
  } else {
    replicate(replicates,
              table(factor(sample(rep(upsu,length=popsize), length(upsu)),
                           levels=unique(psu))))
  }
}

bootweights<-function(strata, psu, replicates=50, fpc=NULL,
                      fpctype=c("population","fraction","correction"),
                      compress=TRUE){

  fpctype<-match.arg(fpctype)

  index<-match(psu,psu[!duplicated(psu)])
  upsu<-unique(psu)

  strata<-as.character(strata)
  weights<-matrix(nrow=length(upsu),ncol=replicates)
  ustrata<-strata[!duplicated(psu)]
  ufpc<-fpc[!duplicated(psu)]
  
  for(s in unique(ustrata)){
    this.stratum<-ustrata==s
    npsu<-length(unique(upsu[this.stratum]))

    if (is.null(fpc))
      weights[this.stratum,]<-bootstratum(upsu[this.stratum],NULL,replicates)
    else {
      this.fpc<-ufpc[this.stratum]
      if (length(unique(this.fpc))>1)
        stop("More than one fpc in stratum",s)
      this.fpc<-this.fpc[1]
      if (fpctype=="population" && this.fpc<npsu)
        stop("Population size smaller than sample size in stratum",s)
      this.fpc <-switch(fpctype,
                   population=this.fpc,
                   fraction=npsu/this.fpc,
                   correction=1-npsu/this.fpc)
      if (this.fpc> 100*npsu)
        warning("Sampling fraction <1% in stratum",s," treated as zero.")
      weights[this.stratum,]<-bootstratum(upsu[this.stratum],
                                          popsize=this.fpc,replicates=replicates)
    }
  }

  ## harmonic mean of stratum sizes
  psu.per.strata<-1/mean(1/table(ustrata))
  
  if (compress){
    rw<-list(weights=weights,index=index)
    class(rw)<-"repweights_compressed"
  } else {
    rw<-weights[index,]
  }

  list(repweights=rw, scale=psu.per.strata/((psu.per.strata-1)*(replicates-1)),
       rscales=rep(1,replicates))
}

subbootweights<-function(strata, psu, replicates=50,
                      compress=TRUE){


  index<-match(psu,psu[!duplicated(psu)])
  upsu<-unique(psu)

  strata<-as.character(strata)
  weights<-matrix(nrow=length(upsu),ncol=replicates)
  ustrata<-strata[!duplicated(psu)]
  
  for(s in unique(ustrata)){
    this.stratum<-ustrata==s
    npsu<-length(unique(upsu[this.stratum]))

    weights[this.stratum,]<-subbootstratum(upsu[this.stratum],replicates)
    
  }

  if (compress){
    rw<-list(weights=weights,index=index)
    class(rw)<-"repweights_compressed"
  } else {
    rw<-weights[index,]
  }

  list(repweights=rw, scale=1/(replicates-1),
       rscales=rep(1,replicates))
}
