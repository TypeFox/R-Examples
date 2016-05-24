
  
hadamard.doubler<-function(H){
    rbind(cbind(H,H),cbind(H,1-H))
  }
  
hadamard<- function(n){
    m<-n-(n %% 4)
    ## hadamard.list, hadamard.sizes in sysdata.rda
    precooked<- which(m < hadamard.sizes & m+4 >=hadamard.sizes)
    if (length(precooked))
      return(hadamard.list[[min(precooked)]])
    if (all(m<hadamard.sizes))
      return(hadamard.list[[1]])
    
    sizes<-hadamard.sizes*2^pmax(0,ceiling(log((n+1)/hadamard.sizes,2)) )
    bestfit<- which.min(sizes-n)
    H<-NULL
    if (sizes[bestfit]-n >4)
        H<-paley(n,sizes[bestfit])
    if (is.null(H)){
      ndoubles<-ceiling(log(sizes/hadamard.sizes, 2))[bestfit]
      H<-hadamard.list[[bestfit]]
      for(i in seq(length=ndoubles))
        H<-hadamard.doubler(H)
    }
    H
  }



jk1weights<-function(psu, fpc=NULL,
                     fpctype=c("population","fraction","correction"),
                     compress=TRUE){
  fpctype<-match.arg(fpctype)
  unq<-unique(psu)
  n<-length(unq)
  if (is.null(fpc))
    fpc<-1
  else {
    fpc<-unique(fpc)
    if (length(fpc)>1) stop("More than one fpc value given")
    if (fpc<0) stop("Negative finite population correction")
    if (fpctype=="population" && fpc<n) stop("Population size smaller than sample size. No can do.")
    fpc <-switch(fpctype, population=(fpc-n)/fpc, fraction=1-fpc, correction=fpc)
 }

  if (compress){
      if(fpc==0 && getOption("survey.drop.replicates")) ## exhaustively sampled strata do not need replicates.
         repweights<-matrix(ncol=0,nrow=length(psu))
      else {
          repweights<-matrix(n/(n-1),n,n)
          diag(repweights)<-0
      }
      repweights<-list(weights=repweights, index=match(psu,unq))
      class(repweights)<-c("repweights_compressed","repweights")
      rval<-list(type="jk1", repweights=repweights,scale=fpc*(n-1)/n)
      rval
  } else {
      if(fpc==0 && getOption("survey.drop.replicates")) ## exhaustively sampled strata do not need replicates.
          return(list(type="jk1",repweights=matrix(ncol=0,nrow=length(psu)), scale=0))
      repweights<-outer(psu, unq, "!=")*n/(n-1)
      class(repweights)<-"repweights"
      rval<-list(type="jk1", repweights=repweights,scale=(fpc*(n-1)/n))
      rval
  }
}




jknweights<-function(strata,psu, fpc=NULL,
                     fpctype=c("population","fraction","correction"),
                     compress=TRUE, lonely.psu=getOption("survey.lonely.psu")){

  sunq<-unique(strata)
  unq<-unique(psu)
  nstrat<-length(sunq)
  n<-length(strata)

  lonely.psu<-match.arg(lonely.psu, c("fail","certainty","remove","adjust","average"))
  
  fpctype<-match.arg(fpctype)
  
  if (is.null(fpc)){
      fpc<-rep(1,nstrat)
      names(fpc)<-as.character(sunq)
      fpctype<-"correction"
  } else if (length(fpc)==n){
      if (length(unique(fpc))>nstrat)
          stop("More distinct fpc values than strata")
      fpc<-sapply(sunq, function(ss) fpc[match(ss,strata)])
      names(fpc)<-as.character(sunq)
  } else if (length(fpc)==1) {
      fpc<-rep(fpc,nstrat)
      names(fpc)<-as.character(sunq)
  } else if (length(fpc)==nstrat){
      nn<-names(fpc)
      if (is.null(nn)) names(fpc)<-as.character(sunq)
      if (!all(names(fpc) %in% as.character(sunq)))
          stop("fpc has names that do not match the stratum identifiers")
  }


  if (compress){
    repweights<-matrix(1,ncol=length(unq),nrow=length(unq))
  } else {
    repweights<-matrix(1,ncol=length(unq), nrow=length(psu))
  }
  counter<-0
  rscales<-numeric(length(unq))
  
  for(ss in as.character(sunq)){
      thisfpc<-fpc[match(ss,names(fpc))]
      theseweights<-jk1weights(psu[strata %in% ss], fpc=thisfpc,
                               fpctype=fpctype,compress=compress)
      nc<-if (compress) NCOL(theseweights$repweights$weights) else NCOL(theseweights$repweights)
      if (nc==1 && thisfpc!=0){
        ## lonely PSUs
        if (lonely.psu=="fail")
          stop("Stratum",ss,"has only one PSU")
        if (lonely.psu=="remove")
          next
        if (lonely.psu=="certainty")
          next
        if (lonely.psu=="average")
          next
        if (lonely.psu=="adjust"){
          nc<-1
          if (compress)
            repweights[, counter+nc]<-ifelse(strata[!duplicated(psu)] %in% ss, 0, nstrat/(nstrat-1))
          else
            repweights[ counter+nc]<-ifelse(strata %in% ss, 0, nstrat/(nstrat-1))
          rscales[counter+nc]<-(nstrat-1)/nstrat
          counter<-counter+nc
          next
        }
      }
      
      if (compress)
        repweights[strata[!duplicated(psu)] %in% ss,counter+seq(length=nc)]<-theseweights$repweights$weights
      else
        repweights[strata %in% ss, counter+seq(length=nc)]<-theseweights$repweights
      
      rscales[counter+seq(length=nc)]<-theseweights$scale
      counter<-counter+nc
  }
  if (counter==0) stop("All strata were exhaustively sampled: you have the whole population")
  scale<-1
  if (counter<length(unq)){
      repweights<-repweights[,1:counter]
      rscales<-rscales[1:counter]
      if (lonely.psu=="average")
        scale<-scale*length(unq)/counter
  }
  if (compress){
    repweights<-list(weights=repweights,index=match(psu,unq))
    class(repweights)<- c("repweights_compressed","repweights")
  } else class(repweights)<-"repweights"
  list(type="jkn", repweights=repweights, rscales=rscales, scale=scale)
}



brrweights<-function(strata,psu, match=NULL, small=c("fail","split","merge"),
                     large=c("split","merge","fail"), fay.rho=0,
                     only.weights=FALSE, compress=TRUE,
                     hadamard.matrix=NULL){

  small<-match.arg(small)
  large<-match.arg(large)

  strata<-as.character(strata)
  
  ssize<-table(strata[!duplicated(psu)])
  if (any(ssize<2) && small=="fail")
    stop("Some strata have fewer than 2 PSUs")
  if (any(ssize>2) && large=="fail")
    stop("Some strata have more than 2 PSUs")

  unq<-which(!duplicated(psu))
  sunq<-strata[unq]
  psunq<-psu[unq]
  weights<-matrix(ncol=2,nrow=length(unq))
  weightstrata<-numeric(length(unq))
  
  if (length(match)==length(strata))
    match<-match[unq]
  if (is.null(match))
    match<-unq  ## default is to match by dataset order
  oo<-order(sunq,match)

  upto <- 0
  
  if(any(ssize==1)){
    smallstrata<-names(ssize)[ssize==1]
    if(small=="split"){
      weights[sunq %in% smallstrata,1]<- 0.5
      weights[sunq %in% smallstrata,2]<- 0.5
      weightstrata[sunq %in% smallstrata]<-1:length(smallstrata)
      upto<-length(smallstrata)
    } else {
      ##small=="merge"
      if (length(smallstrata) > 1){
        weights[oo,][sunq[oo] %in% smallstrata, 1]<-rep(0:1,length.out=length(smallstrata))
        weights[oo,][sunq[oo] %in% smallstrata, 2]<-rep(1:0,length.out=length(smallstrata))
        if(length(smallstrata) %% 2==0)
          weightstrata[oo][sunq[oo] %in% smallstrata]<-rep(1:(length(smallstrata) %/%2), 2)
        else
          weightstrata[oo][sunq[oo] %in% smallstrata]<-c(1,rep(1:(length(smallstrata) %/%2), 2))
        upto<-length(smallstrata) %/% 2
      } else stop("Can't merge with a single small stratum")
    }
  }

  if (any(ssize>2)){
    largestrata<-names(ssize)[ssize>2]
    if (large=="split"){
      if (any(ssize[largestrata] %%2 ==1))
        stop("Can't split with odd numbers of PSUs in a stratum")
      ## make substrata of size 2
      for(ss in largestrata){
        weights[oo,][sunq[oo] %in% ss, 1]<-rep(0:1,length.out=ssize[ss])
        weights[oo,][sunq[oo] %in% ss, 2]<-rep(1:0,length.out=ssize[ss])
        weightstrata[oo][sunq[oo] %in% ss]<-upto+rep(1:(ssize[ss] %/%2),each=2)
        upto<-upto+(ssize[ss] %/% 2)
      }
    } else {
      ## make two substrata.
      halfsize<-ssize[largestrata] %/%2
      otherhalfsize<-ssize[largestrata] - halfsize
      reps<-as.vector(rbind(halfsize,otherhalfsize))
      nlarge<-length(halfsize)
      weights[oo,][sunq[oo] %in% largestrata, 1]<-rep(rep(0:1,nlarge),reps) 
      weights[oo,][sunq[oo] %in% largestrata, 2]<-rep(rep(1:0,nlarge),reps)
      weightstrata[oo][sunq[oo] %in% largestrata]<-upto+rep(1:length(largestrata),ssize[largestrata])
      upto<-upto+length(largestrata)
    }
  }
  if(any(ssize==2)){
    goodstrata<-names(ssize)[ssize==2]
    weights[oo,][sunq[oo] %in% goodstrata, 1]<-rep(0:1,length(goodstrata))
    weights[oo,][sunq[oo] %in% goodstrata, 2]<-rep(1:0,length(goodstrata))
    weightstrata[oo][sunq[oo] %in% goodstrata]<-upto+rep(1:length(goodstrata),each=2)
    upto<-upto+length(goodstrata)
  }
  

  if (is.null(hadamard.matrix)){
    H<-hadamard(upto)
  } else {
    ## the user supplied hadamard.matrix
    ## Check that it is a binary matrix and satifies the
    ## Hadamard determinant property
    if (!is.matrix(hadamard.matrix) || nrow(hadamard.matrix)<upto+1)
      stop("hadamard.matrix must be a matrix of dimension at least nstrata+1")
    values<-unique(as.vector(hadamard.matrix))
    if(length(values)!=2)
      stop("hadamard.matrix has more than two different values")
    H<-ifelse(hadamard.matrix==values[1],0,1)
    if(!is.hadamard(H,full.orthogonal.balance=FALSE))
        stop("hadamard.matrix is not a Hadamard matrix")
  }
  ii<-1:upto
  jj<-1:length(weightstrata)
  sampler<-function(i){
    h<-H[1+ii, i]+1
    col<-h[match(weightstrata,ii)]
    wa<-weights[cbind(jj,col)]
    wb<-weights[cbind(jj,3-col)]
    if (compress)
      wa*(2-fay.rho)+wb*fay.rho
    else
      wa[match(psu,psunq)]*(2-fay.rho)+wb[match(psu,psunq)]*fay.rho
  }

  if (only.weights){
    repweights<-sapply(1:NCOL(H),sampler)
    if (!compress){
      class(repweights)<-"repweights"
      return(repweights)
    }
    repweights<-list(weights=repweights,index=match(psu,psunq))
    class(repweights)<-c("repweights_compressed","repweights")
    repweights
  }  else 
  list(weights=weights, wstrata=weightstrata, strata=sunq, psu=psunq,
       npairs=NCOL(H),sampler=sampler,compress=compress)

}
  

  


##
## Designs with replication weights rather than survey structure.
##

as.svrepdesign<- function(design,type=c("auto","JK1","JKn","BRR","bootstrap","subbootstrap","mrbbootstrap","Fay"),
                          fay.rho=0, fpc=NULL, fpctype=NULL,...,compress=TRUE, mse=getOption("survey.replicates.mse")){

  type<-match.arg(type)

  if (type=="auto"){
    if (!design$has.strata)
      type<-"JK1"
    else
      type<-"JKn"
  }
  selfrep<-NULL
  if (type=="JK1" && design$has.strata)
    stop("Can't use JK1 for a stratified design")
  if (type %in% c("JKn","BRR","Fay") && !design$has.strata)
    stop("Must use JK1 or bootstrap for an unstratified design")
  
  if (is.null(fpc)) {
    fpctype<-"population"
    
    if (is.null(design$fpc) ||
        (inherits(design, "survey.design2") && is.null(design$fpc$popsize))){
      fpc<-NULL
    } else if (type %in% c("Fay","BRR")){
      warning("Finite population correction dropped in conversion")
    } else {
      if (inherits(design,"survey.design2")){
        fpc<-design$fpc$popsize
        if(NCOL(fpc)>1 && type!="mrbbootstrap"){
          fpc<-fpc[,1]
          warning("Finite population corrections after first stage have been dropped")
        }
        if (getOption("survey.drop.replicates")){
          selfrep<-design$fpc$popsize[,1]==design$fpc$sampsize[,1]
        } 
      } else{
        fpc<-design$fpc[,2]
        names(fpc)<-design$fpc[,1]
      }
    }
  } else {
    if (type %in% c("Fay","BRR","subbootstrap")) stop(paste("fpc information cannot be used for type=",type))
    if (is.null(fpctype)) stop("fpctype must be supplied if fpc is supplied")
  }
  if (type=="JK1"){
    ##JK1
    r<-jk1weights(design$cluster[,1], fpc=fpc,fpctype=fpctype, compress=compress)
    repweights<-r$repweights
    scale<-drop(r$scale)
    if (inherits(repweights,"repweights_compressed"))
      rscales<-rep(1, NCOL(repweights$weights))
    else
      rscales<-rep(1, NCOL(repweights))
    
    type<-"JK1"
    pweights<-1/design$prob
  } else if (type %in% c("BRR","Fay")){
    ##BRR
    if (inherits(design,"survey.design2"))
      repweights<-brrweights(design$strata[,1], design$cluster[,1],...,fay.rho=fay.rho,
                             compress=compress,only.weights=TRUE)
    else
      repweights<-brrweights(design$strata, design$cluster[,1],...,fay.rho=fay.rho,
                             compress=compress,only.weights=TRUE)
    pweights<-1/design$prob
    if (length(pweights)==1)
      pweights<-rep(pweights, NROW(design$variables))
    
    if (fay.rho==0)
      type<-"BRR"
    else
      type<-"Fay"

    rscales<-rep(1,ncol(repweights))
    scale<-1/(ncol(repweights)*(1-fay.rho)^2)
    
  } else if (type=="JKn"){
    ##JKn
    if (inherits(design,"survey.design2"))
      r<-jknweights(design$strata[,1],design$cluster[,1], fpc=fpc,
                    fpctype=fpctype, compress=compress)
    else
      r<-jknweights(design$strata,design$cluster[,1], fpc=fpc,
                    fpctype=fpctype, compress=compress)
    pweights<-1/design$prob
    repweights<-r$repweights
    scale<-r$scale
    rscales<-r$rscales
  } else if (type=="bootstrap"){
    ##bootstrap
    if (inherits(design,"survey.design2"))
      r<-bootweights(design$strata[,1],design$cluster[,1], fpc=fpc,
                     fpctype=fpctype, compress=compress,...)
    else
      r<-bootweights(design$strata,design$cluster[,1], fpc=fpc,
                     fpctype=fpctype, compress=compress,...)
    pweights<-1/design$prob
    repweights<-r$repweights
    scale<-r$scale
    rscales<-r$rscales
  }else if (type=="subbootstrap"){
    ##bootstrap
    if (inherits(design,"survey.design2"))
      r<-subbootweights(design$strata[,1],design$cluster[,1],compress=compress,...)
    else
      r<-subbootweights(design$strata,design$cluster[,1],compress=compress,...)
    pweights<-1/design$prob
    repweights<-r$repweights
    scale<-r$scale
    rscales<-r$rscales
  } else if (type=="mrbbootstrap"){
    if (inherits(design,"survey.design2"))
      r<-mrbweights(design$cluster,design$strata,design$fpc,...)
    else
      stop("MRB bootstrap not available for obsolete svydesign objects")
    pweights<-1/design$prob
    repweights<-r$repweights
    scale<-r$scale
    rscales<-r$rscales
  } else stop("Can't happen")

  rval<-list(repweights=repweights, pweights=pweights,
             type=type, rho=fay.rho,scale=scale, rscales=rscales,
             call=sys.call(), combined.weights=FALSE, selfrep=selfrep,mse=mse)
  rval$variables <- design$variables
  class(rval)<-"svyrep.design"
  rval$degf<-degf(rval)
  rval
}



svrepdesign<-function(variables, repweights, weights,data=NULL,...) UseMethod("svrepdesign",data)

svrepdesign.default<-function(variables=NULL,repweights=NULL, weights=NULL,
                              data=NULL,type=c("BRR","Fay","JK1", "JKn","bootstrap","other"),
                              combined.weights=TRUE, rho=NULL, bootstrap.average=NULL,
                              scale=NULL,rscales=NULL,fpc=NULL, fpctype=c("fraction","correction"),
                              mse=getOption("survey.replicates.mse"),...)
{
  
  type<-match.arg(type)
  
  if(type=="Fay" && is.null(rho))
    stop("With type='Fay' you must supply the correct rho")
  
  if (type %in% c("JK1","JKn")  && !is.null(rho))
    warning("rho not relevant to JK1 design: ignored.")
  
  if (type %in% c("other")  && !is.null(rho))
    warning("rho ignored.")

  
  if(is.null(variables))
    variables<-data
    
  if(inherits(variables,"formula")){
    mf<-substitute(model.frame(variables, data=data,na.action=na.pass))
    variables<-eval.parent(mf)
  }
    
  if(inherits(repweights,"formula")){
    mf<-substitute(model.frame(repweights, data=data))
    repweights<-eval.parent(mf)
    repweights<-na.fail(repweights)
  }

  if(is.character(repweights)){##regular expression
    wtcols<-grep(repweights,names(data))
    repweights<-data[,wtcols]
  }
  
  if (is.null(repweights))
    stop("You must provide replication weights")
  
  
  if(inherits(weights,"formula")){
    mf<-substitute(model.frame(weights, data=data))
    weights<-eval.parent(mf)
    weights<-drop(as.matrix(na.fail(weights)))
  }

  if (is.null(weights)){
    warning("No sampling weights provided: equal probability assumed")
    weights<-rep(1,NROW(repweights))
  }

  repwtmn<-mean(apply(repweights,2,mean))
  wtmn<-mean(weights)
  probably.combined.weights<-(repwtmn>5) & (wtmn/repwtmn<5)
  probably.not.combined.weights<-(repwtmn<5) & (wtmn/repwtmn>5)
  if (combined.weights & probably.not.combined.weights)
    warning(paste("Data do not look like combined weights: mean replication weight is", repwtmn," and mean sampling weight is",wtmn))
  if (!combined.weights & probably.combined.weights)
    warning(paste("Data look like combined weights: mean replication weight is", repwtmn," and mean sampling weight is",wtmn))

  if (!is.null(rscales) && !(length(rscales) %in% c(1, ncol(repweights)))){
    stop(paste("rscales has length ",length(rscales),", should be ncol(repweights)",sep=""))
  }
  
  if (type == "BRR")
    scale<-1/ncol(repweights)
  if (type=="Fay")
    scale <-1/(ncol(repweights)*(1-rho)^2)

  if (type=="bootstrap"){
    if(is.null(bootstrap.average))
      bootstrap.average<-1
    if (is.null(scale))
      scale<-bootstrap.average/(ncol(repweights)-1)
    if (is.null(rscales))
      rscales<-rep(1,ncol(repweights))
    }

  if (type=="JK1" && is.null(scale)) {
    if(!combined.weights){
      warning("scale (n-1)/n not provided: guessing from weights")
      scale<-1/max(repweights[,1])
    } else {
      probably.n = ncol(repweights)
      scale<- (probably.n-1)/probably.n
      warning("scale (n-1)/n not provided: guessing n=number of replicates")
    }
  }

  if (type =="JKn" && is.null(rscales))
    if (!combined.weights) {
      warning("rscales (n-1)/n not provided:guessing from weights")
      rscales<-1/apply(repweights,2,max)
    } else stop("Must provide rscales for combined JKn weights")

  if (type=="other" && (is.null(rscales) || is.null(scale))){
    if (is.null(rscales)) rscales<-rep(1,NCOL(repweights))
    if (is.null(scale)) scale<-1
    warning("scale or rscales not specified, set to 1")
  }
  if (is.null(rscales)) rscales<-rep(1,NCOL(repweights))

  if (!is.null(fpc)){
      if (missing(fpctype)) stop("Must specify fpctype")
      fpctype<-match.arg(fpctype)
      if (type %in% c("BRR","Fay")) stop("fpc not available for this type")
      if (type %in% "bootstrap") stop("Separate fpc not needed for bootstrap")
      if (length(fpc)!=length(rscales)) stop("fpc is wrong length")
      if (any(fpc>1) || any(fpc<0)) stop("Illegal fpc value")
      fpc<-switch(fpctype,correction=fpc,fraction=1-fpc)
      rscales<-rscales*fpc
  }
  
  
  rval<-list(type=type, scale=scale, rscales=rscales,  rho=rho,call=sys.call(),
             combined.weights=combined.weights)
  rval$variables<-variables
  rval$pweights<-weights
  if (!inherits(repweights,"repweights"))
    class(rval)<-"repweights"
  rval$repweights<-repweights
  class(rval)<-"svyrep.design"
  rval$degf<-degf(rval)
  rval$mse<-mse
  rval
  
}


print.svyrep.design<-function(x,...){
  cat("Call: ")
  print(x$call)
  if (x$type=="Fay")
    cat("Fay's variance method (rho=",x$rho,") ")
  if (x$type=="BRR")
    cat("Balanced Repeated Replicates ")
  if (x$type=="JK1")
    cat("Unstratified cluster jacknife (JK1) ")
  if (x$type=="JKn")
    cat("Stratified cluster jackknife (JKn) ")
  if (x$type=="bootstrap")
    cat("Survey bootstrap ")
  if (x$type=="mrbbootstrap")
    cat("Multistage rescaled bootstrap ")
  if (x$type=="subbootstrap")
    cat("(n-1) bootstrap ")
  nweights<-ncol(x$repweights)
  cat("with", nweights,"replicates")
  if (!is.null(x$mse) && x$mse) cat(" and MSE variances")
  cat(".\n")
  invisible(x)
}

summary.svyrep.design<-function(object,...){
  class(object)<-c("summary.svyrep.design", class(object))
  object
}

print.summary.svyrep.design<-function(x,...){
  class(x)<-class(x)[-1]
  print(x)
  cat("Variables: \n")
  print(colnames(x))
}



image.svyrep.design<-function(x, ..., col=grey(seq(.5,1,length=30)),
                              type.=c("rep","total")){
  type<-match.arg(type.)
  m<-as.matrix(x$repweights)
  if (type=="total"){
    m<-m*x$pweights
  } 
  
  image(1:NCOL(m), 1:NROW(m), t(m),  col=col, xlab="Replicate", ylab="Observation",...)
  invisible(NULL)
}

"[.svyrep.design"<-function(x, i, j, drop=FALSE){
  if (!missing(i)){
    pwt<-x$pweights
    if (is.data.frame(pwt)) pwt<-pwt[[1]]
    x$pweights<-pwt[i]
    x$repweights<-x$repweights[i,,drop=FALSE]
    if(!is.null(x$selfrep))
        x$selfrep<-x$selfrep[i]
    if (!missing(j))
      x$variables<-x$variables[i,j, drop=FALSE]
    else
      x$variables<-x$variables[i,,drop=FALSE]
    x$degf<-NULL
    x$degf<-degf(x)
  } else {
    x$variables<-x$variables[,j,drop=FALSE]
  }
  x
}


subset.svyrep.design<-function(x,subset,...){
        e <- substitute(subset)
        r <- eval(e, x$variables, parent.frame())
        r <- r & !is.na(r) 
        x<-x[r,]
	x$call<-sys.call(-1)
	x
}

update.svyrep.design<-function(object,...){

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    object$variables[,newnames[j]]<-eval(dots[[j]],object$variables, parent.frame())
  }
  
  object$call<-sys.call(-1)
  object 
}

weights.svyrep.design<-function(object,type=c("replication","sampling","analysis"),...){
  type<-match.arg(type)
  switch(type,
         replication= as.matrix(object$repweights),
         sampling=object$pweights,
         analysis=if(object$combined.weights) as.matrix(object$repweights) else as.matrix(object$repweights)*object$pweights)
}

weights.survey.design<-function(object,...){
  return(1/object$prob)
}


svyquantile.svyrep.design<-svrepquantile<-function(x,design,quantiles,method="linear",
                                                   interval.type=c("probability","quantile"),f=1,
                                                   return.replicates=FALSE,
                                                   ties=c("discrete","rounded"),na.rm=FALSE,...){

  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svyquantile")

  ties<-match.arg(ties)
  interval<-match.arg(interval.type)
  if (design$type %in% c("JK1","JKn") && interval=="quantile")
    warning("Jackknife replicate weights may not give valid standard errors for quantiles")
  if (design$type %in% "other" && interval=="quantile")
    warning("Not all replicate weight designs give valid standard errors for quantiles.")
  if (inherits(x,"formula"))
		x<-model.frame(x,design$variables,na.action=if(na.rm) na.pass else na.fail)
    else if(typeof(x) %in% c("expression","symbol"))
        x<-eval(x, design$variables)


   if (na.rm){
         nas<-rowSums(is.na(x))
       design<-design[nas==0,]
        if (length(nas)>length(design$prob))
          x<-x[nas==0,,drop=FALSE]
        else
          x[nas>0,]<-0
      }
    
  if (NROW(x)<=1){
      rval<-matrix(rep(as.matrix(x),length(quantiles)),ncol=NCOL(x),nrow=length(quantiles),byrow=TRUE)
      dimnames(rval)<-list(paste("q",round(quantiles,2),sep=""), names(x))
      if (getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep))
          vv<-matrix(0,ncol=NCOL(x),nrow=length(quantiles))
      else
          vv<-matrix(NA,ncol=NCOL(x),nrow=length(quantiles))
      dimnames(vv)<-list(paste("q",round(quantiles,2),sep=""), names(x))
      attr(rval,"var")<-vv
      attr(rval,"statistic")<-quantiles
      if (return.replicates)
          rval<-list(mean=rval,replicates=NULL)
      class(rval)<-"svrepstat"
      return(rval)
  }

  
    w<-weights(design,"analysis")

    if (interval=="quantile"){
      ## interval on quantile scale
      if (ties=="discrete")
        computeQuantiles<-function(xx){
          oo<-order(xx)
          
          ws<-weights(design,"sampling")
          cum.ws<-cumsum(ws[oo])/sum(ws)
          rval<-approx(cum.ws,xx[oo],method=method,f=f,
                       yleft=min(xx),yright=max(xx),
                       xout=quantiles,ties=min)$y
          
          cum.w<-apply(w,2,function(wi) cumsum(wi[oo])/sum(wi))
          
          qq<-apply(cum.w, 2,function(cum.wi) approx(cum.wi,xx[oo],method=method,f=f,
                                                     yleft=min(xx),yright=max(xx),
                                                     xout=quantiles,ties=min)$y)
          if (length(quantiles)>1)
            qq<-t(qq)
          else
            qq<-as.matrix(qq)
          ##rval<-colMeans(qq)
          
          rval<-list(quantiles=rval,
                     variances=diag(as.matrix(svrVar(qq,design$scale,design$rscales,mse=design$mse,coef=rval))))
          if (return.replicates)
            rval<-c(rval, list(replicates=qq))
          rval
        } else { ##ties="rounded"
          computeQuantiles<-function(xx){
            ws<-weights(design,"sampling")

            wws<-rowsum(ws,xx,reorder=TRUE)
            uxx<-sort(unique(xx))
            
            cum.wws<-cumsum(wws)/sum(wws)
            rval<-approx(cum.wws,uxx,method=method,f=f,
                     yleft=min(xx),yright=max(xx),
                         xout=quantiles,ties=min)$y
            
            cum.w<-apply(rowsum(w,xx,reorder=TRUE),2,function(wi) cumsum(wi)/sum(wi))
            
            qq<-apply(cum.w, 2,function(cum.wi) approx(cum.wi,uxx,method=method,f=f,
                                                       yleft=min(xx),yright=max(xx),
                                                       xout=quantiles,ties=min)$y)
            if (length(quantiles)>1)
              qq<-t(qq)
            else
              qq<-as.matrix(qq)
            ##rval<-colMeans(qq)
            
            rval<-list(quantiles=rval,
                       variances=diag(as.matrix(svrVar(qq,design$scale,design$rscales,mse=design$mse,coef=rval))))
            if (return.replicates)
              rval<-c(rval, list(replicates=qq))
            rval
          }
        }
    } else {
      ## interval on probability scale, backtransformed.
      if (ties=="discrete"){
        computeQuantiles<-function(xx){
          oo<-order(xx)
          w<-weights(design,"sampling")
          cum.w<- cumsum(w[oo])/sum(w)
          Qf<-approxfun(cum.w,xx[oo],method=method,f=f,
                        yleft=min(xx),yright=max(xx),
                        ties=min)
          
          point.estimates<-Qf(quantiles)
          if(length(quantiles)==1)
            estfun<-as.numeric(xx<point.estimates)
          else
            estfun<-0+outer(xx,point.estimates,"<")
          est<-svymean(estfun,design, return.replicates=return.replicates)
          if (return.replicates)
            q.estimates<-matrix(Qf(est$replicates),nrow=NROW(est$replicates))
          ci<-matrix(Qf(c(coef(est)+2*SE(est), coef(est)-2*SE(est))),ncol=2)
          variances<-((ci[,1]-ci[,2])/4)^2
          rval<-list(quantiles=point.estimates,
                     variances=variances)
          if (return.replicates)
            rval<-c(rval, list(replicates=q.estimates))
          rval
        }
      } else {
        ## ties=rounded
        computeQuantiles<-function(xx){
          w<-weights(design,"sampling")
          ww<-rowsum(w,xx,reorder=TRUE)
          uxx<-sort(unique(xx))
          cum.w<- cumsum(ww)/sum(ww)
          Qf<-approxfun(cum.w,uxx,method=method,f=f,
                        yleft=min(xx),yright=max(xx),
                        ties=min)
          
          point.estimates<-Qf(quantiles)
          if(length(quantiles)==1)
            estfun<-as.numeric(xx<point.estimates)
          else
                estfun<-0+outer(xx,point.estimates,"<")
          est<-svymean(estfun, design, return.replicates=return.replicates)
          if (return.replicates)
            q.estimates<-matrix(Qf(est$replicates),nrow=NROW(est$replicates))
          ci<-matrix(Qf(c(coef(est)+2*SE(est), coef(est)-2*SE(est))),ncol=2)
          variances<-((ci[,1]-ci[,2])/4)^2
          rval<-list(quantiles=point.estimates,
                     variances=variances)
          if (return.replicates)
            rval<-c(rval, list(replicates=q.estimates))
          rval
        }
        
      }
    }

  if (!is.null(dim(x)))
    results<-apply(x,2,computeQuantiles)
  else
    results<-computeQuantiles(x)
  
  rval<-matrix(sapply(results,"[[","quantiles"),ncol=NCOL(x),nrow=length(quantiles),
               dimnames=list(paste("q",round(quantiles,2),sep=""), names(x)))
  vv<-matrix(sapply(results,"[[","variances"),ncol=NCOL(x),nrow=length(quantiles),
             dimnames=list(paste("q",round(quantiles,2),sep=""), names(x)))
  attr(rval,"var")<-vv
  attr(rval, "statistic")<-"quantiles"
  if (return.replicates) {
    reps<-do.call(cbind,lapply(results,"[[","replicates"))
    attr(reps,"scale")<-design$scale
    attr(reps,"rscales")<-design$rscales
    attr(reps,"mse")<-design$mse
    rval<-list(mean=rval, replicates=reps)
  }
  class(rval)<-"svrepstat"
  rval
  
}


svrVar<-function(thetas, scale, rscales,na.action=getOption("na.action"),mse=getOption("survey.replicates.mse"),coef){
  thetas<-get(na.action)(thetas)
  naa<-attr(thetas,"na.action")
  if (!is.null(naa)){
    rscales<-rscales[-naa]
    if (length(rscales))
      warning(length(naa), " replicates gave NA results and were discarded.")
    else
      stop("All replicates contained NAs")
  }
  if (is.null(mse)) mse<-FALSE
 
  if (length(dim(thetas))==2){
    if (mse) {
      meantheta<-coef
    } else {
      meantheta<-colMeans(thetas[rscales>0,,drop=FALSE])
    }
    v<-crossprod( sweep(thetas,2, meantheta,"-")*sqrt(rscales))*scale
  }  else {
    if (mse){
      meantheta<-coef
    } else {
      meantheta<-mean(thetas[rscales>0])
    }
    v<- sum( (thetas-meantheta)^2*rscales)*scale
  }
  attr(v,"na.replicates")<-naa
  attr(v,"means")<-meantheta
  return(v)
}


svyvar.svyrep.design<-svrepvar<-function(x, design, na.rm=FALSE, rho=NULL,
                                         return.replicates=FALSE,...,estimate.only=FALSE){

  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svyvar")
  
  if (inherits(x,"formula"))
    x<-model.frame(x,design$variables,na.action=na.pass)
  else if(typeof(x) %in% c("expression","symbol"))
    x<-eval(x, design$variables)

  wts<-design$repweights
  scale<-design$scale
  rscales<-design$rscales
  if (!is.null(rho)) .NotYetUsed("rho")
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    if(any(nas>0)){
      design<-design[nas==0,]
      x<-x[nas==0,,drop=FALSE]
      wts<-wts[nas==0,,drop=FALSE]
    }
  }
  
  if (design$combined.weights)
    pw<-1
  else
    pw<-design$pweights

  n<-NROW(x)
  p<-NCOL(x)
  v<-function(w){
    xbar<-colSums(as.vector(w)*pw*x)/sum(as.vector(w)*pw)
    xdev<-sweep(x,2,xbar,"-")
    x1<-matrix(rep(xdev,p),ncol=p*p)
    x2<-xdev[,rep(1:p,each=p),drop=FALSE]
    (n/(n-1))*colSums(x1*x2*as.vector(w)*pw)/sum(as.vector(w)*pw)
  }

  if (design$combined.weights)
    rval<-v(design$pweights)
  else
    rval<-v(rep(1,length(design$pweights)))

  rval<-matrix(rval, ncol=p)
  dimnames(rval)<-list(colnames(x),colnames(x))
  if (estimate.only) return(rval)
  
  repvars<-apply(wts,2, v)
  
  repvars<-drop(t(repvars))
  attr(rval,"var")<-svrVar(repvars, scale, rscales,mse=design$mse, coef=rval)
  attr(rval, "statistic")<-"variance"
  if (return.replicates){
    attr(repvars,"scale")<-design$scale
    attr(repvars,"rscales")<-design$rscales
    attr(repvars,"mse")<-design$mse
    rval<-list(variance=rval, replicates=repvars)
  }
  class(rval)<-c("svrepvar","svrepstat")
  rval

}


print.svrepvar<-function (x,  covariance=FALSE, ...) 
{
    if (is.list(x)) x<-x[[1]]
    vv <- as.matrix(attr(x, "var"))
    if (covariance){
      nms<-outer(rownames(x),colnames(x),paste,sep=":")
      m<-cbind(as.vector(x), sqrt(diag(vv)))
      rownames(m)<-nms
    } else{
      ii <- which(diag(sqrt(length(x)))>0)
      m <- cbind(x[ii], sqrt(diag(vv))[ii])
      if(length(ii)==1) rownames(m)<-rownames(x)
    }
    colnames(m) <- c(attr(x, "statistic"), "SE")
    printCoefmat(m)
}

as.matrix.svrepvar<-function(x,...) if (is.list(x)) unclass(x[[1]]) else unclass(x)


svymean.svyrep.design<-svrepmean<-function(x,design, na.rm=FALSE, rho=NULL,
                                           return.replicates=FALSE,deff=FALSE,...)
{
  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svymean")
  if (!inherits(design,"svyrep.design")) stop("design is not a replicate survey design")
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,design$variables,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  } else {
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else {
          if(is.data.frame(x) && any(sapply(x,is.factor))){
              xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
              cols<-sapply(xx,NCOL)
              scols<-c(0,cumsum(cols))
              cn<-character(sum(cols))
              for(i in 1:length(xx))
                  cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
              x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
              for(i in 1:length(xx)){
                  x[,scols[i]+1:cols[i]]<-xx[[i]]
              }
              colnames(x)<-cn
          }
      }
  } 
  
  x<-as.matrix(x)
  
  if (na.rm){
    nas<-rowSums(is.na(x))
    design<-design[nas==0,]
    x<-x[nas==0,,drop=FALSE]
  }
  
  wts<-design$repweights
  scale<-design$scale
  rscales<-design$rscales
  if (!is.null(rho)) .NotYetUsed("rho")
  
  if (!design$combined.weights)
    pw<-design$pweights
  else
    pw<-1
  
  rval<-colSums(design$pweights*x)/sum(design$pweights)

  if (getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep)){
    v<-matrix(0,length(rval),length(rval))
    repmeans<-NULL
  } else {
  if (inherits(wts, "repweights_compressed")){
    repmeans<-matrix(ncol=NCOL(x), nrow=ncol(wts$weights))
    for(i in 1:ncol(wts$weights)){
      wi<-wts$weights[wts$index,i]
      repmeans[i,]<-t(colSums(wi*x*pw)/sum(pw*wi))
    }
  } else {
    repmeans<-matrix(ncol=NCOL(x), nrow=ncol(wts))
    for(i in 1:ncol(wts)){
      repmeans[i,]<-t(colSums(wts[,i]*x*pw)/sum(pw*wts[,i]))
    }
  }
  repmeans<-drop(repmeans)
  v <- svrVar(repmeans, scale, rscales,mse=design$mse, coef=rval)
}
  attr(rval,"var") <-v
  attr(rval, "statistic")<-"mean"
  if (return.replicates){
    attr(repmeans,"scale")<-design$scale
    attr(repmeans,"rscales")<-design$rscales
    attr(repmeans,"mse")<-design$mse
    rval<-list(mean=rval, replicates=repmeans)
  }
  if (is.character(deff) || deff){
      nobs<-length(design$pweights)
      npop<-sum(design$pweights)
      vsrs<-unclass(svyvar(x,design,na.rm=na.rm, return.replicates=FALSE,estimate.only=TRUE))/length(design$pweights)
      if (deff!="replace")
        vsrs<-vsrs*(npop-nobs)/npop
      attr(rval,"deff") <- v/vsrs
  }
  class(rval)<-"svrepstat"
  rval
}



svytotal.svyrep.design<-svreptotal<-function(x,design, na.rm=FALSE, rho=NULL,
                                             return.replicates=FALSE, deff=FALSE,...)
{
 if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svytotal")
 if (!inherits(design,"svyrep.design")) stop("design is not a replicate survey design")
 
 if (inherits(x,"formula")){
     ## do the right thing with factors
     mf<-model.frame(x,design$variables,na.action=na.pass)
     xx<-lapply(attr(terms(x),"variables")[-1],
                function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
     x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
     scols<-c(0,cumsum(cols))
     for(i in 1:length(xx)){
         x[,scols[i]+1:cols[i]]<-xx[[i]]
     }
     colnames(x)<-do.call("c",lapply(xx,colnames))
 }  else{
      if(typeof(x) %in% c("expression","symbol"))
          x<-eval(x, design$variables)
      else {
          if(is.data.frame(x) && any(sapply(x,is.factor))){
              xx<-lapply(x, function(xi) {if (is.factor(xi)) 0+(outer(xi,levels(xi),"==")) else xi})
              cols<-sapply(xx,NCOL)
              scols<-c(0,cumsum(cols))
              cn<-character(sum(cols))
              for(i in 1:length(xx))
                  cn[scols[i]+1:cols[i]]<-paste(names(x)[i],levels(x[[i]]),sep="")
              x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
              for(i in 1:length(xx)){
                  x[,scols[i]+1:cols[i]]<-xx[[i]]
              }
              colnames(x)<-cn
          }
      }
  }
  
 
 
  x<-as.matrix(x)
 
 if (na.rm){
     nas<-rowSums(is.na(x))
     design<-design[nas==0,]
     x<-x[nas==0,,drop=FALSE]
 }
  
 wts<-design$repweights
 scale<-design$scale
 rscales<-design$rscales
 if (!is.null(rho)) .NotYetUsed("rho")
 
 if (!design$combined.weights)
     pw<-design$pweights
  else
    pw<-1
 
 rval<-colSums(design$pweights*x)

 if (is.character(deff) || deff){
     nobs<-length(design$pweights)
     npop<-sum(design$pweights)
     vsrs<-unclass(svyvar(x,design,na.rm=na.rm, return.replicates=FALSE,estimate.only=TRUE))*sum(design$pweights)^2/nobs
     if (deff!="replace")
         vsrs<-vsrs*(npop-nobs)/npop
 }

 if (getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep)){
   v<-matrix(0,nrow=NROW(rval),ncol=NROW(rval))
   repmeans<-NULL
 } else {
   if (inherits(wts, "repweights_compressed")){
     if (getOption("survey.drop.replicates") && !is.null(design$selfrep)){
       wts$index<-wts$index[!design$selfrep]
       x<-x[!design$selfrep,,drop=FALSE]
       pw<-pw[!design$selfrep]
     }
     repmeans<-matrix(ncol=NCOL(x), nrow=ncol(wts$weights))
     for(i in 1:ncol(wts$weights)){
       wi<-wts$weights[wts$index,i]
       repmeans[i,]<-t(colSums(wi*x*pw))
     }
 } else {
   if (getOption("survey.drop.replicates") && !is.null(design$selfrep)){
     wts<-wts[!design$selfrep,,drop=FALSE]
     x<-x[!design$selfrep,,drop=FALSE]
     pw<-pw[!design$selfrep]
   }
   repmeans<-matrix(ncol=NCOL(x), nrow=ncol(wts))
   for(i in 1:ncol(wts)){
     repmeans[i,]<-t(colSums(wts[,i]*x*pw))
   }
 }
   repmeans<-drop(repmeans)
   v <- svrVar(repmeans, scale, rscales,mse=design$mse,coef=rval)
 }
attr(rval,"var") <- v
attr(rval, "statistic")<-"total"
if (return.replicates){
  attr(repmeans,"scale")<-design$scale
  attr(repmeans,"rscales")<-design$rscales
  attr(repmeans,"mse")<-design$mse 
  rval<-list(mean=rval, replicates=repmeans)
}
 
if (is.character(deff) || deff)
  attr(rval,"deff") <- v/vsrs
class(rval)<-"svrepstat"
rval
}



svycoxph.svyrep.design<-function(formula, design, subset=NULL,...,return.replicates=FALSE,na.action,
                                 multicore=getOption("survey.multicore")){
  require(survival)
  subset<-substitute(subset)
  subset<-eval(subset, design$variables, parent.frame())
  if (!is.null(subset))
    design<-design[subset,]
  if (multicore && !require(parallel,quietly=TRUE))
    multicore<-FALSE
  
  data<-design$variables
  
  
  g<-match.call()
  g$design<-NULL
  g$return.replicates<-NULL
  g$weights<-quote(.survey.prob.weights)
  g[[1]]<-quote(coxph)
  g$x<-TRUE
  
  scale<-design$scale
  rscales<-design$rscales
  
  pwts<-design$pweights/sum(design$pweights)
  if (is.data.frame(pwts)) pwts<-pwts[[1]]
  
  if (!all(all.vars(formula) %in% names(data))) 
    stop("all variables must be in design= argument")
  .survey.prob.weights<-pwts
  full<-with(data,eval(g))
  
  nas<-attr(full$model, "na.action")
  
  betas<-matrix(ncol=length(coef(full)),nrow=ncol(design$repweights))
  
  wts<-design$repweights
  
  if (!design$combined.weights){
    pw1<-pwts
    rwt<-pw1/mean(pw1)
  } else{
    rwt<-1/mean(as.vector(wts[,1]))
    pw1<-rwt
  }
  
  if (length(nas))
    wts<-wts[-nas,]
  beta0<-coef(full)

  ## coxph doesn't allow zero weights
  EPSILON<-1e-10

  if(full$method %in% c("efron","breslow")){
    if (attr(full$y,"type")=="right")
      fitter<-coxph.fit
    else if(attr(full$y,"type")=="counting")
      fitter<-survival::agreg.fit
    else stop("invalid survival type")
  } else fitter<-survival::agexact.fit
                 
##  g$init<-beta0
## for(i in 1:ncol(wts)){
##    .survey.prob.weights<-as.vector(wts[,i])*pw1+EPSILON
##    betas[i,]<-with(data,coef(eval(g)))
##  }
  if (multicore){
    betas<-do.call(rbind, mclapply(1:ncol(wts), function(i){
      fitter(full$x, full$y, full$strata, full$offset,
             coef(full), coxph.control(),
             as.vector(wts[,i])*pw1+EPSILON,
             full$method,  names(full$resid))$coef
    }))
  }else{
    for(i in 1:ncol(wts)){
      betas[i,]<-fitter(full$x, full$y, full$strata, full$offset,
                        coef(full), coxph.control(),
                        as.vector(wts[,i])*pw1+EPSILON,
                        full$method,  names(full$resid))$coef
      
    }
  }
  
  if (length(nas))
    design<-design[-nas,]
  
  v<-svrVar(betas,scale, rscales, mse=design$mse, coef=beta0)
  
  full$var<-v
  if (return.replicates){
    attr(betas,"scale")<-design$scale
    attr(betas,"rscales")<-design$rscales
    attr(betas,"mse")<-design$mse
    full$replicates<-betas
  }
  full$naive.var<-NULL
  full$wald.test<-coef(full)%*%solve(full$var,coef(full))
  full$loglik<-c(NA,NA)
  full$rscore<-NULL
  full$score<-NA
  full$degf.residual<-degf(design)+1-length(coef(full)[!is.na(coef(full))])
  
  class(full)<-c("svrepcoxph","svycoxph",class(full))
  full$call<-match.call()
  full$printcall<-sys.call(-1)
  full$survey.design<-design
  
  full
}

svrepglm<-svyglm.svyrep.design<-function(formula, design, subset=NULL, ...,
                                         rho=NULL, return.replicates=FALSE, na.action,
                                         multicore=getOption("survey.multicore")){

  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svyglm")
  
  subset<-substitute(subset)
  subset<-eval(subset, design$variables, parent.frame())
  if (!is.null(subset))
    design<-design[subset,]

  if(multicore && !require("parallel",quietly=TRUE))
    multicore<-FALSE
  
  data<-design$variables
  

  g<-match.call()
  formula<-eval.parent(formula)
  environment(formula)<-environment()
  g$formula<-formula
  g$data<-quote(data)
  g$design<-NULL
  g$var<-g$rho<-g$return.replicates<-g$multicore<-NULL
  g$weights<-quote(.survey.prob.weights)
  g[[1]]<-quote(glm)      
  g$model<-TRUE
  g$x<-TRUE
  g$y<-TRUE
  
      scale<-design$scale
      rscales<-design$rscales
      if (!is.null(rho)) .NotYetUsed(rho)
      
      pwts<-design$pweights/sum(design$pweights)
      if (is.data.frame(pwts)) pwts<-pwts[[1]]
      
      if (!all(all.vars(formula) %in% names(data))) 
	stop("all variables must be in design= argument")
      .survey.prob.weights<-pwts
      full<-with(data,eval(g))
  
      full$naive.cov<-summary(full)$cov.unscaled
  
      nas<-attr(full$model, "na.action")

      if(getOption("survey.drop.replicates") && !is.null(design$selfrep) && all(design$selfrep)){

          v<-matrix(0,ncol=length(coef(full)),nrow=length(coef(full)))
          betas<-NULL

      } else {
          betas<-matrix(ncol=length(coef(full)),
                        nrow=ncol(design$repweights))
          
          if (!design$combined.weights)
              pw1<-pwts
          else
              pw1<-rep(1,length(pwts))
          wts<-design$repweights
          if (length(nas)){
              wts<-wts[-nas,]
              pw1<-pw1[-nas]
          }
          XX<-full$x
          YY<-full$y
          beta0<-coef(full)
          if (!all(is.finite(beta0))) stop(paste("Infinite/NA values in estimate (",paste(beta0,collapse=","),")"))
          if(is.null(full$offset))
          offs<-rep(0,nrow(XX))
          else
              offs<-full$offset
          incpt<-as.logical(attr(terms(full),"intercept"))
          fam<-full$family
          contrl<-full$control
          if (multicore){
            betas<-do.call(rbind,mclapply(1:ncol(wts), function(i){
              wi<-as.vector(wts[,i])*pw1
              glm.fit(XX, YY, weights = wi/sum(wi),
                      start =beta0,
                      offset = offs,
                      family = fam, control = contrl,
                      intercept = incpt)$coefficients
              
            }))
          } else {
            for(i in 1:ncol(wts)){
              wi<-as.vector(wts[,i])*pw1
              betas[i,]<-glm.fit(XX, YY, weights = wi/sum(wi),
                                 start =beta0,
                                 offset = offs,
                                 family = fam, control = contrl,
                                 intercept = incpt)$coefficients
            }
          }
          v<-svrVar(betas,scale, rscales,mse=design$mse,coef=beta0)
  }

  full$x<-NULL
  full$df.residual<-degf(design)+1-length(coef(full)[!is.na(coef(full))])
  
  if (length(nas))
      design<-design[-nas,]

  full$cov.unscaled<-v
  if (return.replicates){
    attr(betas,"scale")<-design$scale
    attr(betas,"rscales")<-design$rscales
    attr(betas,"mse")<-design$mse
    full$replicates<-betas
  }  
  class(full)<-c("svrepglm", "svyglm", class(full))
  full$call<-sys.call(-1)
  if(!("formula" %in% names(full$call))) {
    if (is.null(names(full$call)))
      i<-1
    else
      i<-min(which(names(full$call)[-1]==""))
    names(full$call)[i+1]<-"formula"
  }
  full$survey.design<-design
  full
}


print.summary.svyglm<-function (x, digits = max(3, getOption("digits") - 3),
                                symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
  ##if (!exists("printCoefmat")) printCoefmat<-print.coefmat

  cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "") 

    cat("Survey design:\n")
    print(x$survey.design$call)
   
        if (!is.null(df <- x$df) && (nsingular <- df[3] - df[1])) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- is.na(x$coefficients[,1])) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n",  "Number of Fisher Scoring iterations: ", 
        x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}



predict.svrepglm <- function(object, newdata=NULL, total=NULL,
                           type = c("link", "response","terms"),
                           se.fit=(type!="terms"),
                           vcov=FALSE, return.replicates=!is.null(object$replicates),...){
    if(is.null(newdata))
      newdata<-model.frame(object$survey.design)
    type<-match.arg(type)
    if (type=="terms")
      return(predterms(object,se=se.fit,...))
    tt<-delete.response(terms(formula(object)))
    mf<-model.frame(tt,data=newdata)
    mm<-model.matrix(tt,mf)
    if (!is.null(total) && attr(tt,"intercept")){
        mm[,attr(tt,"intercept")]<-mm[,attr(tt,"intercept")]*total
    }
    eta<-drop(mm %*% coef(object))
    d<-drop(object$family$mu.eta(eta))
    eta<-switch(type, link=eta, response=object$family$linkinv(eta))
    if(se.fit){
        if(vcov){
            vv<-mm %*% vcov(object) %*% t(mm)
            attr(eta,"var")<-switch(type,
                                    link=vv,
                                    response=d*(t(vv*d)))
        } else {
            ## FIXME make this more efficient
            vv<-drop(rowSums((mm %*% vcov(object)) * mm))
            attr(eta,"var")<-switch(type,
                                    link=vv,
                                    response=drop(d*(t(vv*d))))
        }
    }
    attr(eta,"statistic")<-type

    if (return.replicates){
      if (is.null(object$replicates)) {
        warning("replicates are not present in the fit")
      } else{
        pred.replicates<-t(apply(object$replicates,1, function(beta){
          etai<-drop(mm %*% beta)
          switch(type, link=etai, response=object$family$linkinv(etai))
        }))
        attr(pred.replicates,"scale")<-attr(object$replicates,"scale")
        attr(pred.replicates,"rscales")<-attr(object$replicates,"rscales")
        attr(pred.replicates,"mse")<-attr(object$replicates,"mse")
        eta<-list(eta,replicates=pred.replicates)
      }
    }

    class(eta)<-"svrepstat"
    eta
  }
    

    

svyratio.svyrep.design<-svrepratio<-function(numerator=formula,denominator, design,
                                             na.rm=FALSE,formula,covmat=FALSE,
                                             return.replicates=FALSE,deff=FALSE,...){

  if (!exists(".Generic"))
    .Deprecated("svyratio")
  
  if (!inherits(design, "svyrep.design")) stop("design must be a svyrepdesign object")
  
  if (inherits(numerator,"formula"))
    numerator<-model.frame(numerator,design$variables, na.action=na.pass)
  else if(typeof(numerator) %in% c("expression","symbol"))
    numerator<-eval(numerator, design$variables)
  if (inherits(denominator,"formula"))
    denominator<-model.frame(denominator,design$variables, na.action=na.pass)
  else if(typeof(denominator) %in% c("expression","symbol"))
    denominator<-eval(denominator, design$variables)
  
  nn<-NCOL(numerator)
  nd<-NCOL(denominator)
  
  all<-cbind(numerator,denominator)
  nas<-!complete.cases(all)
  if (na.rm==TRUE){
      design<-design[!nas,]
      all<-all[!nas,,drop=FALSE]
      numerator<-numerator[!nas,,drop=FALSE]
      denominator<-denominator[!nas,,drop=FALSE]
  }
  allstats<-svymean(all, design, return.replicates=TRUE)
  
  rval<-list(ratio=outer(allstats$mean[1:nn], allstats$mean[nn+1:nd], "/"))
  
  if (is.null(allstats$replicates)){
    ##only self-representing strata.
    vars<-matrix(0,nrow=nn,ncol=nd)
  }else {
    vars<-matrix(nrow=nn,ncol=nd)
    if (deff) deffs<-matrix(nrow=nn,ncol=nd)
    for(i in 1:nn){
      for(j in 1:nd){
        vars[i,j]<-svrVar(allstats$replicates[,i]/allstats$replicates[,nn+j],
                          design$scale, design$rscales,mse=design$mse,coef=rval$ratio[i,j])
        if (deff)
          deffs[i,j]<-deff(svytotal(numerator[,i]-rval[i,j]*denominator[,j],design))
      }
    }
  }
  if (covmat){
      if (is.null(allstats$replicates))
          vcovmat<-matrix(0,nn*nd,nn*nd)
      else
          vcovmat<-as.matrix(svrVar(allstats$replicates[,rep(1:nn,nd)]/allstats$replicates[,nn+rep(1:nd,each=nn)],
                          design$scale, design$rscales,mse=design$mse,coef=as.vector(rval$ratio)))
      rownames(vcovmat)<-names(numerator)[rep(1:nn,nd)]
      colnames(vcovmat)<-names(denominator)[rep(1:nd,each=nn)]
      rval$vcov<-vcovmat
  }
  if (return.replicates) {
    reps<-allstats$replicates[, rep(1:nn, nd)]/allstats$replicates[, nn + rep(1:nd, each = nn)]
    attr(reps,"scale")<-design$scale
    attr(reps,"rscales")<-design$rscales
    attr(reps,"mse")<-design$mse
    rval$replicates<-reps
  }
  rval$var<-vars
  attr(rval,"call")<-sys.call()
  if (deff) attr(rval,"deff")<-deffs
  class(rval)<-"svyratio"
  rval
}

vcov.svyratio <- function(object, ...){
    covmat<-object$vcov
    if (is.null(covmat)){
        covmat<-matrix(NaN,length(object$var),length(object$var))
        diag(covmat)<-as.vector(object$var)
    }
    nms<-as.vector(outer(rownames(object$ratio),colnames(object$ratio),paste,sep="/"))
    dimnames(covmat)<-list(nms,nms)
    covmat
}

residuals.svrepglm<-function(object,type = c("deviance", "pearson", "working", 
    "response", "partial"),...){
	type<-match.arg(type)
	if (type=="pearson"){
   	   y <- object$y
	   mu <- object$fitted.values
    	   wts <- object$prior.weights
	   r<-(y - mu) * sqrt(wts)/(sqrt(object$family$variance(mu))*sqrt(object$survey.design$pweights/sum(object$survey.design$pweights)))
	   if (is.null(object$na.action)) 
        	r
    	   else 
	        naresid(object$na.action, r)
	} else 
		NextMethod()

}

logLik.svrepglm<-function(object,...){
   stop("svrepglm not fitted by maximum likelihood.")
}


withReplicates<-function(design, theta,  ..., return.replicates=FALSE){
  UseMethod("withReplicates",design)
}

withReplicates.svrepvar<-function(design, theta, ...,return.replicates=FALSE){
  if (is.null(reps<-design$replicates)) stop("object does not contain replicate estimates")

  p<-sqrt(NCOL(reps))
    if (is.function(theta)){
        full<-theta(design[[1]],...)
        thetas<-drop(t(apply(reps,1,
                             function(rr) theta(matrix(rr,p,p), ...))))
    } else{
        full<-eval(theta, list(.replicate=design[[1]]))
        thetas<-drop(t(apply(reps,1,
                             function(rr) eval(theta, list(.replicate=matrix(rr,p,p))))))
    }

  v<-svrVar(thetas, attr(reps,"scale"), attr(reps,"rscales"), mse=attr(reps,"mse"), coef=full)
  
  attr(full,"var")<-v
  attr(full,"statistic")<-"theta"
  
  if (return.replicates){
    attr(thetas,"scale")<-attr(reps,"scale")
    attr(thetas,"rscales")<-attr(reps,"rscales")
    attr(thetas,"mse")<-attr(reps,"mse")
    rval<-list(theta=full, replicates=thetas)
  }  else {
    rval<-full
  }
  class(rval)<-"svrepstat"
  rval
  
  
}

withReplicates.svrepstat<-function(design, theta, ..., return.replicates=FALSE){
  if (is.null(reps<-design$replicates)) stop("object does not contain replicate estimates")

  reps<-as.matrix(reps)
  
  if (is.function(theta)){
        full<-theta(design[[1]],...)
        thetas<-drop(t(apply(reps,1,theta, ...)))
      } else{
        full<-eval(theta, list(.replicate=design[[1]]))
        thetas<-drop(t(apply(reps,1,
                             function(rr) eval(theta, list(.replicate=rr)))))
      }
  
  v<-svrVar(thetas, attr(reps,"scale"), attr(reps,"rscales"), mse=attr(reps,"mse"), coef=full)
  
  attr(full,"var")<-v
  attr(full,"statistic")<-"theta"
  
  if (return.replicates){
    attr(thetas,"scale")<-attr(reps,"scale")
    attr(thetas,"rscales")<-attr(reps,"rscales")
    attr(thetas,"mse")<-attr(reps,"mse")
    rval<-list(theta=full, replicates=thetas)
  }  else {
    rval<-full
  }
  class(rval)<-"svrepstat"
  rval 
}


withReplicates.svyrep.design<-function(design, theta, rho=NULL,...,
                         scale.weights=FALSE,
                         return.replicates=FALSE){
    wts<-design$repweights
    scale<-design$scale
    rscales<-design$rscales
    if (!is.null(rho)) .NotYetUsed("rho")
    
    if (scale.weights)
      pwts<-design$pweights/sum(design$pweights)
    else
      pwts<-design$pweights
    
  if (inherits(wts,"repweights_compressed")){
      if (scale.weights)
          wts$weights<-sweep(wts$weights,2,drop(colSums(wts$weights)),"/")
  } else {
      if (scale.weights)
          wts<-sweep(wts,2, drop(colSums(wts)),"/")
  }

    rpwts<-if (design$combined.weights) 1 else pwts
    data<-design$variables
    
    if (is.function(theta)){
        full<-theta(pwts,data,...)
        thetas<-drop(t(apply(wts,2,
                             function(ww) theta(as.vector(ww)*rpwts, data, ...))))
    } else{
        .weights<-pwts
        full<-with(data, eval(theta))
        thetas<-drop(t(apply(wts,2,
                             function(.weights) {.weights<-as.vector(.weights)*rpwts
                                                 with(data, eval(theta))})))
    }
      
  v<-svrVar(thetas, scale, rscales,mse=design$mse, coef=full)

    attr(full,"var")<-v
    attr(full,"statistic")<-"theta"

    if (return.replicates)
      rval<-list(theta=full, replicates=thetas)
    else
      rval<-full
    class(rval)<-"svrepstat"
    rval
  }

coef.svrepstat<-function(object,...){
  if (is.list(object)) object<-object[[1]]
  attr(object,"statistic")<-NULL
  attr(object,"deff")<-NULL
  attr(object,"var")<-NULL
  unclass(object)
}

vcov.svrepstat<-function (object, ...) 
{
  nms <- names(coef(object))
  if (is.list(object)) 
    object <- object[[1]]
  v <- as.matrix(attr(object, "var"))
  
  if (length(object) == NCOL(v)) {
    dimnames(v) <- list(nms, nms)
    v
  }
  else if (length(object) == length(v)) {
    dnms <- dimnames(coef(object))
    vmat <- matrix(NA, nrow = length(object), ncol = length(object))
    diag(vmat) <- as.vector(v)
    nms <- as.vector(outer(dnms[[1]], dnms[[2]], paste, sep = ":"))
    dimnames(vmat) <- list(nms, nms)
    vmat
  }
}



as.data.frame.svrepstat<-function(x,...){
  if (is.list(x)) {
    x<-x[[1]]
    class(x)<-"svrepstat"
  }
  rval<-data.frame(statistic=coef(x),SE=SE(x))
  names(rval)[1]<-attr(x,"statistic")
  if (!is.null(attr(x,"deff")))
    rval<-cbind(rval,deff=deff(x))
  rval
}

SE<-function(object,...){
  UseMethod("SE")
}

SE.default<-function(object,...){
  sqrt(diag(vcov(object,...)))
}

SE.svrepstat<-function(object,...){
  if (is.list(object)){
    object<-object[[1]]
  }
  vv<-as.matrix(attr(object,"var"))
  if (!is.null(dim(object)) && length(object)==length(vv))
    sqrt(vv)
  else
    sqrt(diag(vv))
}

print.svrepstat<-function(x,...){
  if (is.list(x)){
    x<-x[[1]]
  }
  vv<-attr(x,"var")
  deff<-attr(x, "deff")
  if (!is.null(dim(x)) && length(x)==length(vv)){
    cat("Statistic:\n")
    prmatrix(x)
    cat("SE:\n")
    print(sqrt(vv))
    if (!is.null(deff)){
      cat("Design Effect:\n")
      printCoefmat()
    }
  } else if(length(x)==NCOL(vv)){
    m<-cbind(x,sqrt(diag(as.matrix(vv))))
    if (is.null(deff))
      colnames(m)<-c(attr(x,"statistic"),"SE")
    else {
      m<-cbind(m,deff(x))
      colnames(m)<-c(attr(x,"statistic"),"SE","DEff")
    }
    printCoefmat(m)
  } else {stop("incorrect structure of svrepstat object")}

  naa<-attr(vv,"na.replicates")
  if (!is.null(naa))
    cat("Note: NA results discarded for",length(naa),"replicates (",naa,")\n")
}

summary.svrepglm<-function (object, correlation = FALSE, df.resid=NULL,...) 
{
    Qr <- object$qr
    est.disp <- TRUE
    if (is.null(df.resid))
      df.r <- object$df.residual
    else
      df.r<-df.resid
    presid<-resid(object,"pearson")
    dispersion<- sum(  object$survey.design$pweights*presid^2,na.rm=TRUE)/sum(object$survey.design$pweights)
    coef.p <- coef(object)
    covmat<-vcov(object)
    dimnames(covmat) <- list(names(coef.p), names(coef.p))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
    }
    else if (df.r > 0) {
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", 
            "Pr(>|t|)"))
    }
    else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p), dn)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        aic = object$aic, coefficients = coef.table, dispersion = dispersion, 
        df = c(object$rank, df.r,NCOL(Qr$qr)), cov.unscaled = covmat, 
        cov.scaled = covmat))
    if (correlation) {
        dd <- sqrt(diag(covmat))
        ans$correlation <- covmat/outer(dd, dd)
    }

    ans$aliased<-is.na(ans$coef)
    ans$survey.design<-list(call=object$survey.design$call,
                            type=object$survey.design$type)
    class(ans) <- c("summary.svyglm","summary.glm")
    return(ans)
}


svytable.svyrep.design<-svreptable<-function(formula, design,
                                             Ntotal=sum(weights(design, "sampling")),
                                             round=FALSE,...){

  if (!exists(".Generic",inherits=FALSE))
    .Deprecated("svytable")
  
   weights<-design$pweights
   if (is.data.frame(weights)) weights<-weights[[1]]
   ## unstratified or unadjusted.
   if (is.null(Ntotal) || length(Ntotal)==1){
       ff<-eval(substitute(lhs~rhs,list(lhs=quote(weights), rhs=formula[[2]])))
       tbl<-xtabs(ff, data=design$variables,...)
       if (!is.null(Ntotal)) {
           tbl<-tbl*sum(Ntotal)/sum(tbl)
       }
       if (round)
           tbl<-round(tbl)
       class(tbl) <- c("svytable", class(tbl))
       attr(tbl, "call")<-match.call()
       return(tbl)
   }
   ## adjusted and stratified
   ff<-eval(substitute(lhs~strata+rhs,list(lhs=quote(weights),
                                           rhs=formula[[2]],
                                           strata=quote(design$strata))))
   tbl<-xtabs(ff, data=design$variables,...)
   ss<-match(sort(unique(design$strata)), Ntotal[,1])
   dm<-dim(tbl)
   layer<-prod(dm[-1])
   tbl<-sweep(tbl,1,Ntotal[ss, 2]/apply(tbl,1,sum),"*")
   tbl<-apply(tbl, 2:length(dm), sum)
   if (round)
       tbl<-round(tbl)
   class(tbl)<-c("svytable", "xtabs","table")
   attr(tbl, "call")<-match.call()

   tbl
}


postStratify<-function(design,strata, population, partial=FALSE,...){
  UseMethod("postStratify")
}



postStratify.svyrep.design<-function(design, strata, population,
                                     partial=FALSE,compress=NULL,...){

  if(inherits(strata,"formula")){
    mf<-substitute(model.frame(strata, data=design$variables,na.action=na.fail))
    strata<-eval.parent(mf)
  }
  strata<-as.data.frame(strata)
  if (is.null(compress))
    compress<-inherits(design$repweights, "repweights_compressed")

  sampletable<-xtabs(design$pweights~.,data=strata)
  sampletable<-as.data.frame(sampletable)

  if (inherits(population,"table"))
    population<-as.data.frame(population)
  else if (is.data.frame(population))
    population$Freq <- as.vector(population$Freq)
  else
    stop("population must be a table or dataframe")

  if (!all(names(strata) %in% names(population)))
    stop("Stratifying variables don't match")
  nn<- names(population) %in% names(strata)
  if (sum(!nn)!=1)
    stop("stratifying variables don't match")

  names(population)[which(!nn)]<-"Pop.Freq"
  
  both<-merge(sampletable, population, by=names(strata), all=TRUE)

  samplezero <- both$Freq %in% c(0, NA)
  popzero <- both$Pop.Freq %in% c(0, NA)
  both<-both[!(samplezero & popzero),]
  
  if (any(onlysample<- popzero & !samplezero)){
    print(both[onlysample,])
    stop("Strata in sample absent from population. This Can't Happen")
  }
  if (any(onlypop <- samplezero & !popzero)){
    if (partial){
      both<-both[!onlypop,]
      warning("Some strata absent from sample: ignored")
    } else {
      print(both[onlypop,])
      stop("Some strata absent from sample: use partial=TRUE to ignore them.")
    }
  } 

  reweight<-both$Pop.Freq/both$Freq
  both$label <- do.call("interaction", list(both[,names(strata)]))
  designlabel <- do.call("interaction", strata)
  index<-match(designlabel, both$label)

  oldpw<-design$pweights
  design$pweights<-design$pweights*reweight[index]

  
  if (design$combined.weights){
    replicateFreq<- rowsum(as.matrix(design$repweights),
                           match(designlabel, both$label),
                           reorder=TRUE)
    repreweight<-  both$Pop.Freq/replicateFreq
    design$repweights <- as.matrix(design$repweights)*repreweight[index]
  } else { 
    replicateFreq<- rowsum(as.matrix(design$repweights)*oldpw,
                           match(designlabel, both$label),
                           reorder=TRUE)
    repreweight<- both$Pop.Freq/replicateFreq
    design$repweights <- as.matrix(design$repweights)* (repreweight/reweight)[index,]
  }

  if (compress) design$repweights<-compressWeights(design$repweights)
  
  design$call<-sys.call(-1)
  if(!is.null(design$degf)){
    design$degf<-NULL
    design$degf<-degf(design)
  }
  design
}


rake<-function(design, sample.margins, population.margins,
               control=list(maxit=10, epsilon=1, verbose=FALSE),
                 compress=NULL){

    if (!missing(control)){
        control.defaults<-formals(rake)$control
        for(n in names(control.defaults))
            if(!(n %in% names(control)))
                control[[n]]<-control.defaults[[n]]
    }

    is.rep<-inherits(design,"svyrep.design")

    if (is.rep && is.null(compress))
      compress<-inherits(design$repweights,"repweights_compressed")

    if (is.rep) design$degf<-NULL
    
    if (length(sample.margins)!=length(population.margins))
        stop("sample.margins and population.margins do not match.")

    nmar<-length(sample.margins)
    
    if (control$epsilon<1) 
        epsilon<-control$epsilon*sum(weights(design,"sampling"))
    else
        epsilon<-control$epsilon

    
    
    strata<-lapply(sample.margins, function(margin)
                   if(inherits(margin,"formula")){
                     mf<-model.frame(margin, data=design$variables,na.action=na.fail)
                   }
                   )
    

    allterms<-unlist(lapply(sample.margins,all.vars))
    ff<-formula(paste("~", paste(allterms,collapse="+"),sep=""))
    oldtable<-svytable(ff, design)
    if (control$verbose)
        print(oldtable)

    oldpoststrata<-design$postStrata
    iter<-0
    converged<-FALSE
    while(iter < control$maxit){
        ## we don't want to accumulate more poststrata with each iteration
        design$postStrata<-NULL
        
        for(i in 1:nmar){
            design<-postStratify(design, strata[[i]],
                                 population.margins[[i]],
                                 compress=FALSE)
        }
        newtable<-svytable(ff, design)
        if (control$verbose)
            print(newtable)

        delta<-max(abs(oldtable-newtable))
        if (delta<epsilon){
            converged<-TRUE
            break
        }
        oldtable<-newtable
        iter<-iter+1
    }

    ## changed in 3.6-3 to allow the projections to be iterated
    ## in svyrecvar
    rakestrata<-design$postStrata
    if(!is.null(rakestrata)){
      class(rakestrata)<-"raking"
      design$postStrata<-c(oldpoststrata, list(rakestrata))
    }
    
    design$call<-sys.call()

    if (is.rep && compress)
      design$repweights<-compressWeights(design$repweights)
    if(is.rep)
      design$degf<-degf(design)
    
    if(!converged)
        warning("Raking did not converge after ", iter, " iterations.\n")

    return(design)
        
}




## degrees of freedom for repweights design
degf<-function(design,...) UseMethod("degf")

degf.svyrep.design<-function(design,tol=1e-5,...){
  if (!inherits(design,"svyrep.design"))
    stop("Not a survey design with replicate weights")
  rval<-design$degf ##cached version
  if(is.null(rval))
    rval<-qr(weights(design,"analysis"), tol=1e-5)$rank-1
  rval
}

degf.survey.design2<-function(design,...){
  inset<- weights(design,"sampling")!=0
  length(unique(design$cluster[inset, 1])) - length(unique(design$strata[inset, 1]))
}

degf.twophase<-function(design,...){
  degf(design$phase2)
}

dim.svyrep.design<-function(x) dim(x$variables)
