##
## Constructing cov(R_i,R_j)/pi^*_ij, or \check{\check{\Delta}}_ij in Sarndal's notation
## We use this form because it can be sparse and because it is easy to combine
## multistage and multiphase sampling.
##
## The routines to compute the variances are in ht.R
##

pi2Dcheck<-function(pmat,tolerance=min(pmat)/1e4){
    rval<-(pmat-outer(diag(pmat),diag(pmat)))/pmat
    rval[abs(rval)<tolerance]<-0
    as(rval,"sparseMatrix")
}


## Overton's approximation for PPS
overton2Dcheck<-function(prob,strat=rep(1,length(prob))){
    fbar<-outer(prob,prob,"+")/2
    n<-ave(strat,strat,FUN=length)
    rval<- 1- (n-fbar)/(n-1) 
    rval[!outer(strat,strat,"==") | fbar==1]<-0
    diag(rval)<-(1-diag(fbar))
    as(rval,"sparseMatrix")
}

multi.overton2Dcheck<-function(id,strata,prob){
    nstage<-ncol(id)
    rval<-vector("list",nstage)
    for(stage in 1:nstage){
        uid<-!duplicated(id[,stage])
        rval[[stage]]<-list(id=id[,stage],
                            dcheck=overton2Dcheck(prob[uid,stage],
                              strata[uid,stage])
                           )
      }
    rval
}

## truncated Hartley-Rao approximation
HRDcheck<-function(prob,strat=rep(1,length(prob)),p2bar){
    fbar<-outer(prob,prob,"+")
    n<-ave(strat,strat,FUN=length)
    rval<- 1- (n-fbar+p2bar)/(n-1) 
    rval[!outer(strat,strat,"==") | fbar==1]<-0
    diag(rval)<-(1-prob)
    as(rval,"sparseMatrix")
}

multi.HRDcheck<-function(id,strata,prob,p2bar){
    nstage<-ncol(id)
    rval<-vector("list",nstage)
    for(stage in 1:nstage){
        uid<-!duplicated(id[,stage])
        rval[[stage]]<-list(id=id[,stage],
                            dcheck=HRDcheck(prob[uid,stage],
                              strata[uid,stage],p2bar[[stage]][strata[uid,stage]])
                           )
      }
    rval
}

## truncated Hartley-Rao approximation, using sample estimate mean(p) for sum(p^2/n)

multi.HR1Dcheck<-function(id,strata,prob){
    nstage<-ncol(id)
    rval<-vector("list",nstage)
    for(stage in 1:nstage){
        uid<-!duplicated(id[,stage])
        rval[[stage]]<-list(id=id[,stage],
                            dcheck=HRDcheck(prob[uid,stage],
                              strata[uid,stage],
                              ave(prob[uid,stage], strata[uid,stage])
                              )
                            )
      }
    rval
  }

##not used yet
combine_stages<-function(Dcheck1,Dcheck2){
  as(-Dcheck1*Dcheck2+Dcheck1+Dcheck2,"sparseMatrix")
}

make_pps_covmat<-function(design,method){##FIXME
  require("Matrix",quietly=TRUE) || stop("These designs require the Matrix package")
  if (method=="overton")
    multi.overton2Dcheck(design$cluster, design$strata, design$allprob)
  else stop("method",method,"not recognized")

}

image.pps<-function(x,...){
  require(Matrix)
  Matrix::image(x$dcheck[[1]]$dcheck,...)
}

##
pps_design<-function(method, ids,strata=NULL, probs=NULL, fpc=NULL,
                     subset, data,call=sys.call(),variance="HT",...){
  UseMethod("pps_design")
}
pps_design.character<-function(method,ids,strata=NULL, probs=NULL, fpc=NULL, variables=variables,
                   subset, data,call=sys.call(),variance="HT",...){

  if (length(ids[[2]])>1 && method!="brewer") stop("Multistage PPS sampling not supported with this method")
  rval<-svydesign(ids=ids,strata=strata,weights=NULL,
                probs=probs,fpc=fpc,data=data,pps="other")

  deltacheck<-make_pps_covmat(rval, method)

  rval$dcheck=deltacheck
  rval$variance<-variance

  rval$call<-call
  class(rval) <- c("pps","survey.design")
  rval
}

ppsmat<-function(jointprob, tolerance=0.0001){
  if ((!is.matrix(jointprob)) || !(NROW(jointprob)==NCOL(jointprob)))
    stop("jointprob must be a square matrix")
  rval<-list(pij=jointprob, tolerance=tolerance,call=sys.call())
  class(rval)<-"ppsmat"
  rval
}

print.ppsmat<-function(x,...) {
  cat("PPS: Joint probability matrix: ")
  print(x$call)
  invisible(x)
}

pps_design.ppsmat<-function(method,ids,strata=NULL, probs=NULL, fpc=NULL,variables=variables,
                   subset, data,call=sys.call(),variance="HT",...){

  if (length(ids[[2]])>1) stop("Multistage PPS sampling not supported")
  rval<-svydesign(ids=ids,strata=strata,weights=NULL,variables=variables,
                probs=probs,fpc=fpc,data=data,pps="other")

  deltacheck<-pi2Dcheck(method$pij,method$tolerance)
  rval$variance<-variance

  rval$dcheck<-list(list(id=1:nrow(method$pij), dcheck=deltacheck))

  rval$call<-call
  class(rval) <- c("pps","survey.design")
  rval
}

HR<-function(psum=NULL, strata=NULL){
  if (is.null(psum)) { ## estimate
    rval<-list(pbar=NULL,call=sys.call())
  } else if (is.data.frame(strata) || is.matrix(strata)){ #redundant
    pbar<-lapply(1:NCOL(strata), function(i){
      psum[!duplicated(strata[,i]),i]})
    strata<-lapply(1:NCOL(strata), function(i){
      strata[!duplicated(strata[,i]),i]})
    rval<-list(pbar=pbar, strata=strata,call=sys.call())
  } else if (is.null(strata) && is.numeric(psum) && length(psum)==1){
    ## single number
    rval<-list(pbar=list(psum),strata=list(1), call=sys.call())
  } else{ ## non-redundant list
    rval<-list(pbar=psum, strata=strata,call=sys.call())
  }
  class(rval)<-"HR"
  rval
}

print.HR<-function(x,...) {
  cat("PPS: Hartley-Rao correction: ")
  print(x$call)
  invisible(x)
}
pps_design.HR<-function(method,ids,strata=NULL, probs=NULL, fpc=NULL,
                   subset, data,call=sys.call(),variables=variables,variance="HT",...){

  if (length(ids[[2]])>1) stop("Multistage PPS sampling not supported with this method")
  rval<-svydesign(ids=ids,strata=strata,weights=NULL,
                probs=probs,fpc=fpc,data=data,pps="other")

  if (is.null(method$pbar))    ## sample estimate of sum(p^2/n)
    deltacheck<-multi.HR1Dcheck(rval$cluster,rval$strata,rval$allprob)
  else 
    deltacheck<-multi.HRDcheck(rval$cluster,rval$strata,rval$allprob, method$pbar)

  rval$dcheck=deltacheck
  rval$variance<-variance

  rval$call<-call
  class(rval) <- c("pps","survey.design")
  rval
}
print.pps<-function(x,...){
  cat("Sparse-matrix design object:\n ")
  print(x$call)
}

summary.pps<-function(object,...){
  class(object)<-"summary.pps"
  object
}

print.summary.pps<-function(x,...,varnames=TRUE){
  cat("Two-phase sparse-matrix design:\n ")
  print(x$call)
  cat("Sampling probabilities:\n")
  print(summary(x$prob))
  if (varnames){
    cat("Data variables:\n")
    print(names(x$variables))
  }
  invisible(x)
}




ppsvar<-function(x,design){
  postStrata<-design$postStrata
  est<-design$variance ##Yates-Grundy or Horvitz-Thompson
  if (!is.null(postStrata)){
    for (psvar in postStrata){
      if (inherits(psvar, "greg_calibration")) {
        if (psvar$stage==0){
          ## G-calibration at population level
          y<-qr.resid(psvar$qr,y/psvar$w)*psvar$w
        } else {
          ## G-calibration within clusters
          stop("calibration within clusters not yet available for PPS designs")
        }
      } else {
        ## ordinary post-stratification
        psw<-attr(psvar, "weights")
        postStrata<-as.factor(psvar)
        psmeans<-rowsum(y/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
        x<- y-psmeans[match(psvar,sort(unique(psvar))),]*psw
      }
    }
  }
  dcheck<-design$dcheck
  if (length(dcheck)!=1) stop("Multistage not implemented yet")
  rval<-switch(est,HT=htvar.matrix(rowsum(x,dcheck[[1]]$id,reorder=FALSE),dcheck[[1]]$dcheck),
               YG=ygvar.matrix(rowsum(x,dcheck[[1]]$id,reorder=FALSE),dcheck[[1]]$dcheck),
               stop("can't happen"))
  rval
}

svytotal.pps<-function(x,design, na.rm=FALSE, deff=FALSE,...){
    
    
    if (inherits(x,"formula")){
        ## do the right thing with factors
        mf<-model.frame(x,model.frame(design), na.action=na.pass)
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
    if(length(nas)>length(design$prob))
      x<-x[nas==0,,drop=FALSE]
    else
      x[nas>0,]<-0
  }
    
    N<-sum(1/design$prob)
    total <- colSums(x/as.vector(design$prob),na.rm=na.rm)
    class(total)<-"svystat"
    attr(total, "var")<-v<-ppsvar(x/design$prob,design)
    attr(total,"statistic")<-"total"
    
    if (is.character(deff) || deff){
      nobs<-NROW(design$cluster)
      if (deff=="replace")
        vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2*(N-nobs)/N
      else
        vsrs<-svyvar(x,design,na.rm=na.rm)*sum(weights(design))^2
    attr(total, "deff")<-v/vsrs
    }
    
    
  return(total)
  }

svymean.pps<-function(x,design, na.rm=FALSE,deff=FALSE,...){
  
  if (inherits(x,"formula")){
    ## do the right thing with factors
    mf<-model.frame(x,model.frame(design) ,na.action=na.pass)
    xx<-lapply(attr(terms(x),"variables")[-1],
               function(tt) model.matrix(eval(bquote(~0+.(tt))),mf))
    cols<-sapply(xx,NCOL)
    x<-matrix(nrow=NROW(xx[[1]]),ncol=sum(cols))
    scols<-c(0,cumsum(cols))
    for(i in 1:length(xx)){
      x[,scols[i]+1:cols[i]]<-xx[[i]]
    }
    colnames(x)<-do.call("c",lapply(xx,colnames))
  }
  else {
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
    if (any(nas>0))
      design<-design[nas==0,]
    x[nas>0,]<-0
  }
  
  pweights<-1/design$prob
  psum<-sum(pweights)
  average<-colSums(x*pweights/psum)
  x<-sweep(x,2,average)
  v<-ppsvar(x*pweights/psum,design)
  attr(average,"var")<-v
  attr(average,"statistic")<-"mean"
  class(average)<-"svystat"
  if (is.character(deff) || deff){
      nobs<-nrow(design)
      if(deff=="replace"){
        vsrs<-svyvar(x,design,na.rm=na.rm)/(nobs)
      } else {
        if(psum<nobs) {
          vsrs<-NA*v
          warning("Sample size greater than population size: are weights correctly scaled?")
        } else{
          vsrs<-svyvar(x,design,na.rm=na.rm)*(psum-nobs)/(psum*nobs)
        }
      }
      attr(average, "deff")<-v/vsrs
  }
  
  return(average)
}



svyratio.pps<-function(numerator=formula, denominator, design, separate=FALSE,na.rm=FALSE,formula,...){

    if (separate){
      strats<-sort(unique(design$strata[,1]))
      rval<-list(ratios=lapply(strats,
                   function(s) {
                     tmp<-svyratio(numerator, denominator,
                                   subset(design, design$phase2$strata[,1] %in% s),
                                   separate=FALSE,...)
                     attr(tmp,"call")<-bquote(Stratum==.(s))
                     tmp}))
      names(rval$ratios)<-strats
   
      class(rval)<-c("svyratio_separate")
      rval$call<-sys.call()
      rval$strata<-strats
      return(rval)
    }
  
    if (inherits(numerator,"formula"))
        numerator<-model.frame(numerator,model.frame(design),na.action=na.pass)
    else if(typeof(numerator) %in% c("expression","symbol"))
        numerator<-eval(numerator, design$variables)
    if (inherits(denominator,"formula"))
        denominator<-model.frame(denominator,model.frame(design),na.action=na.pass)
    else if(typeof(denominator) %in% c("expression","symbol"))
        denominator<-eval(denominator, model.frame(design))

    nn<-NCOL(numerator)
    nd<-NCOL(denominator)

    all<-cbind(numerator,denominator)
    nas<-!complete.cases(all)
    if (na.rm){
        design<-design[!nas,]
        all<-all[!nas,,drop=FALSE]
        numerator<-numerator[!nas,,drop=FALSE]
        denominator<-denominator[!nas,,drop=FALSE]
    }
    allstats<-svytotal(all, design) 
    rval<-list(ratio=outer(allstats[1:nn],allstats[nn+1:nd],"/"))


    vars<-matrix(ncol=nd,nrow=nn)
    for(i in 1:nn){
      for(j in 1:nd){
        r<-(numerator[,i]-rval$ratio[i,j]*denominator[,j])/sum(denominator[,j]/design$prob)
        vars[i,j]<-ppsvar(r*1/design$prob, design)
      }
    }
    colnames(vars)<-names(denominator)
    rownames(vars)<-names(numerator)
    rval$var<-vars
    attr(rval,"call")<-sys.call()
    class(rval)<-"svyratio"
    rval
    
  }


"[.pps"<-function (x,i, ..., drop=TRUE){
  if (!missing(i)){ 
    ## Set weights to zero:  don't try to save memory
    ## There should be an easier way to complement a subscript..
    if (is.logical(i) && any(!i)){
      ## logical indexing: use !
      x$prob[!i]<-Inf
      x$dcheck<-lapply(x$dcheck, function(m) {m$dcheck[!i,!i]<-0; m})
    } else if (is.numeric(i) && length(i)){
      ## numeric indexing: use -
      x$prob[-i]<-Inf
      x$dcheck<-lapply(x$dcheck, function(m) {m$dcheck[-i,-i]<-0;m})
    } else if (is.character(i)) {
      ##character indexing: use brute force and ignorance
      tmp<-x$prob[i,]
      x$prob<-rep(Inf, length(x$prob))
      x$prob[i,]<-tmp
      x$dcheck<-lapply(x$dcheck, function(m) {n<-Matrix(ncol(m$dcheck),ncol(m$dcheck)); n[i,i]<-m$dcheck[i,i]; m$dcheck<-n;m})
    }
    index<-is.finite(x$prob)
    psu<-!duplicated(x$cluster[index,1])
    tt<-table(x$strata[index,1][psu])
    if(any(tt==1)){
      warning(sum(tt==1)," strata have only one PSU in this subset.")
    }
  } else {
    
  }
  x
}

degf.pps<-function(design,...) {
    inset <- weights(design, "sampling") != 0
    length(unique(design$cluster[inset, 1])) - length(unique(design$strata[inset,1]))
}





postStratify.pps<-function(design, ...) {
	stop("postStratify not yet implemented for these pps designs. Use calibrate()")
}
