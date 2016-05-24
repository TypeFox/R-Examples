rsnpset<-function(Y,delta=NULL,G,X=NULL,snp.sets,
                  score=c("cox", "binomial", "gaussian") ,B=0, r.method="monte carlo",
                  v.method="empirical", v.permute=FALSE, ret.rank=FALSE,
                  pinv.check=FALSE, pinv.method="specdecomp", pinv.tol=7.8e-8){              
  score <- match.arg(score)
  if(B > 0 & !(r.method %in% c("permutation", "monte carlo"))) {
    stop('r.method must be either "permutation" or "monte carlo"')
  }
  if(B > 0 & r.method=="monte carlo" & (v.permute!=FALSE | ret.rank!=FALSE)){
    message('\nNote: Setting v.permute=FALSE and ret.rank=FALSE for r.method="monte carlo".\n')
    v.permute <- FALSE
    ret.rank <- FALSE
  }
  if(score=="cox") {
    if(!all(Y>=0)) {
      stop("Y should be greater than or equal to 0")
    }
    if(!all(delta %in% c(0,1))) {
      stop("delta should be 0 or 1")
    }
  }
  delta<-as.numeric(delta) # always convert, for passing to rcpp
  
	if(score!="gaussian" & !is.null(X)) {
		stop('Covariates are only supported for score="gaussian"')
	}
  if(B > 0 & r.method=="permutation" & !is.null(X)) {
		stop("Permutation resampling is inappropriate for models including covariates")
	}
  
  if(score=="binomial" & length(unique(Y))!=2) {
    stop("Y should have two levels")
  } else if(score=="binomial" & length(unique(Y))==2) {
		Y<-as.numeric(c(0,1)[factor(Y)])
	}
  Y<-as.numeric(Y)
  if(any(is.na(Y))) {
    stop("Y should be numeric")
  }

  if(is.null(colnames(G))){
    stop("Column names of G cannot be NULL")
  }
  if(any(is.na(colnames(G)))){
    stop("Column names of G cannot be NA")
  }

	if(is.null(X)) {
		X<-matrix(0,0,0)
	} else {
    if(any(is.na(X))){ 
      stop("X should be numeric")
    }
    if(any(!is.numeric(X))){ 
      stop("X should be numeric")
    }
    if(any(apply(X, 2, sd)==0)){
      stop("X contains one or more columns of identical values")
    }
		X <- cbind(rep(1, nrow(X)), X)
	}
    
  allGeneIDs <- colnames(G)
  geneIndexSets<-vector("list",length(snp.sets))
  names(geneIndexSets)<-names(snp.sets)
  for(i in 1:length(snp.sets)) {
     geneIndexSets[[i]]<-na.omit(fmatch(snp.sets[[i]],allGeneIDs))
  }
  geneIndexSets<-geneIndexSets[lapply(geneIndexSets,length)>0]

  if (is.null(getDoParName())) {
      registerDoSEQ()
  }
  result<-foreach(i=0:B, .packages="RSNPset") %dorng%{ 
      .Call("rsnpsetRcpp",Y,delta,G,X,geneIndexSets,i,score,v.method,r.method,pinv.check,pinv.tol,PACKAGE="RSNPset")
  }
  
  attributes(result) <- NULL
  
  names(result)<-paste0("Replication.",0:(length(result)-1))
  names(result)[1]<-"Observed"
  
  class(result)<-"RSNPset"
  attr(result,"n")<-length(Y)
  attr(result,"KSub")<-length(snp.sets)
  attr(result,"KAna")<-length(geneIndexSets)
  attr(result,"mSub")<-sapply(snp.sets, length)
  attr(result,"mAna")<-sapply(geneIndexSets, length)
  attr(result,"r.method")<-r.method
  attr(result,"B")<-B
  attr(result,"ret.rank")<-ret.rank
  attr(result,"v.permute")<-v.permute
  attr(result,"pinv.tol")<-pinv.tol
  
  if(pinv.check==FALSE) {
    attr(result,"pinv.check")<-NA
  }
  else {
    attr(result,"pinv.check")<- lapply(result,function(x) x[,3:7])
  }
  
  for(i in 1:(B+1)) {
    rownames(result[[i]])<-names(geneIndexSets)
      
    if(ret.rank==FALSE && i>1) {
      result[[i]]<-result[[i]][,1,drop=FALSE]
    }
    else {
      result[[i]]<-result[[i]][,1:2]
    }
  }
  result[[1]]["m"]<-sapply(geneIndexSets, length)
  return( result )
}
