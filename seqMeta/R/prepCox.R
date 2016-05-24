#' @rdname prepScores
#' @export
prepCox <- function(Z, formula, SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", data=parent.frame(), verbose =FALSE){
  #require(survival)
  env <- environment()
  if(is.null(SNPInfo)){ 
    warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
    aggregateBy = "SKATgene"
  } else {
    SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
  }
  
  #fit null model: 
  nullmodel <- coxph(formula=formula,data=data)
  nullmodel$strata <- eval(parse(text=rownames(attr(nullmodel$terms, "factors"))[attr(nullmodel$terms, "specials")$strata]), envir=data) # necessary for stratified analysis - 2014-10-07 - HC
  X<-stats::model.matrix(nullmodel,data)
  rn<-row.names(stats::model.frame(nullmodel,data=data))
  nullcoef<-stats::coef(nullmodel)
  
  ##match snps in Z with master list in SNPInfo file 
  mysnps <- colnames(Z)
  
  SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
  which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
  ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])
  #which.snps <- match(mysnps[which.snps.Z],SNPInfo[,snpNames])
  
  nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
  if(nsnps == 0){ 
    stop("no column names in Z match SNP names in the SNP Info file!")
  }
  n = nrow(Z)
  
  if(verbose){
    cat("\n Calculating signed LRTs... Progress:\n")
    pb <- utils::txtProgressBar(min = 0, max = nsnps, style = 3)
    pb.i <- 0
  }
  
  ##fit individual betas/se's
  maf0 <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
  maf0[is.nan(maf0)] <- -1
  
  maf <- maf0[ZtoSI]
  names(maf) <- SNPInfo[,snpNames]
  
  #zlrt <- numeric(nrow(SNPInfo))
  zlrt <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
    if(any(is.na(z))){
      if(all(is.na(z))) z <- rep(0,length(z))
      mz <- mean(z, na.rm=TRUE)
      z[is.na(z)] <- mz
    }
    if (verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(nsnps/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))
    }
    model<- coxlr.fit(cbind(z,X), nullmodel$y, nullmodel$strata, NULL,
                      init=c(0,nullcoef),coxph.control(iter.max=100),NULL,"efron",rn)
    return(sign(stats::coef(model)[1])*sqrt(2*diff(model$loglik)))
  })[ZtoSI]
  names(zlrt) <- SNPInfo[,snpNames]
  zlrt[is.na(zlrt)] <- 0
  if(verbose) close(pb)
  
  
  #deal with monomorphic SNPs
  zlrt[maf == 0] <- 0
  
  #differentiate missing from monomorphic:
  maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1
  
  #split into genes
  zlrt 	<- 	split(zlrt, SNPInfo[,aggregateBy])
  maf 	<- 	split(maf, SNPInfo[,aggregateBy])
  
  ngenes <- length(unique(SNPInfo[,aggregateBy]))
  if(verbose){
    cat("\n Calculating covariance... Progress:\n")
    pb <- utils::txtProgressBar(min = 0, max = ngenes, style = 3)
    pb.i <- 0
  }
  
  ##get covariance matrices:
  re <- as.list(by(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
    inds <- match(snp.names,colnames(Z))
    mcov <- matrix(0,length(snp.names),length(snp.names))
    if(length(stats::na.omit(inds)) > 0){
      Z0 <- as.matrix(Z[,stats::na.omit(inds),drop=FALSE])
      if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
        if(all(is.na(z))) z <- rep(0,length(z))
        mz <- mean(z, na.rm=TRUE)
        z[is.na(z)] <- mz
        z
      })
      zvar <- apply(Z0,2,stats::var)
      mod1 <- coxlr.fit(cbind(Z0[,zvar !=0],X), nullmodel$y, nullmodel$strata, NULL,
                        init=c(rep(0,ncol(Z0[,zvar !=0,drop=FALSE])),nullcoef),coxph.control(iter.max=0),NULL,"efron",rn)
      mcov[which(!is.na(inds))[zvar !=0], which(!is.na(inds))[zvar !=0]] <- if(ncol(X) == 0) mod1$var_i
      	else mod1$var_i[1:sum(zvar !=0),1:sum(zvar !=0),drop=FALSE] - 
      	mod1$var_i[1:sum(zvar !=0),(1+sum(zvar !=0)):(ncol(X)+sum(zvar !=0)),drop=FALSE] %*% crossprod(ginv_s(
      	mod1$var_i[(1+sum(zvar !=0)):(ncol(X)+sum(zvar !=0)),(1+sum(zvar !=0)):(ncol(X)+sum(zvar !=0)),drop=FALSE]), 
      	mod1$var_i[(1+sum(zvar !=0)):(ncol(X)+sum(zvar !=0)),1:sum(zvar !=0),drop=FALSE])
    }
    rownames(mcov) <- colnames(mcov) <- snp.names
    if(verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(ngenes/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))		  
    }
    return(Matrix::forceSymmetric(Matrix(mcov,sparse=TRUE)))
  }),simplify = FALSE)
  
  ##aggregate
  for(k in 1:length(re)){
    re[[k]] <- list("scores" = zlrt[[k]]*sqrt(diag(re[[k]])), "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = 1) 
  }
  if(verbose) close(pb)
  
  class(re) <- "seqMeta"
  return(re)
}


ginv_s <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xe <- eigen(X,symmetric=TRUE)
  Positive <- Xe$values > max(tol * Xe$values[1L], 0)
  if (all(Positive)) 
    #       Xe$v %*% (1/Xe$d * t(Xe$u))
    Xe$vectors %*% (1/Xe$values * t(Xe$vectors)) # 2014-10-07 - HC
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xe$vectors[, Positive, drop = FALSE] %*% ((1/Xe$values[Positive]) * 
                                                   t(Xe$vectors[, Positive, drop = FALSE]))
}
