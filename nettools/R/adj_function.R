## Normal correlation (default pearson)
Adjcor <- function(x,method='pearson',...){
    Adj <- abs(cor(x,method=method))
    diag(Adj) <- 0
    return(Adj)
}

## WGCNA method with cor^P
AdjWGCNA <- function(x,P,...){
  ## Default P = 6
  if (P!=6)
    warning(paste("WGCNA computed with P =",P,"!!!"))
  ## Adj <- WGCNA::unsignedAdjacency(x,
  ##                                 power = P,
  ##                                 corFnc = 'cor')
  Adj <- unsignedAdjacency(x,
                           power = P,
                           corFnc = 'cor')
  diag(Adj) <- 0
  return(Adj)
}

## Function for computing the FDR value on the null hypothesis
fdrrun <- function(x,idx,FUN='cor', cl=NULL, ...){
  FUN <- match.fun(FUN)
  measure <- 1
  if (!is.na(match("measure",names(list(...))))){
    measure <- list(...)[["measure"]]
    if (length(measure)==0)
      measure <- 1
  }
  
  cormat <- matrix(Inf,ncol=dim(x)[2],nrow=dim(x)[2])
  myfun <- function(x, y, idx, tmpfun, measure){
    idi <- idx[x,]
    if (idi[1]>idi[2]){
      ## cat(idi, "\n")
      return(abs(tmpfun(sample(y[,idi[2]]),y[,idi[1]])[[measure]]))
    }else{return(NA)}
  }
  if(nrow(idx)>2){
    if (is.null(cl)){
      aa <- sapply(seq(dim(idx)[1]),myfun, y=x, idx=idx, tmpfun=FUN, measure)
    } else {
      ## NB need to load package on all the cores?!??!?!
      aa <- parSapply(cl, seq(dim(idx)[1]),myfun, y=x, idx=idx, tmpfun=FUN, measure)
    }
    idx <- cbind(idx,aa)
    tmpid <- subset(idx,!is.na(idx[,3]))
    cormat[tmpid[,c(1,2)]] <- cormat[tmpid[,c(2,1)]] <- tmpid[,3]
  } else {
    stop("Not conformable array")
  }
  invisible(cormat)
  
}

## WGCNA FDR
AdjWGCNAFDR <- function(x,FDR,P,...){
  Adj <- AdjWGCNA(x,P=P)
  idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))
  
  ## check for multicore
  if (!is.na(match("n.cores", names(list(...))))){
    n.cores <- list(...)[["n.cores"]]
    if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
      if (n.cores > detectCores())
        n.cores <- detectCores() - 1
      cl <- makeCluster(n.cores)
    } else {
      cl <- NULL
    }
  } else {
    cl <- NULL
  }
  
  for (i in seq(1/FDR)){
    if (nrow(idx) > 2){
      cormat <- fdrrun(x,idx,FUN=cor,cl)
      idx <- which(Adj>cormat,arr.ind=TRUE)
    }
  }
  adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                     dimnames=list(rownames(Adj),colnames(Adj)))
  if(dim(idx)[1]>0){
    for (i in seq(dim(idx)[1])){
      adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
    }
  }
  return(adjfinal)
}

## Bicor
Adjbicor <- function(x,...){
    ## Adj <- WGCNA::bicor(x)
    Adj <- abs(bicor(x))
    diag(Adj) <- 0
    return(Adj)
}

## Bicor FDR
AdjbicorFDR <- function(x,FDR,P,...){
    Adj <- Adjbicor(x)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))

    ## Check for multicore
    if (!is.na(match("n.cores", names(list(...))))){
      n.cores <- list(...)[["n.cores"]]
      if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
        if (n.cores > detectCores())
          n.cores <- detectCores() - 1
        cl <- makeCluster(n.cores)
      } else {
        cl <- NULL
      }
    } else {
      cl <- NULL
    }
    
    for (i in seq(1/FDR)){
      if (nrow(idx) > 2){
        cormat <- fdrrun(x,idx,FUN='bicor', cl, ...)
        idx <- which(Adj>cormat,arr.ind=TRUE)
      }
    }
    adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                       dimnames=list(rownames(Adj),colnames(Adj)))
    if(dim(idx)[1]>0){
        for (i in seq(dim(idx)[1])){
            adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
        }
    }
    return(adjfinal)
}

## TOM
AdjTOM <- function(x,P,...){
  if (P!=6)
    warning(paste("WGCNA computed with P =",P,"!!!"))
  ## Adj <- WGCNA::TOMsimilarityFromExpr(datExpr=x,
  ##                                     power=P,
  ##                                     corType="pearson",
  ##                                     TOMType="unsigned",
  ##                                     verbose=0)
  Adj <- TOMsimilarityFromExpr(datExpr=x,
                               power=P,
                               corType="pearson",
                               TOMType="unsigned",
                               verbose=0)
  diag(Adj) <- 0
  return(Adj)
}

## Aracne
AdjARACNE <- function(x,...){
    ##anet <- minet::minet(x,method="aracne",estimator="spearman")
    Adj <- minet(x,method="aracne",estimator="spearman")
    diag(Adj) <- 0
    return(Adj)
}

## CLR normalized
AdjCLR <- function(x,...){
    ## cnet <- minet::minet(x,method="clr",estimator="spearman")
    Adj <- minet(x,method="clr",estimator="spearman")
    diag(Adj) <- 0
    return(Adj)
}

## MINE
AdjMINE <- function(x,measure,alpha,C, var.thr=1e-10,...){
  Adj <- mine(x,alpha=alpha, C=C,var.thr=var.thr)[[measure]]
  if (!is.null(Adj)){
    diag(Adj) <- 0
    return(Adj)
  } else {
    stop("Invalid measure for method MINE")
  }
}

## MINEFDR
AdjMINEFDR <- function(x,measure,alpha,C,FDR,...){
    Adj <- AdjMINE(x,measure,alpha,C,...)
    idx <- as.matrix(expand.grid(seq(dim(Adj)[1]),seq(dim(Adj)[2])))

    ## check for multicore
    if (!is.na(match("n.cores", names(list(...))))){
      n.cores <- list(...)[["n.cores"]]
      if (is.numeric(n.cores) && n.cores > 1 && detectCores()>1){
        if (n.cores > detectCores())
          n.cores <- detectCores() - 1
        cl <- makeCluster(n.cores)
      } else {
        cl <- NULL
      }
    } else {
      cl <- NULL
    }
    
    for (i in seq(1/FDR)){
      if (nrow(idx) > 2){
        cormat <- fdrrun(x,idx,FUN='mine',measure=measure, cl,...)
        idx <- which(Adj>cormat,arr.ind=TRUE)
      }
    }
    adjfinal <- matrix(0,ncol=dim(Adj)[2],nrow=dim(Adj)[1],
                       dimnames=list(rownames(Adj),colnames(Adj)))
    if(nrow(idx)>2){
        for (i in seq(dim(idx)[1])){
            adjfinal[idx[i,1],idx[i,2]] <- Adj[idx[i,1],idx[i,2]]
        }
    }
    return(adjfinal)
}


## DTWMIC
AdjDTWMIC <- function(x,DP,...){
    Adj <- matrix(0,nrow=ncol(x),ncol=ncol(x))
    for(i in 1:(ncol(x)-1)){
        for(j in (i+1):ncol(x)){
            d <- 1/(1+dtw(x[,i],x[,j],distance.only=TRUE)$normalizedDistance)
            m <- mine(x[,i],x[,j])$MIC
            Adj[i,j] <- Adj[j,i]<- sqrt(d**2/2+m**2/2)**DP
        }
    }
    rownames(Adj) <- colnames(Adj) <- colnames(x)
    return(Adj)
}

## Function for check the variance by features
checkvar <- function(x, tol=1e-5, ...){
  
  ## Compute the variance by columns
  feat.var <- apply(x,2,var)
  
  ## Get the indexes of the feature with low variance
  idx <- which(feat.var<tol)
  if (length(idx) == 0L){
    return (NULL)
  } else {
    return(idx)
  }
}
