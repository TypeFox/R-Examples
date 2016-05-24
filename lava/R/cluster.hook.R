cluster.post.hook <- function(x,...) {
  if (class(x)[1]=="multigroupfit") {
    if (is.null(x$cluster)) return(NULL)
    if (any(unlist(lapply(x$cluster,is.null)))) return(NULL)
    allclusters <- unlist(x$cluster)
    uclust <- unique(allclusters)
    K <- length(uclust)
    G <- x$model$ngroup
    S0 <- lapply(score(x,indiv=TRUE), function(x) { x[which(is.na(x))] <- 0; x })
    S <- matrix(0,length(pars(x)),nrow=K)
    aS <- c() ##matrix(0,S[[1]]
    for (i in uclust) {
      for (j in seq_len(G)) {
        idx <- which(x$cluster[[j]]==i)
        if (length(idx)>0)
          S[i,] <- S[i,] + colSums(S0[[j]][idx,,drop=FALSE])
      }
    }
    J <- crossprod(S)
    I <- information(x,type="hessian",...)
    iI <- Inverse(I)
    asVar <- iI%*%J%*%iI
    x$vcov <- asVar
    return(x)
  }

  ## lvmfit:
  if (!is.null(x$cluster)) {
    uclust <- unique(x$cluster)
    K <- length(uclust)
    S <- score(x,indiv=TRUE) #,...)
    I <- information(x,type="hessian") #,...)
    iI <- Inverse(I)

    S0 <- matrix(0,ncol=ncol(S),nrow=K)
    count <- 0
    ##    J1 <- matrix(0,ncol=ncol(S),nrow=ncol(S))
    for (i in uclust) {
      count <- count+1
      S0[count,] <- colSums(S[which(x$cluster==i),,drop=FALSE])
      ##      J1 <- J1+tcrossprod(S0[count,])
    };
    p <- ncol(S)
    ## adj1 <- 1
    adj1 <- K/(K-1) 
    ## adj1 <- K/(K-p) ## Mancl & DeRouen, 2001
      
    J <- adj1*crossprod(S0)
    col3 <- sqrt(diag(iI)); ## Naive se
    nn <- c("Estimate","Robust SE", "Naive SE", "P-value")
    asVar <- iI%*%J%*%iI
  } else {
    asVar <- x$vcov
  }
  diag(asVar)[diag(asVar)==0] <- NA
  mycoef <- x$opt$estimate
  x$vcov <- asVar
  SD <- sqrt(diag(asVar))
  Z <- mycoef/SD
  pval <- 2*(pnorm(abs(Z),lower.tail=FALSE))
  if (is.null(x$cluster)) {
    col3 <- Z
    nn <-  c("Estimate","Std. Error", "Z-value", "P-value")
  }
  newcoef <- cbind(mycoef, SD, col3, pval);
  nparall <- index(x)$npar + ifelse(x$control$meanstructure, index(x)$npar.mean,0)
  if (!is.null(x$expar)) {
    nparall <- nparall+length(x$expar)
  }
  mycoef <- matrix(NA,nrow=nparall,ncol=4)
  mycoef[x$pp.idx,] <- newcoef
  colnames(mycoef) <- nn
  mynames <- c()
  if (x$control$meanstructure) {
    mynames <- vars(x)[index(x)$v1==1]
  }
  if (index(x)$npar>0) {
    mynames <- c(mynames,paste0("p",seq_len(index(x)$npar)))
  }
  if (!is.null(x$expar)) {
    mynames <- c(mynames,names(x$expar))
  }

  rownames(mycoef) <- mynames
  x$coef <- mycoef
  return(x)
}
