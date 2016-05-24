cvobs <- function(n,ntest,ntrain,Kfold,type= c("random","timeseries","consecutive", "interleaved"),npermut,seed) {
  tabletype <- c("random","timeseries","consecutive", "interleaved")
  type <- match.arg(type,tabletype)
  if (missing(ntest)) ntest <- NULL
  if (missing(ntrain)) ntrain <- NULL
  if (missing(npermut)) npermut <- NULL
  if (missing(Kfold)) Kfold <- FALSE
  if (is.logical(Kfold)&&(!Kfold)) {
    if (type%in%c("interleaved","consecutive")) stop("When NOT using Kfold cross-validation, type must be 'random' or 'timeseries'\n")
  }
  if (missing(seed)) seed <- NULL
  if (type=="timeseries") {
    if (is.null(ntest)&is.null(ntrain)) 
      stop("argument ntest or ntrain is missing\n")
    if (!is.null(ntest)) {
      if (is.null(npermut)) stop("number of permutation is missing\n")
      res <- as.list(rep(0,npermut))
      names(res) <- paste("V",1:npermut,sep="")
      for (j in 1:npermut) {
        res[[j]] <- seq(n-ntest+1-npermut+j,n-npermut+j,by=1)
      }
      attr(res, "type") <- "timeseries"
      return(res)
    }
    if (!is.null(ntrain)) {
      ntest <- n-ntrain-npermut+1
      if (is.null(npermut)) stop("number of permutation is missing\n")
      res <- as.list(rep(0,npermut))
      names(res) <- paste("V",1:npermut,sep="")
      for (j in 1:npermut) {
        res[[j]] <- seq(n-ntest+1-npermut+j,n-npermut+j,by=1)
      }
      attr(res, "timeseries")
      attr(res,"incomplete") <- 0
      return(res)
    }
  } else {
  if (is.logical(Kfold)) {
    if ((Kfold)&(is.null(ntest)&is.null(ntrain)))
      stop("argument(s) Kfold and/or ntest and/or ntrain are missing\n")
    if (Kfold&(!is.null(ntest))) {
      if (!is.numeric(ntest)) stop("ntest is not numeric\n")
      Kfold <- n%/%ntest
    }
    if (Kfold&(!is.null(ntrain))) {
      if (!is.numeric(ntrain)) stop("ntrain is not numeric\n")
      Kfold <- n%/%(n-ntrain)
    }
  }
  if (is.numeric(Kfold))   {
    if (n < Kfold) 
      stop("More segments (Kfold) than observations")
    length.seg <- ceiling(n/Kfold)
    incomplete <- Kfold * length.seg - n
    complete <- Kfold - incomplete  
    type <- match.arg(type)
    if (!is.null(seed)) set.seed(seed)
    switch(type, random = {
      inds <- matrix(c(sample(1:n), rep(NA, incomplete)), nrow = length.seg, 
                     byrow = TRUE)
    }, consecutive = {
      if (complete < Kfold) {
        inds <- cbind(matrix(1:(length.seg * complete), nrow = length.seg), 
                      rbind(matrix((length.seg * complete + 1):n, nrow = length.seg - 
                                   1), NA))
      }
      else {
        inds <- matrix(1:n, nrow = length.seg)
      }
    }, interleaved = {
      inds <- matrix(c(1:n, rep(NA, incomplete)), nrow = length.seg, 
                     byrow = TRUE)
    })
    res <- lapply(as.data.frame(inds), function(x) c(na.omit(x)))
    attr(res, "incomplete") <- incomplete
    attr(res, "type") <- if (length.seg == 1) 
      "leave-one-out"
    else type
  }
  if ((!is.numeric(Kfold))&(type=="random")) {
    if (is.null(npermut)) stop("number of random draw of test set (npermut) is unknown\n")
    if (!is.null(seed)) set.seed(seed)
    if (is.null(ntest)&is.null(ntrain)) stop("ntrain or ntest is missing\n")
    if (is.null(ntest)) ntest <- n-ntrain
    if (is.null(ntrain)) ntrain <- n-ntest
    if (n!= (ntrain+ntest)) stop("number of observations not equal the sum of those  in test set plus those in test set\n")
    res <- as.list(rep(0,npermut))
    names(res) <- paste("V",1:npermut,sep="")
    for (j in 1:npermut) {
      res[[j]] <- sample(1:n,ntest)
    }
    attr(res, "type") <- "randperm"
    attr(res,"incomplete") <- 0
  }
}
  return(res)
}

iterchoiceAcv <- function(X,y,bx,df,kernelx,ddlmini,ntest,ntrain,Kfold,type,npermut,seed,Kmin,Kmax,criterion,fraction) {
  choixssecv2 <- function(k,sel,SSx,y,valpr,tPADmdemiY,DdemiPA,ddlmin,index0) {
    sse <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
      }
      prov <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      prov1 <- matrix(prov*as.vector(tPADmdemiY[[j]]),nj,1)
      Yrescv <- SSx[[j]]%*%DdemiPA[[j]]%*%prov1
      sse <- sse+sum((y[sel[[j]]]-Yrescv)^2)
    }
    return(sse)
  }
  choixsapcv2 <- function(k,sel,SSx,y,valpr,tPADmdemiY,DdemiPA,ddlmin,index0) {
    sap <- 0
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        nj <- length(y[-(sel[[j]][1]:n)])
      } else {
        nj <- length(y[-sel[[j]]])
      }
      prov <- rev(sumvalpr(k,nj,rev(valpr[[j]]),nj-index0[j]+1,nj-ddlmin[j]+1))
      prov1 <- matrix(prov*as.vector(tPADmdemiY[[j]]),length(y[-sel[[j]]]),1)
      Yrescv <- SSx[[j]]%*%DdemiPA[[j]]%*%prov1
      sap <- sap+sum(abs((y[sel[[j]]]-Yrescv)/y[sel[[j]]]))
    }
    return(sap)
  }
  
  n <- nrow(X)
  sel <- cvobs(n,ntest,ntrain,Kfold,type,npermut,seed)
  res <- list(minimum=0,objective=.Machine$double.xmax)
  objectif <- .Machine$double.xmax
  fraction <- sort(unique(c(fraction[(fraction<Kmax)&(fraction>Kmin)],Kmin,Kmax)))
  dep <- length(fraction)
    tPADmdemiY <- as.list(rep(0,length(sel)))
    valpr <- DdemiPA <- SSx <- tPADmdemiY
    index0 <- ddlmin <- nj <- rep(0,length(sel))
    for (j in 1:length(sel)) {
      if (attr(sel,"type")=="timeseries") {
        XA <- X[-(sel[[j]][1]:n),,drop=FALSE]
        YA <- y[-(sel[[j]][1]:n)]
      } else  {
        XA <- X[-sel[[j]],,drop=FALSE]
        YA <- y[-sel[[j]]]
      }  
      if (is.null(bx)&(!is.null(df))) {
        bx <- bwchoice(X=XA,objectif=df,kernelx,itermax=100)
      }
      nj[j] <- length(YA)
      listeA <- calcA(X=XA,bx=bx,kernelx=kernelx)
      listeA.eig <- eigen(listeA$A,symmetric=TRUE)
      if (any(zapsmall(listeA.eig$values-1,digits=9)==0)) {
        ddlmin[j] <-  sum(zapsmall(listeA.eig$values-1,digits=9)==0)
      } else ddlmin[j] <- NA
      if (any(zapsmall(listeA.eig$values,digits=9)==0)) {
        index0[j] <-  which(zapsmall(listeA.eig$values,digits=9)==0)[1]
      } else index0[j] <- NA
      tPADmdemiY[[j]] <- t(listeA.eig$vectors*(1/listeA$Ddemi))%*%YA
      valpr[[j]]  <- 1-listeA.eig$values
      DdemiPA[[j]] <- (listeA$Ddemi*listeA.eig$vectors)
      SSx[[j]] <- kernelSx(kernelx=kernelx,X=XA,bx=bx,X[sel[[j]],,drop=FALSE])
    }
    if (criterion=="rmse") {
      repeat {
        mini <- fraction[dep-1]
        objectif <- choixssecv2(mini,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (objectif<.Machine$double.xmax) {
          bb <- fraction[dep]
          break 
        }
        dep <- dep-1
        if (dep==1) stop(paste("decrease Kmax below",fraction[2]))
      }
      phi <- (sqrt(5) - 1)/2
      repeat {
        x1 <- mini + (1-phi)*(bb-mini)
        objectif <- choixssecv2(x1,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (objectif<.Machine$double.xmax) break else {
          bb <- x1
        }
      }
      fraction <- c(fraction[1:(dep-1)],bb)
      for (i in 1:(length(fraction)-1)) {
        res1 <- optimize(choixssecv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (res1$objective<res$objective) res <- res1
      }
    } else {
      repeat {
        mini <- fraction[dep-1]
        objectif <- choixsapcv2(mini,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (objectif<.Machine$double.xmax) {
          bb <- fraction[dep]
          break 
        }
        dep <- dep-1
        if (dep==1) stop(paste("decrease Kmax below",fraction[2]))
      }
      phi <- (sqrt(5) - 1)/2
      repeat {
        x1 <- mini + (1-phi)*(bb-mini)
        objectif <- choixsapcv2(x1,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (objectif<.Machine$double.xmax) break else {
          bb <- x1
        }
      }
      fraction <- c(fraction[1:(dep-1)],bb)
      for (i in 1:(length(fraction)-1)) {
        res1 <- optimize(choixsapcv2,lower=fraction[i],upper=fraction[i+1],tol=0.5,sel=sel,SSx=SSx,y=y,valpr=valpr,tPADmdemiY=tPADmdemiY,DdemiPA=DdemiPA,ddlmin=ddlmin,index0=index0)
        if (res1$objective<res$objective) res <- res1
      }
    }
  return(list(iter=round(res$minimum),objective=res$objective/sum(unlist(lapply(sel,length)))))
}


