netSI <- function(x,indicator="all", d='HIM', adj.method='cor', 
                  method="montecarlo", k=3, h=20, n.cores=NULL,save=FALSE,
                  verbose=TRUE, ...){

  ## Get the function call parameters
  ## NB the parameters should be evaluated with eval()
  Call <- match.call()
  
  ##add a check so that an unexisting parameter cannot be passed
  id.Call <- match( names(Call),c("x", "indicator", "d", "adj.method" , 
                                  "method","k","h","n.cores","save","verbose",
                                  "FDR","P","measure","alpha","C","DP","var.thr",
                                  "components","n.nodes","ga","sseed","rho"), nomatch=0)
  if(sum(id.Call[-1]==0)==1){
    warning("The parameter '",names(Call)[which(id.Call==0)[2]],"' will be ignored",call.=FALSE)
  }
  if(sum(id.Call[-1]==0)>1){
    msg <- "The following parameters will be ignored:\n"
    for(i in which(id.Call==0)[-1]){
      msg <- paste(msg,"'",names(Call)[i],"'\n",sep="")
    }
    warning(msg,call.=FALSE)
  }
  
  id.Call <- match(c("x", "indicator", "d", "adj.method","method","k","h","n.cores"), 
                   names(Call), nomatch=0)
  if(id.Call[1]==0){
    stop("A dataset must be provided",call.=FALSE)
  }else{
    ## m <- eval(temp, parent.frame())
    ## x <- eval(Call$x)
    if(!(is.matrix(x) | is.data.frame(x))){
      stop("x must be a matrix or a data frame", call.=FALSE)
    }
  }
  
  ## Choose the indicators
  if(id.Call[2]!=0){
    indicator <- eval(Call$indicator)
    if(!is.character(indicator))
      stop("indicator must be one of 'S','SI','Sw','Sd','all'.", call.=FALSE)
  }
  
  INDICATORS <- c('S','SI','Sw','Sd',"all")
  indicator <- pmatch(indicator,INDICATORS)
  
  if(is.na(indicator))
    stop("invalid indicator", call. =FALSE)
  if(indicator == -1)
    stop("ambiguous indicator", call. =FALSE)
  
  DISTANCE <- c("HIM","IM","H",
                "ipsen","Ipsen","IpsenMikhailov","Ipsen-Mikhailov",
                "hamming","Hamming")
  
  d <- pmatch(d, DISTANCE)
  if(d==2L | (d>=4  & d<=7L))
    d <- 2L
  if(d==3L |(d>=8L & d<=9L))
    d <- 3L
  
  if(is.na(d))
    stop("invalid distance", call. =FALSE)
  if(d == -1)
    stop("ambiguous distance", call. =FALSE)
  
  d <- c("HIM","IM","H")[d]
  
  ##check on save and verbose
  if(!is.logical(save))
    stop("save must be TRUE or FALSE",call.=FALSE)
  
  if(!is.logical(verbose))
      stop("verbose must be TRUE or FALSE", call.=FALSE)
  
  ## It can be passed through ...
  if (is.null(Call$sseed)){
    sseed <- 0
  } else {
    sseed <- eval(Call$sseed)}
  set.seed(sseed)
  
  ## Pass parameter gamma to netdist functions
  ## if(!is.null(Call$ga)){
  ##   ga <- eval(Call$ga)
  ## } else {
  ##   ga <- Call$ga}

  ## Pass parameter components to netdist functions
  if(!is.null(Call$components)){
    components <- eval(Call$components)
    if(d=="HIM" & components==TRUE){
      warning(paste("components parameter will be ignored. \n",
                    "The stability indicators will be computed just for",
                    d, "distance.\n",
                    "For computing them for Hamming or Ipsen-Mikhailov", 
                    "distance use dist=H or dist=IM"), call.=FALSE)
      components <- FALSE
    }
  }
  
  ## n.cores parameter
  ## n.cores <- eval(Call$n.cores)
  
  ## Check availability of cores, otherwise set a default
  if(is.null(n.cores)){
    if(detectCores()>=2){
      n.cores <- detectCores()-1
      cl <- makeCluster(n.cores)
      warning("The computation has been automatically parallelized", call.=FALSE)
    } else {
      cl <- NULL
    }
  } else {
    if (n.cores==1){
      cl <- NULL
    } else {
      if (n.cores<detectCores()){
        cl <- makeCluster(n.cores)
      } else {
        if(detectCores()>=2){
          n.cores <- detectCores()-1
          cl <- makeCluster(n.cores)
          warning("The computation has been automatically parallelized", call.=FALSE)
        } 
      }
    }
  }
  
  ## Get the dimension of the input matrix
  ddim <- nrow(x)
  
  ## Get the resampling indexes
  if(verbose==TRUE) cat("computing resampling...\n")
  idxs <-  resamplingIDX(ddim,method=method, k=k, h=h)

  ## length of the list for optimization purposes
  ADJcv <- vector("list",length=length(idxs))
  
  if(verbose==TRUE) cat("computing adjacency matrices...\n")
  
  ## Compute the adjacency matrices for each resampling index
  if(!is.null(cl)){
    ## Parallel computation
    ADJcv <- parLapply(cl=cl,X=idxs,fun=function(x,DAT,method,...){
      ss <- DAT[x,]
      tmp <- mat2adj(ss, method=method, ...)
      return(tmp)
    },DAT=x,method=adj.method,...)
  } else {
    ## One core computation
    ADJcv <- lapply(X=idxs,FUN=function(x,DAT,method, ...){
      ss <- DAT[x,]
      tmp <- mat2adj(ss, method=method, ...)
      return(tmp)
    },DAT=x,method=adj.method,...)
  }
  
  ##computing the adjacency matrix on the whole dataset
  ADJall <- mat2adj(x=x,method=adj.method,...)
  
  ## Compute the stability indicators
  netsi <- list()
  if(indicator==1L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator S...\n")
    netsi[["S"]] <- netsiS(ADJall, ADJcv, d=d, cl=cl, ...)
  }
  if(indicator==3L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sw...\n")
    netsi[["Sw"]] <- netsiSw(ADJcv, cl=cl)
  }
  if(indicator==4L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator Sd...\n")
    netsi[["Sd"]] <- netsiSd(ADJcv, cl=cl)
  }
  if (!is.null(cl))
    stopCluster(cl)

  if(indicator==2L | indicator==5L){
    if(verbose==TRUE) cat("computing stability indicator SI...\n")
    netsi[["SI"]] <- netsiSI(ADJcv, d=d, n.cores=n.cores, ...)
  }

  results <- list("S"=mean(netsi[["S"]]),
                  "SI"=mean(netsi[["SI"]]),
                  "Sw"=apply(netsi[["Sw"]], 1,compute.indicator),
                  "Sd"=apply(netsi[["Sd"]], 2,compute.indicator)
                  )
  
  if(save==TRUE){
    results$call <- Call
    results$ADJ <- ADJall
    results$ADJlist <- ADJcv
    results$S_boot <- netsi[["S"]]
    results$SI_boot <- netsi[["SI"]]
    results$Sw_boot <- netsi[["Sw"]]
    results$Sd_boot <- netsi[["Sd"]]
  }
  
  return(results)
}


## Stability indicator S
netsiS <- function(g, H, d, cl, ...){
  DIST <- c("HIM","IM","H")
  type <- pmatch(d,DIST)
  type <- DIST[type]
  if(!is.null(cl)){
    s <- parLapply(cl=cl,X=H,fun=function(x,g,type, ...){
      res <- netdist(g,x,d=type, n.cores=1, ...)[[type]]
      return(res)
    }, g=g, type=type, ...)
  
  }else{
    s <- lapply(X=H,FUN=function(x,g,type, ...){
      res <- netdist(g,x,d=type, n.cores=1, ...)[[type]]
      return(res)
    },g=g,type=type, ...)
  }
  return(unlist(s))
}

## Stability indicator SI
netsiSI <- function(H, d, ...){
  DIST <- c("HIM","IM","H")
  type <- pmatch(d,DIST)
  type <- DIST[type]
  s <- netdist(H,d=type, ...)[[1]]
  return(s[upper.tri(s)])
}

## Degree stability
netsiSd <- function(H,cl){
  if (length(H)){
    n <- ncol(H[[1]])
  } else {
    stop("No adjacency matrix computed",call.=FALSE)
  }
  
  if (!is.null(cl)){
    ## Parallel computation
    dd <- parLapply(cl=cl, X=H, rowSums)
  } else {
    ## One core computation
    dd <- lapply(H, rowSums)
  }
  dd <- matrix(unlist(dd),ncol=n, byrow=TRUE)
  colnames(dd) <- colnames(H[[1]])
  return(dd)
}

## Edge stability
netsiSw <- function(H,cl){
  if (length(H))
    n <- nrow(H[[1]]) else stop("List of adjacency matrices do not exist")
  
  com <- combn(1:n, 2)
  
  if (!is.null(cl)){
    ## Parallel computation
    tmp <- parLapply(cl, H,
                     function(x,com){
                       sapply(1:ncol(com),
                              function(y, x, allcom){
                                x[allcom[1,y],allcom[2,y]]},
                              x=x, allcom=com)
                     },
                     com=com)
  } else {
    ## One core computation
    tmp <- lapply(H,
                  function(x,com){
                    sapply(1:ncol(com),
                           function(y, x, allcom){
                             x[allcom[1,y],allcom[2,y]]},
                           x=x, allcom=com)
                  },
                  com=com)
  }
  
  ## Set up the results in a matrix
  m <- matrix(unlist(tmp), ncol=length(H))
  rownames(m) <- paste(com[1,],com[2,],sep="-")
  colnames(m) <- names(H)
  
  return(m)
}

## Function for the computation of resampling indexes
resamplingIDX <- function(N,method="montecarlo", k=3, h=20){
  
  ## Check of resampling methods
  METHODS <- c('montecarlo','LOO','kCV')
  method <- pmatch(method, METHODS)
  
  if(is.na(method))
    stop("invalid resampling method")
  if(method == -1)
    stop("ambiguous resampling method")
  
  ## Montecarlo
  if (method==1L)
    take <- lapply(1:h,function(x){tmp <- sample(seq(N),floor(N*(1-(1/k)))) 
                                   return(tmp)})
  
  ## Leave One Out
  if (method==2L){
    if (k!=1)
      warning("h is different than 1 but method is set to LOO (Leave One Out cross-validation schema).\nh will be ignored.")
    h <- N
    take <- lapply(1:N,function(x,allid){return(allid[which(allid!=x)])},allid=1:N)
  }
  
  ## K-fold cross-validation
  if (method==3L){
    if (k>=N){
      stop("Number of fold bigger than samples in the dataset!!!")
    }else{
      take <- list()
      for (H in seq(h)){
        tmpcv <- sample(seq(N),N)
        s <- split(tmpcv,cut(seq(N),k))
        names(s) <- NULL
        take <- c(take,lapply(s,function(x,idx){tmp <- idx[!is.element(idx,x)]
                                                return(tmp)},idx=tmpcv))
      }
    }
  }
  
  ## return a list with indexes
  return(take)
}

compute.indicator <- function(x){
  ##  Compute the indicator value as Range/Mean over resamplings
  tmp <- range(x)
  mm <- mean(x, na.rm=TRUE)
  if (all.equal(mm,0)==TRUE){
    rr <- NA
  } else {
    rr <- (tmp[2] - tmp[1])/mm
  }
  return(rr)
}
