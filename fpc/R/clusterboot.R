
clujaccard <- function(c1,c2,zerobyzero=NA){
#  print(c1)
#  print(c2)
  if (sum(c1)+sum(c2)-sum(c1 & c2)==0) out <- zerobyzero
  else
    out <- sum(c1 & c2)/(sum(c1)+sum(c2)-sum(c1 & c2))
  out
}


noisemclustCBI <- function(data, G=NULL, k=NULL, emModelNames=NULL,
                           nnk=0, hcmodel = NULL,
                         Vinv = NULL, summary.out=FALSE, ...){
#  require(mclust)
#  require(prabclus)
  if (!is.null(k)) G <- k
  data <- as.matrix(data)
#  print(str(data))
  if (nnk > 0) {
    noise <- as.logical(1 - NNclean(data, nnk)$z)
# print(noise)
    if (!is.null(hcmodel)) 
      hcPairs <- hc(modelName = hcmodel, data = data)
    if (is.null(Vinv) & is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),...)
    if (!is.null(Vinv) & is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),
                      Vinv=Vinv,...)
    if (is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),...)
    if (!is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),
                      Vinv=Vinv,...)
  }
  else {
        if (!is.null(hcmodel)) {
            hcPairs <- hc(modelName = hcmodel, data = data)
            c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs),
                      Vinv=Vinv,...)           
        }
        else
            c1 <- mclustBIC(data, G, emModelNames,
                      Vinv=Vinv,...)           
        noise <- rep(0, nrow(data))
  }
# print(c1)
  sc1 <- summary(c1,data)
# print(str(sc1))  
  sc1c <- sc1$classification
  cl <- list()
  if (sc1$G==0)
    sc1c <- rep(0,nrow(data))
  nc <- nccl <- max(sc1c)
  if (sum(sc1c==0)>0){
    nc <- nccl+1
    sc1c[sc1c==0] <- nc
  }
  for (i in 1:nc)
      cl[[i]] <- sc1c == i
  if (summary.out)
    out <- list(result = c1, nc = nc, nccl = nccl, clusterlist = cl, 
        partition = sc1c, nnk = nnk, initnoise = as.logical(noise),
                mclustsummary=sc1,
        clustermethod = "mclustBIC")
  else
    out <- list(result = c1, nc = nc, nccl = nccl, clusterlist = cl, 
        partition = sc1c, nnk = nnk, initnoise = as.logical(noise), 
        clustermethod = "mclustBIC")
  out
}

distnoisemclustCBI <- function(dmatrix, G=NULL, k=NULL, emModelNames=NULL,
                           nnk=0, hcmodel = NULL,
                         Vinv = NULL, mdsmethod="classical",
                            mdsdim=4, summary.out=FALSE,
                               points.out=FALSE, ...){
  dmatrix <- as.matrix(dmatrix)
  if (!is.null(k)) G <- k
  n <- ncol(dmatrix)
#  require(MASS)
#  require(prabclus)
#  require(mclust)
  if (mdsmethod != "classical") {
    mindm <- min(dmatrix[dmatrix > 0])/10
    dmatrix[dmatrix<mindm] <- mindm
  }
  data <- switch(mdsmethod, classical = cmdscale(dmatrix, k = mdsdim), 
        kruskal = isoMDS(dmatrix, k = mdsdim)$points, sammon =
                sammon(dmatrix, k = mdsdim)$points)
  data <- as.matrix(data) 
  if (nnk > 0) {
    noise <- as.logical(1 - NNclean(data, nnk)$z)
# print(noise)
    if (!is.null(hcmodel)) 
      hcPairs <- hc(modelName = hcmodel, data = data)
    if (is.null(Vinv) & is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),...)
    if (!is.null(Vinv) & is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),
                      Vinv=Vinv,...)
    if (is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),...)
    if (!is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),
                      Vinv=Vinv,...)
  }
  else {
        if (!is.null(hcmodel)) {
            hcPairs <- hc(modelName = hcmodel, data = data)
            c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs),
                      Vinv=Vinv,...)           
        }
        else
            c1 <- mclustBIC(data, G, emModelNames,
                      Vinv=Vinv,...)           
        noise <- rep(0, nrow(data))
  }
# print(c1)
  sc1 <- summary(c1,data)
# print(str(sc1))  
  sc1c <- sc1$classification
  cl <- list()
  if (sc1$G==0)
    sc1c <- rep(0,nrow(data))
  nc <- nccl <- max(sc1c)
  if (sum(sc1c==0)>0){
    nc <- nccl+1
    sc1c[sc1c==0] <- nc
  }
  for (i in 1:nc)
      cl[[i]] <- sc1c == i
  out <- list(result = c1, nc = nc, nccl = nccl, clusterlist = cl, 
        partition = sc1c, nnk = nnk, initnoise = as.logical(noise), 
        clustermethod = "mclustBIC")
  if (summary.out)
    out$mclustsummary <- sc1
  if (points.out)
    out$points <- data
  out
}

mergenormCBI <- function(data, G=NULL, k=NULL, emModelNames=NULL, nnk=0,
                         hcmodel = NULL,
                         Vinv = NULL, mergemethod="bhat",
                         cutoff=0.1,...){
#  require(mclust)
#  require(prabclus)
  if (!is.null(k)) G <- k  
  if (nnk > 0) {
    noise <- as.logical(1 - NNclean(data, nnk)$z)
# print(noise)    
    if (!is.null(hcmodel)) 
      hcPairs <- hc(modelName = hcmodel, data = data)
    if (is.null(Vinv) & is.null(hcmodel)){ 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),...)
#      print("mclust done")
    }
    if (!is.null(Vinv) & is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames, initialization=list(noise=noise),
                      Vinv=Vinv,...)
    if (is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),...)
    if (!is.null(Vinv) & !is.null(hcmodel)) 
      c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs, noise=noise),
                      Vinv=Vinv,...)
  }
  else {
        if (!is.null(hcmodel)) {
            hcPairs <- hc(modelName = hcmodel, data = data)
            c1 <- mclustBIC(data, G, emModelNames,
                      initialization=list(hcPairs=hcPairs),
                      Vinv=Vinv,...)           
        }
        else
            c1 <- mclustBIC(data, G, emModelNames,
                      Vinv=Vinv,...)           
        noise <- rep(0, nrow(data))
  }
  sc1 <- summary(c1,data)
#  print("sum")
#  print(max(sc1$classification))
  jsc1 <- mergenormals(data,sc1,method=mergemethod,cutoff=cutoff,...)
#  print(jsc1$clusternumbers)
#  print("merge")
  renumcl <- jsc1$clustering
#  cln <- 1
#  for (i in 1:max(jsc1$clustering)){
#     if(sum(jsc1$clustering==i)>0){
#       renumcl[jsc1$clustering==i] <- cln
#       cln <- cln+1
#     }
#   }
  cl <- list()
  nc <- nccl <- max(jsc1$clustering)
  if (sum(jsc1$clustering==0)>0)
    nc <- nc+1
  for (i in 1:nccl)
    cl[[i]] <- renumcl == i
  if (nc>nccl)
    cl[[nc]] <- renumcl == 0
  out <- list(result = jsc1, nc = nc, nccl = nccl, clusterlist = cl, 
        partition = renumcl, nnk = nnk, initnoise = as.logical(noise), 
        clustermethod = "mclust/mergenormals")
  out
}


hclustCBI <- function(data,k,cut="number",method,scaling=TRUE,noisecut=0,...){
  if(!identical(scaling,FALSE))
    sdata <- scale(data,center=TRUE,scale=scaling)
  else
    sdata <- data
  n <- nrow(data)
  noise <- FALSE
  c1 <- hclust(dist(sdata),method=method)
  if (cut=="number")
    partition <- cutree(c1,k=k)
  else
    partition <- cutree(c1,h=k)
  cl <- list()
  nc <- max(partition)
  clsizes <- numeric(0)
  for (i in 1:nc) clsizes[i] <- sum(partition==i)
  ncn <- sum(clsizes>noisecut)
  if (ncn<nc){
    noise <- TRUE
    newcln <- (1:nc)[clsizes>noisecut]
    nc <- ncn+1
    newpart <- rep(nc,n)
    for (i in 1:ncn)
      newpart[partition==newcln[i]] <- i
    partition <- newpart
  }
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,noise=noise,
              nc=nc,clusterlist=cl,partition=partition,
              clustermethod="hclust/cutree")
  out
}

hclusttreeCBI <- function(data,minlevel=2,method,scaling=TRUE,...){
  if(!identical(scaling,FALSE))
    sdata <- scale(data,center=TRUE,scale=scaling)
  else
    sdata <- data
  c1 <- hclust(dist(sdata),method=method)
  n <- nrow(data)
  clist <- list()
  for (i in 1:n){
    clist[[i]] <- rep(FALSE,n)
    clist[[i]][i] <- TRUE
  }
  clcount <- n
  for (j in 1:(n-2)){
    clcount <- clcount+1
    if (c1$merge[j,1]<0) clist1 <- clist[[-c1$merge[j,1]]]
    else clist1 <- clist[[n+c1$merge[j,1]]]
    if (c1$merge[j,2]<0) clist2 <- clist[[-c1$merge[j,2]]]
    else clist2 <- clist[[n+c1$merge[j,2]]]
    clist[[clcount]] <- clist2 | clist1
  }
  clusterlist <- list()
  if (minlevel==1){
    clusterlist <- clist
    nc <- clcount
  }
  else{
    for (j in (n+minlevel-1):clcount)
      clusterlist[[j-minlevel-n+2]] <- clist[[j]]
    nc <- clcount-minlevel-n+2
  }
#  print(nc)
#  print(sc1)
  out <- list(result=c1,nc=nc,clusterlist=clusterlist,partition=cutree(c1,2),
              clustermethod="hclust, full tree")
  out
}

# for clara & pam, data called x
claraCBI <- function(data,k,usepam=TRUE,diss=inherits(data,"dist"),...){
  if (usepam)    
    c1 <- pam(data,k=k,diss=diss,...)
  else
    c1 <- clara(data,k=k,...)
  partition <- c1$clustering
  cl <- list()
  nc <- k
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              clustermethod="clara/pam")
  out
}

speccCBI <- function(data,k,...){
#  require(kernlab)
  data <- as.matrix(data)
  options(show.error.messages = FALSE)
  c1 <- try(specc(data,centers=k,...))
  options(show.error.messages = TRUE)
  if (class(c1)=="try-error"){
    partition <- rep(1,nrow(data))
    cat("Function specc returned an error, probably a one-point cluster.\n All observations were classified to cluster 1.\n")
  }
  else
    partition <- c1@.Data
  cl <- list()
  nc <- k
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              clustermethod="spectral")
  out
}

# tclustCBI <- function(data,k,trim=0.05,...){
#   if(require(tclust)){
#     data <- as.matrix(data)
#     c1 <- tclust(data,k=k,alpha=trim,...)
#     sc1c <- c1$cluster
#     cl <- list()
#     nc <- nccl <- max(sc1c)
#     if (sum(sc1c==0)>0){
#       nc <- nccl+1
#       sc1c[sc1c==0] <- nc
#     }
#     for (i in 1:nc)
#       cl[[i]] <- sc1c == i
#     out <- list(result=c1,nc=nc,nccl=nccl,clusterlist=cl,partition=sc1c,
#               clustermethod="tclust")
#     out
#   }
#   else
#     warning("tclust could not be loaded")    
# }

trimkmeansCBI <- function(data,k,scaling=TRUE,trim=0.1,...){
    c1 <- trimkmeans(data,k=k,scaling=scaling,trim=trim,...)
    partition <- c1$classification
    cl <- list()
    nc <- k+1
    nccl <- k
#  print(nc)
#  print(sc1)
    for (i in 1:nc)
      cl[[i]] <- partition==i
    out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              nccl=nccl,clustermethod="trimkmeans")
    out      
}

kmeansCBI <- function(data,krange,k=NULL,scaling=FALSE,runs=1,criterion="ch",...){
  if (!is.null(k)) krange <- k
  if(!identical(scaling,FALSE))
    sdata <- scale(data,center=TRUE,scale=scaling)
  else
    sdata <- data
  c1 <- kmeansruns(sdata,krange,runs=runs,criterion=criterion,...)
  partition <- c1$cluster
  cl <- list()
  nc <- krange
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              clustermethod="kmeans")
  out
}

pamkCBI <- function (data, krange = 2:10,k=NULL,
                     criterion="asw", usepam=TRUE,
                     scaling = TRUE, diss = inherits(data,"dist"), ...) 
{
#    require(cluster)
    if (!is.null(k)) krange <- k
    c1 <- pamk(data, krange = krange, criterion=criterion, usepam=usepam,
               scaling = scaling, diss = diss, 
        ...)
    partition <- c1$pamobject$clustering
    cl <- list()
    nc <- c1$nc
#    print(nc)
    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition, 
        clustermethod = "pam/estimated k",criterion=criterion)
    out
}


# for dbscan, data called x
dbscanCBI <- function(data,eps,MinPts,diss=inherits(data,"dist"),...){
  if (diss)
    c1 <- dbscan(data,eps,MinPts,method="dist",...)
  else
    c1 <- dbscan(data,eps,MinPts,...)
#  plot(c1, data)
  partition <- c1$cluster
  cl <- list()
  nccl <- max(partition)
  partition[partition==0] <- nccl+1
  nc <- max(partition)
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,nc=nc,nccl=nccl,clusterlist=cl,partition=partition,
              clustermethod="dbscan")
  out
}

mahalCBI <- function(data,clustercut=0.5,...){
  c1 <- fixmahal(data,...)
  partition <- rep(1,nrow(data))
  cl <- fpclusters(c1)
  nc <- length(cl)
  if (nc>0){
    for (i in 1:nc)
      cl[[i]] <- as.integer(cl[[i]]>clustercut)
  }
  else{
    nc <- 1
    cl[[1]] <- as.integer(c1$g[[1]]>clustercut)
  }
  out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              clustermethod="fixmahal")
  out
}


clusterboot <- function(data,B=100,
                        distances=(class(data)=="dist"),
                        bootmethod="boot",
                        bscompare=TRUE, multipleboot=FALSE,
                        jittertuning=0.05, noisetuning=c(0.05,4),
                        subtuning=floor(nrow(data)/2),
                        clustermethod,noisemethod=FALSE,
                        count=TRUE,
                        showplots=FALSE,dissolution=0.5,
                        recover=0.75,seed=NULL,...){
  sumlogic <- function(x,y,relation="eq")
    switch (relation,
          eq = sum(x==y, na.rm=TRUE),
          s = sum(x<y, na.rm=TRUE),
          l = sum(x>y, na.rm=TRUE),
          se = sum(x<=y, na.rm=TRUE),
          le = sum(x>=y, na.rm=TRUE))

  if (!is.null(seed)) set.seed(seed)
  invisible(distances)
  data <- as.matrix(data)
  if (distances & showplots) dpoints <- cmdscale(data)
  n <- nrow(data)
  p <- ncol(data)
  cod <- cov(data)
  md <- colMeans(data)
  lb <- length(bootmethod)
  c1 <- clustermethod(data,...)
#  print(noisemethod)
#  print(str(c1))
  if (noisemethod){
    if (c1$nccl==0)
      stop("No clusters, only noise estimated!")
  }
  else
    c1$nccl <- c1$nc
  bootresult <- jitterresult <- noiseresult <-
    bojitresult <- subsetresult <- matrix(0,nrow=c1$nc,ncol=B)
  if (("jitter" %in% bootmethod) | ("bojit" %in% bootmethod)){
    jsd <- numeric(0)
    ecd <- eigen(cod, symmetric=TRUE)
    ecd$values[ecd$values<0] <- 0
    ecd$values[is.na(ecd$values)] <- 0
    rotdata <- data %*% solve(t(ecd$vectors))
    for (i in 1:p){
      sx <- sort(rotdata[,i])
      dx <- sx[2:n]-sx[1:(n-1)]
      dx <- dx[dx>0]
      jsd[i] <- quantile(dx,jittertuning)
    }
  }
  if ("noise" %in% bootmethod){
    ecd <- eigen(cod, symmetric=TRUE)
    ecd$values[ecd$values<0] <- 0
  }
  if (showplots){
    if (distances)
      plot(dpoints,pch=sapply(c1$partition,toString),col=c1$partition)
    else 
      plot(data,pch=sapply(c1$partition,toString),col=c1$partition)
  }
#  matrixlist <- list()
  for (l in 1:lb){
    for (i in 1:B){
      if (count) cat(bootmethod[l],i,"\n")
      if (bootmethod[l]=="boot"){
        bsamp <- sample(n,n,replace=TRUE)
        if (!multipleboot) bsamp <- unique(bsamp)
        if (distances)
          mdata <- data[bsamp,bsamp]
        else
          mdata <- data[bsamp,]
      }
      if (bootmethod[l]=="subset"){
        bsamp <- sample(n,subtuning,replace=FALSE)
        if (distances)
          mdata <- data[bsamp,bsamp]
        else
          mdata <- data[bsamp,]
      }
      if (bootmethod[l]=="jitter"){
        jnoise <- matrix(0,ncol=p,nrow=n)
        for (j in 1:p)
          jnoise[,j] <- rnorm(n,sd=jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data+jnoise
        bsamp <- 1:n
      }
      if (bootmethod[l]=="bojit"){
        bsamp <- sample(n,n,replace=TRUE)
        jnoise <- matrix(0,ncol=p,nrow=n)
        for (j in 1:p)
          jnoise[,j] <- rnorm(n,sd=jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data[bsamp,]+jnoise
      }
      if (bootmethod[l]=="noise"){
        noiseind <- as.logical(rbinom(n,1,noisetuning[1]))
        nn <- sum(noiseind)
        jnoise <- matrix(0,ncol=p,nrow=nn)
        for (j in 1:p)
          jnoise[,j] <- runif(nn,min=-noisetuning[2]*sqrt(ecd$values[j]),
                              max=noisetuning[2]*sqrt(ecd$values[j]))
        jnoise <- t(t(jnoise %*% t(ecd$vectors))+md) 
        mdata <- data
        mdata[noiseind,] <- jnoise
        bsamp <- (1:n)[!noiseind]
      }
# print(mdata)
      bc1 <- clustermethod(mdata,...)
# print("clustermethod done")
      if (showplots){
        if (distances)
          plot(dpoints[bsamp,],pch=sapply(bc1$partition,toString),
               col=bc1$partition)
        else
          plot(mdata,pch=sapply(bc1$partition,toString),col=bc1$partition)
      }
      if (noisemethod){
        effnc1 <- c1$nccl
        effnb1 <- bc1$nccl
      }
      else{
        effnc1 <- c1$nc
        effnb1 <- bc1$nc
      }
      for (j in 1:effnc1){
        maxgamma <- 0
        if (effnb1>0){
          for (k in 1:effnb1){
            if (multipleboot){
              if (bscompare) ncases <- 1:n
              else{
                ncases <- 1
#                print(bsamp)
#                print(str(bsamp))
                m <- 2 # for (m in 2:n)
                if (m<=n){
                  if (!(bsamp[m] %in% bsamp[1:(m-1)]))
                    ncases <- c(ncases,m)
                  m <- m+1
                }
              }
            }
            else ncases <- 1:length(bsamp)
#          print(ncases)
#            if (j==c1$nc){
#          print(c1$clusterlist[[j]][bsamp][ncases])
#          print(bc1$clusterlist[[k]][ncases])
#        }
#            print(bc1$nc)
#            print(bc1$nccl)
#            print(bc1$partition)
#            print(str(bc1$clusterlist))
            cg <- switch(bootmethod[l],
              boot=clujaccard(c1$clusterlist[[j]][bsamp][ncases],
                bc1$clusterlist[[k]][ncases],zerobyzero=0),
              bojit=clujaccard(c1$clusterlist[[j]][bsamp][ncases],
                bc1$clusterlist[[k]][ncases],zerobyzero=0),
              subset=clujaccard(c1$clusterlist[[j]][bsamp],
                bc1$clusterlist[[k]],zerobyzero=0),
              jitter=clujaccard(c1$clusterlist[[j]],bc1$clusterlist[[k]],zerobyzero=0),
              noise=clujaccard(c1$clusterlist[[j]][!noiseind],
                 bc1$clusterlist[[k]][!noiseind],zerobyzero=0))
            if (cg>maxgamma) maxgamma <- cg
          }
        }
        if (bootmethod[l]=="boot") bootresult[j,i] <- maxgamma
        if (bootmethod[l]=="subset") subsetresult[j,i] <- maxgamma
        if (bootmethod[l]=="bojit") bojitresult[j,i] <- maxgamma
        if (bootmethod[l]=="jitter") jitterresult[j,i] <- maxgamma
        if (bootmethod[l]=="noise") noiseresult[j,i] <- maxgamma
      }
      if (noisemethod){
        if (c1$nc>c1$nccl){
          j <- c1$nc
          if (bc1$nc>bc1$nccl)
            maxgamma <- switch(bootmethod[l],
              boot=clujaccard(c1$clusterlist[[c1$nc]][bsamp][ncases],
                bc1$clusterlist[[bc1$nc]][ncases],zerobyzero=0),
              bojit=clujaccard(c1$clusterlist[[c1$nc]][bsamp][ncases],
                bc1$clusterlist[[bc1$nc]][ncases],zerobyzero=0),
              subset=clujaccard(c1$clusterlist[[c1$nc]][bsamp],
                bc1$clusterlist[[bc1$nc]],zerobyzero=0),
              jitter=clujaccard(c1$clusterlist[[c1$nc]],bc1$clusterlist[[bc1$nc]],zerobyzero=0),
              noise=clujaccard(c1$clusterlist[[c1$nc]][!noiseind],
                 bc1$clusterlist[[bc1$nc]][!noiseind],zerobyzero=0))
          else
            maxgamma <- 0
          if (bootmethod[l]=="boot") bootresult[j,i] <- maxgamma
          if (bootmethod[l]=="subset") subsetresult[j,i] <- maxgamma
          if (bootmethod[l]=="bojit") bojitresult[j,i] <- maxgamma
          if (bootmethod[l]=="jitter") jitterresult[j,i] <- maxgamma
          if (bootmethod[l]=="noise") noiseresult[j,i] <- maxgamma
        }
      }      
    }
  }

  if (!("boot" %in% bootmethod))
    bootresult <- bootmean <- bootbrd <- bootrecover <- NULL
  else{
              bootmean=apply(bootresult,1,mean, na.rm=TRUE)
              bootbrd=apply(bootresult,1,sumlogic,y=dissolution,relation="se")
              bootrecover=apply(bootresult,1,sumlogic,y=recover,relation="l")
       }
  if (!("jitter" %in% bootmethod))
    jitterresult <- jittermean <- jitterbrd <- jitterrecover <- NULL
  else{
    jittermean=apply(jitterresult,1,mean, na.rm=TRUE)
    jitterbrd=apply(jitterresult,1,sumlogic,y=dissolution,relation="se")
    jitterrecover=apply(jitterresult,1,sumlogic,y=recover,relation="l")
  }
  if (!("subset" %in% bootmethod)) subsetresult <- subsetmean <-
    subsetbrd <- subsetrecover <- NULL
  else{
    subsetmean=apply(subsetresult,1,mean, na.rm=TRUE)
    subsetbrd=apply(subsetresult,1,sumlogic,y=dissolution,relation="se")
    subsetrecover=apply(subsetresult,1,sumlogic,y=recover,relation="l")
  }
  if (!("noise" %in% bootmethod)) noiseresult <- noisemean <-
    noisebrd <- noiserecover <- NULL
  else{
    noisemean=apply(noiseresult,1,mean, na.rm=TRUE)
    noisebrd=apply(noiseresult,1,sumlogic,y=dissolution,relation="se")
    noiserecover=apply(noiseresult,1,sumlogic,y=recover,relation="l")
  }
  if (!("bojit" %in% bootmethod)) bojitresult <- bojitmean <-
    bojitbrd <- bojitrecover <- NULL
  else{
              bojitmean=apply(bojitresult,1,mean, na.rm=TRUE)
              bojitbrd=apply(bojitresult,1,sumlogic,y=dissolution,relation="se")
              bojitrecover=apply(bojitresult,1,sumlogic,y=recover,
                relation="l")
  }
  if (showplots){
    if (distances)
      plot(dpoints,pch=sapply(c1$partition,toString),col=c1$partition)
    else 
      plot(data,pch=sapply(c1$partition,toString),col=c1$partition)
  }
  out <- list(result=c1,partition=c1$partition,
              nc=c1$nc, nccl=c1$nccl,
              clustermethod=c1$clustermethod,
              B=B,
              noisemethod=noisemethod,
              bootmethod=bootmethod,
              multipleboot=multipleboot,
              dissolution=dissolution,
              recover=recover,
              bootresult=bootresult,
              bootmean=bootmean,
              bootbrd=bootbrd,
              bootrecover=bootrecover,
              jitterresult=jitterresult,
              jittermean=jittermean,
              jitterbrd=jitterbrd,
              jitterrecover=jitterrecover,
              subsetresult=subsetresult,
              subsetmean=subsetmean,
              subsetbrd=subsetbrd,
              subsetrecover=subsetrecover,
              bojitresult=bojitresult,
              bojitmean=bojitmean,
              bojitbrd=bojitbrd,
              bojitrecover=bojitrecover,
              noiseresult=noiseresult,
              noisemean=noisemean, 
              noisebrd=noisebrd,
              noiserecover=noiserecover)
  class(out) <- "clboot"
  out
}

print.clboot <- function(x, statistics=c("mean","dissolution","recovery"),...){
  cat("* Cluster stability assessment *\n")
  cat("Cluster method: ",x$clustermethod,"\n")
  cat("Full clustering results are given as parameter result\n")
  cat("of the clusterboot object, which also provides further statistics\n")
  cat("of the resampling results.\n")
  cat("Number of resampling runs: ",x$B,"\n\n")
  cat("Number of clusters found in data: ",x$nc,"\n\n")
  if (x$noisemethod & x$nc>x$nccl)
   cat("The last cluster corresponds to noise/outliers\n\n")   
  for (i in 1:length(x$bootmethod)){
    if (x$bootmethod[i]=="boot"){
      cat(" Clusterwise Jaccard bootstrap ")
      if (!x$multipleboot) cat("(omitting multiple points) ")
      if ("mean" %in% statistics){
        cat("mean:\n")
        print(x$bootmean)
      }
      if ("dissolution" %in% statistics){
        cat("dissolved:\n")
        print(x$bootbrd)
      }
      if ("recovery" %in% statistics){
        cat("recovered:\n")
        print(x$bootrecover)
      }
    }
    if (x$bootmethod[i]=="subset"){
      cat(" Clusterwise Jaccard subsetting ")
      if ("mean" %in% statistics){
        cat("mean:\n")
        print(x$subsetmean)
      }
      if ("dissolution" %in% statistics){
        cat("dissolved:\n")
        print(x$subsetbrd)
      }
      if ("recovery" %in% statistics){
        cat("recovered:\n")
        print(x$subsetrecover)
      }
    }
    if (x$bootmethod[i]=="jitter"){
      cat(" Clusterwise Jaccard jittering ")
      if ("mean" %in% statistics){
        cat("mean:\n")
        print(x$jittermean)
      }
      if ("dissolution" %in% statistics){
        cat("dissolved:\n")
        print(x$jitterbrd)
      }
      if ("recovery" %in% statistics){
        cat("recovered:\n")
        print(x$jitterrecover)
      }
    }
    if (x$bootmethod[i]=="noise"){
      cat(" Clusterwise Jaccard replacement by noise ")
      if ("mean" %in% statistics){
        cat("mean:\n")
        print(x$noisemean)
      }
      if ("dissolution" %in% statistics){
        cat("dissolved:\n")
        print(x$noisebrd)
      }
      if ("recovery" %in% statistics){
        cat("recovered:\n")
        print(x$noiserecover)
      }
    }
    if (x$bootmethod[i]=="bojit"){
      cat(" Clusterwise Jaccard bootstrap/jittering ")
      if ("mean" %in% statistics){
        cat("mean:\n")
        print(x$bojitmean)
      }
      if ("dissolution" %in% statistics){
        cat("dissolved:\n")
        print(x$bojitbrd)
      }
      if ("recovery" %in% statistics){
        cat("recovered:\n")
        print(x$bojitrecover)
      }
    }
  }
  invisible(x)
}

plot.clboot <- function(x,xlim=c(0,1),breaks=seq(0,1,by=0.05),...){
  par(mfrow=c(x$nc,length(x$bootmethod)))
  for (j in 1:x$nc)
    for (i in 1:length(x$bootmethod)){
      if (x$bootmethod[i]=="boot")
        hist(x$bootresult[j,],xlim=xlim,breaks=breaks,
             xlab="Jaccard similarity",
             main=paste("Bootstrap, cluster ",j))
      if (x$bootmethod[i]=="subset")
        hist(x$subsetresult[j,],xlim=xlim,breaks=breaks,
             xlab="Jaccard similarity",
             main=paste("Subsetting, cluster ",j))
      if (x$bootmethod[i]=="jitter")
        hist(x$jitterresult[j,],xlim=xlim,breaks=breaks,
             xlab="Jaccard similarity",
             main=paste("Jittering, cluster ",j))
      if (x$bootmethod[i]=="noise")
        hist(x$noiseresult[j,],xlim=xlim,breaks=breaks,
             xlab="Jaccard similarity",
             main=paste("Replacement by noise, cluster ",j))
      if (x$bootmethod[i]=="bojit")
        hist(x$bojitresult[j,],xlim=xlim,breaks=breaks,
             xlab="Jaccard similarity",
             main=paste("Bootstrap/jittering, cluster ",j))
    }
  par(mfrow=c(1,1))
  invisible()
}

  


disttrimkmeansCBI <- function(dmatrix,k,scaling=TRUE,trim=0.1,
                           mdsmethod="classical",
                            mdsdim=4,...){
    dmatrix <- as.matrix(dmatrix)
    n <- ncol(dmatrix)
  #  require(MASS)
    if (mdsmethod != "classical") {
      mindm <- min(dmatrix[dmatrix > 0])/10
      dmatrix[dmatrix<mindm] <- mindm 
    }
    data <- switch(mdsmethod, classical = cmdscale(dmatrix, k = mdsdim), 
          kruskal = isoMDS(dmatrix, k = mdsdim)$points, sammon =
                  sammon(dmatrix, k = mdsdim)$points)
    c1 <- trimkmeans(data,k=k,scaling=scaling,trim=trim,...)
    partition <- c1$classification
    cl <- list()
    nccl <- k
    nc <- k+1
  #  print(nc)
  #  print(sc1)
    for (i in 1:nc)
      cl[[i]] <- partition==i
    out <- list(result=c1,nc=nc,nccl=nccl,clusterlist=cl,partition=partition,
                clustermethod="trimkmeans plus MDS")
    out
}


disthclustCBI <- function(dmatrix,k,cut="number",method,noisecut=0,...){
  n <- nrow(as.matrix(dmatrix))
  c1 <- hclust(as.dist(dmatrix),method=method)
  noise <- FALSE
  if (cut=="number")
    partition <- cutree(c1,k=k)
  else
    partition <- cutree(c1,h=k)
  nc <- max(partition)
  clsizes <- numeric(0)
  for (i in 1:nc) clsizes[i] <- sum(partition==i)
  ncn <- sum(clsizes>noisecut)
  if (ncn<nc){
    noise <- TRUE
    newcln <- (1:nc)[clsizes>noisecut]
    nc <- ncn+1
    newpart <- rep(nc,n)
    for (i in 1:ncn)
      newpart[partition==newcln[i]] <- i
    partition <- newpart
  }
  cl <- list()
#  print(nc)
#  print(sc1)
  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,noise=noise,nc=nc,nccl=ncn,
              clusterlist=cl,partition=partition,
              clustermethod="hclust")
  out
}

disthclusttreeCBI <- function(dmatrix,minlevel=2,method,...){
  n <- nrow(as.matrix(dmatrix))
  c1 <- hclust(as.dist(dmatrix),method=method)
  clist <- list()
  for (i in 1:n){
    clist[[i]] <- rep(FALSE,n)
    clist[[i]][i] <- TRUE
  }
  clcount <- n
  for (j in 1:(n-2)){
    clcount <- clcount+1
    if (c1$merge[j,1]<0) clist1 <- clist[[-c1$merge[j,1]]]
    else clist1 <- clist[[n+c1$merge[j,1]]]
    if (c1$merge[j,2]<0) clist2 <- clist[[-c1$merge[j,2]]]
    else clist2 <- clist[[n+c1$merge[j,2]]]
    clist[[clcount]] <- clist2 | clist1
  }
  clusterlist <- list()
  if (minlevel==1){
    clusterlist <- clist
    nc <- clcount
  }
  else{
    for (j in (n+minlevel-1):clcount)
      clusterlist[[j-minlevel-n+2]] <- clist[[j]]
    nc <- clcount-minlevel-n+2
  }
#  print(nc)
#  print(sc1)
  out <- list(result=c1,nc=nc,clusterlist=clusterlist,partition=cutree(c1,2),
              clustermethod="hclust, full tree")
  out
}


# clustering<0: to be predicted
classifnp <- function(data,clustering,
                      method="centroid",cdist=NULL,
                      centroids=NULL,nnk=1){
  data <- as.matrix(data)
  k <- max(clustering)
  p <- ncol(data)
  n <- nrow(data)
  topredict <- clustering<0
  if(method=="averagedist"){
    if(is.null(cdist))
      cdist <- as.matrix(dist(data))
    else
      cdist <- as.matrix(cdist)
    prmatrix <- matrix(0,ncol=k,nrow=sum(topredict))
#    browser()
    for (j in 1:k)
      prmatrix[,j] <- rowMeans(as.matrix(cdist[topredict,clustering==j]))
    clpred <- apply(prmatrix,1,which.min)
    clustering[topredict] <- clpred
  }
  if(method=="centroid"){
#    require(class)
    if(is.null(centroids)){
      centroids <- matrix(0,ncol=p,nrow=k)
      for (j in 1:k)
        centroids[j,] <- colMeans(as.matrix(data[clustering==j,]))
#      browser()
    }
    clustering[topredict] <- knn1(centroids,data[topredict,],1:k)
  }
  if(method=="qda"){
    qq <- try(qda(data[!topredict,],
              grouping=as.factor(clustering[!topredict])),silent=TRUE)
    if (identical(attr(qq,"class"),"try-error"))
      qq <- lda(data[!topredict,],
              grouping=as.factor(clustering[!topredict]))        
    clustering[topredict] <- as.integer(predict(qq,data[topredict,])$class)
  }
  if(method=="knn")
    clustering[topredict] <- as.integer(knn(data[!topredict,],
                                            data[topredict,],
                                            as.factor(clustering[!topredict]),
                                            k=nnk))
  clustering
}
    
# method centroids doesn't allow centroids=NULL
classifdist <- function(cdist,clustering,
                      method="averagedist",
                      centroids=NULL,nnk=1){
  k <- max(clustering)
  n <- nrow(data)
  cdist <- as.matrix(cdist)
  topredict <- clustering<0
  if(method=="averagedist"){
    prmatrix <- matrix(0,ncol=k,nrow=sum(topredict))
    for (j in 1:k)
      prmatrix[,j] <- rowMeans(as.matrix(cdist[topredict,clustering==j]))
    clpred <- apply(prmatrix,1,which.min)
    clustering[topredict] <- clpred
  }
  if(method=="centroid")
    clustering[topredict] <- apply(cdist[topredict,centroids,drop=FALSE],1,which.min)
  if(method=="knn"){
    cdist[topredict,topredict] <- max(cdist)+1
    if (nnk==1){
      bestobs <- apply(cdist[topredict,,drop=FALSE],1,which.min)
      clustering[topredict] <- clustering[bestobs]
    }
    else{
      for(i in (1:n)[topredict]){
        bestobs <- order(cdist[i,])[1:k]
        clasum <- numeric(0)
        for (j in 1:k)
          clasum[j] <- sum(clustering[bestobs]==j)
        clustering[i] <- which.max(clasum)
      }
    }
  }
  clustering
}
    
nselectboot <- function (data, B = 50, distances = inherits(data, "dist"), clustermethod = NULL, 
    classification = "averagedist", krange = 2:10, count = FALSE, 
    nnk = 1, ...) 
{
    dista <- distances
    data <- as.matrix(data)
    if (classification == "average") {
        if (dista) 
            dmat <- data
        else dmat <- as.matrix(dist(data))
    }
    stab <- matrix(0, nrow = B, ncol = max(krange))
    n <- nrow(data)
    for (k in krange) {
        if (count) 
            cat(k, " clusters\n")
        for (i in 1:B) {
            if (count) 
                print(i)
            d1 <- sample(n, n, replace = TRUE)
            d2 <- sample(n, n, replace = TRUE)
            if (dista) {
                dmat1 <- data[d1, d1]
                dmat2 <- data[d2, d2]
            }
            else {
                dmat1 <- data[d1, ]
                dmat2 <- data[d2, ]
            }
#            print("start clustermethod")
            clm1 <- clustermethod(dmat1, k = k, ...)
#            cl1 <- clm1$partition
#            print("done 1")
            clm2 <- clustermethod(dmat2, k = k, ...)
#            cl2 <- clm2$partition
#            print("done 2")
            cj1 <- cj2 <- rep(-1, n)
            cj1[d1] <- clm1$partition
            cj2[d2] <- clm2$partition
#            browser()
            if (dista) {
                cj1 <- classifdist(data, cj1, method = classification, 
                  centroids = clm1$result$medoids, nnk = nnk)
                cj2 <- classifdist(data, cj2, method = classification, 
                  centroids = clm2$result$medoids, nnk = nnk)
            }
            else {
                centroids <- NULL
                if (classification == "centroid") {
                  if (identical(clustermethod, kmeansCBI)){ 
                    centroids1 <- clm1$result$centers
                    centroids2 <- clm2$result$centers
                  }
                  if (identical(clustermethod, claraCBI)){ 
                    centroids1 <- clm1$result$medoids
                    centroids2 <- clm2$result$medoids
                  }
                }
#                print("classifnp")
                cj1 <- classifnp(data, cj1, method = classification, 
                  centroids = centroids1, nnk = nnk)
#                print("done 1")
                cj2 <- classifnp(data, cj2, method = classification, 
                  centroids = centroids2, nnk = nnk)
#                plot(data,col=cj1,pch=cj1,main=paste(k,"clusters a -",i))
#                plot(data,col=cj2,pch=cj2,main=paste(k,"clusters b -",i))
#                print("done 2")
#                browser()
            }
#            print(" for n loop")
#            j <- 1
            ctable <- table(cj1,cj2)
            nck1 <- rowSums(ctable)
            stab[i,k] <- sum(nck1^2-rowSums(ctable^2))
#            browser()
#            for(j in 1:k)
#             for(j in 1:n){
#               cj1e <- cj1==cj1[j]
#               cj1u <- !cj1e
#               cj2e <- cj2==cj2[j]
#               cj2u <- !cj2e
#               stab[i,k] <- stab[i,k]+sum(cj1e & cj2u)+sum(cj2e & cj1u)
#              j <- j+1 }
#            print("done")
        } # for i
    } # for k
    stab <- stab/n^2
    stabk <- rep(NA, max(krange))
#   browser()
    for (k in krange) stabk[k] <- mean(stab[, k])
    kopt <- which.min(stabk)
    out <- list(kopt = kopt, stabk = stabk, stab = stab)
}
 
        
# nselectboot <- function (data, B = 50, distances = inherits(data, "dist"), clustermethod = NULL, 
#     classification = "averagedist", krange = 2:10, count = FALSE, 
#     nnk = 1, ...) 
# {
#     dista <- distances
#     data <- as.matrix(data)
#     if (classification == "average") {
#         if (dista) 
#             dmat <- data
#         else dmat <- as.matrix(dist(data))
#     }
#     stab <- matrix(0, nrow = B, ncol = max(krange))
#     n <- nrow(data)
#     for (k in krange) {
#         if (count) 
#             cat(k, " clusters\n")
#         for (i in 1:B) {
#             if (count) 
#                 print(i)
#             d1 <- sample(n, n, replace = TRUE)
#             d2 <- sample(n, n, replace = TRUE)
#             if (dista) {
#                 dmat1 <- data[d1, d1]
#                 dmat2 <- data[d2, d2]
#             }
#             else {
#                 dmat1 <- data[d1, ]
#                 dmat2 <- data[d2, ]
#             }
#             clm1 <- clustermethod(dmat1, k = k, ...)
#             clm2 <- clustermethod(dmat2, k = k, ...)
#             cj1 <- cj2 <- rep(-1, n)
#             cj1[d1] <- clm1$partition
#             cj2[d2] <- clm2$partition
#             if (dista) {
#                 cj1 <- classifdist(data, cj1, method = classification, 
#                   centroids = clm1$result$medoids, nnk = nnk)
#                 cj2 <- classifdist(data, cj2, method = classification, 
#                   centroids = clm2$result$medoids, nnk = nnk)
#             }
#             else {
#                 centroids <- NULL
#                 if (classification == "centroid") {
#                   if (identical(clustermethod, kmeansCBI)){ 
#                     centroids1 <- clm1$result$centers
#                     centroids2 <- clm2$result$centers
#                   }
#                   if (identical(clustermethod, claraCBI)){ 
#                     centroids1 <- clm1$result$medoids
#                     centroids2 <- clm2$result$medoids
#                   }
#                 }
#                 cj1 <- classifnp(data, cj1, method = classification, 
#                   centroids = centroids1, nnk = nnk)
#                 cj2 <- classifnp(data, cj2, method = classification, 
#                   centroids = centroids2, nnk = nnk)
#             }
#             j <- 1
#             if (j<=n){
#               cj1e <- cj1==cj1[j]
#               cj1u <- !cj1e
#               cj2e <- cj2==cj2[j]
#               cj2u <- !cj2e
#               stab[i,k] <- stab[i,k]+sum(cj1e & cj2u)+sum(cj2e & cj1u)
#               j <- j+1
#             }
#         }
#     }
#     stab <- stab/n^2
#     stabk <- rep(NA, max(krange))
#     for (k in krange) stabk[k] <- mean(stab[, k])
#     kopt <- which.min(stabk)
#     out <- list(kopt = kopt, stabk = stabk, stab = stab)
# }
