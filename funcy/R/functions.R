#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

initClust <- function(data, k, init, seed, nrep){
    nc <- dim(data)[1]
    if(init=="kmeans"){
        res <- suppressMessages(stepFlexclust(x=data, k=k, nrep=nrep, seed=seed, FUN=cclust))
        dist <- distEuclidean(data, res@centers)
        distSum <- res@clusinfo[,2]
        clusters <- res@cluster
        centers <- res@centers
        return(list(dist = dist, clusters = clusters, distSum =
                        distSum, centers=centers))
        
    }else if (init=="random"){
        set.seed(seed)
        T = t(rmultinom(nc,1,c(rep(1/k,k))))
        clusters <- apply(T, 1, function(x) which(x==max(x)))
        
    }else if (init=="hclust"){
        T   = matrix(0,nc,k)
        clusters = cutree(hclust(dist(data),method='ward'),k)
    }
    index <- rep(0,k)
    for (l in 1:k){
        ind1 <- which(clusters==l)
        T[ind1,l] = 1
        ctrs <- matrix(0,nrow=k,ncol=max(nc))
        partM <- as.matrix(dist(data[ind1,], upper=TRUE, diag=FALSE))
        if(length(ind1)==1) l <- 1
        else s <- which(rowSums(partM)==min(rowSums(partM)))[1]
        index[l] <- ind1[s]
    }
    centers <- data[index,]
    dist <- as.matrix(dist(data))[index,]
    distSum <- rowSums(dist)
    
    return(list(dist = dist, clusters = clusters, distSum = distSum, centers=centers))
}


##relabels the output with the higher class number according to the
##output with the lower class number. Output with lower class number
##must be stored first. 
relabel <- function(cl1, cl2, ctr1, ctr2){
    n1 <- length(cl1); n2 <- length(cl2)
    if(n1 > n2)
        stop("please put configuration with less observations on first place")
    t1 <- dim(ctr1)[1]; t2 <- dim(ctr2)[1]
    if(n1 < n2){
        indxN <- sample(1:n2, n1)
        cl2Old <- cl2
        ctr2Old <- ctr2
        cl2 <- cl2[indxN]
        indxT <- sample(1:t2, t1)
        ctr2 <- ctr2[sort(indxT),]
    }else{
        cl2Old <- cl2
        ctr2Old <- ctr2
    }
    
    clNrMin <- min(cl1,cl2)
    clNrMax <- max(cl1,cl2)
    
    dist <- ctrDist(ctr1, ctr2)
    max <- resort(dist)
    
    code <- paste(apply(max[c(2,1),],2,function(x)
        paste(x,collapse=" = ")), collapse=";")
    cl2New <- Recode(cl2Old, code)
    rownames(max) <- c("old cl2","new cl2")
    ctr2New <-ctr2Old[,max[2,]]
    
    return(list(recode=max, clusters=cl2New, centers=ctr2New))
}


ctrDist <- function(ctr1, ctr2){
    n1 <- dim(ctr1)[2]
    n2 <- dim(ctr2)[2]
    
    tb <- matrix(0, nrow=n1, ncol=n2)
    rownames(tb) <- 1:n1
    colnames(tb) <- 1:n2
    dist <- tb
    for(i in 1:n1)
        for(j in 1:n2)
            dist[i,j] <- sum((ctr1[,i]-ctr2[,j])^2)
    return(dist)
}


resort <- function(dist, func=min){
    clNrMin <- min(dim(dist))
    clNrMax <- max(dim(dist))
    if(is.null(rownames(dist)))
        rownames(dist) <- 1:clNrMin
    if(is.null(colnames(dist)))
        colnames(dist) <-  1:clNrMax
    
    names1 <- rownames(dist)
    names2 <- colnames(dist)
    
    max <- matrix(0, nrow=2, ncol=clNrMax)
    
    ind <- list()
    for(i in 1: (clNrMin)){
        if(is.null(dim(dist))) dist <- t(dist)
        ind[[i]] <- which(dist==func(dist, na.rm=TRUE), arr.ind = TRUE)
        if(dim(ind[[i]])[1]>1) ind[[i]] <- ind[[i]][1,]
        max[1,i] <- names1[ind[[i]][1]]
        max[2,i] <- names2[ind[[i]][2]]
        dist <- dist[-ind[[i]][1],-ind[[i]][2]]
        names1 <- names1[-ind[[i]][1]]
        names2 <- names2[-ind[[i]][2]]
    }
    if(!length(names1)==0)
        max[1,(clNrMax-length(names1)+1):clNrMax] <- names1
    if(!length(names2)==0)
        max[2,(clNrMax-length(names2)+1):clNrMax] <- names2
    mode(max) <- "numeric"
    if(sum(max[1,]==0)!=0)
        max[1,which(max[1,]==0)] <- max[2,!max[2,]%in%max[1,]]
    max <- max[,order(max[1, ]),drop=FALSE]
}


dist2centers <- function(data, centers){
    chf <- checkFormat(data)
    reg <-chf$reg
    data <- chf$data
    if(reg){
        nc <- dim(data)[2]
        nt <- dim(data)[1]
        res <- t(sapply(1:nc, function(x) colSums((data[,x]-centers)^2))/nt)
    }else{
        dataNew <- formatFuncy(data, format="Format3")
        res <- with(dataNew, makeSparse(Yin, Tin, isobs,
                                        t_all))
        out_indx <- match(sort(unique(data[,3])), dataNew$t_all)
        data <- as.matrix(res$longYin)
        isobs <- as.matrix(res$longIsobs)
        if(dim(centers)[2]==1)
            sumFct <- sum
        else
            sumFct <- colSums
        
        nc <- dim(data)[1]
        
        if(!is.null(time))
            centers <- centers[out_indx,]
        res <- t(sapply(1:nc, function(x){
            dist <- (data[x,]-centers)[isobs[x,]==1,]^2
            if(is.null(dim(dist)) | is.null(dim(centers)))
                return(dist/sum(isobs[x,]))
            else
                return(colSums(dist)/sum(isobs[x,]))
        }
                        )
                 )
        
    } 
    return(res)
}


makeCoeffs <- function(data, base=NULL, reg, dimBase, grid=NULL, pert, baseType){
    
    if(is.null(base)){
        tempBase <- makeBasis(baseType, grid, dimBase)$phi
        base <- svd(tempBase)$u
    }else{
        base <- base[,1:dimBase]
    }
    if(baseType=="fourier")
        dimBase <- dim(base)[2]
    
    if(reg){
        coeffs <- t((solve(t(base) %*% base + pert *
                               diag(dimBase))
                     %*%t(base))%*%t(data))
        fullBase <- base
    }else{
        curveIndx <- data[,1]
        timeIndx <- match(data[,3],grid)
        n <- max(curveIndx)
        fullBase <- base[timeIndx,  ]
        coeffs <- matrix(0,nrow=n,ncol=sum(dimBase))
        for (i in 1:n){
            if(is.null(dim(base)[1]))
                base <- t(t(base))
            basei <- fullBase[curveIndx==i,]
            yi <- data[curveIndx==i,2] 
            if(length(yi)>1){
                coeffs[i,] <- solve(t(basei) %*% basei + pert * diag(dimBase)) %*% t(basei) %*%yi
            }else{
                coeffs[i,] <- ((basei) * basei + pert )^(-1) * (basei)*yi
            }
        }
    }
    return(list(coeffs=coeffs, base=base, fullBase=fullBase, dimBase=dimBase))
}


makeBasis <- function(basis, time, nbasis){
    m <- length(time)
    switch(basis,
           "splines"={
               bObj <-  create.bspline.irregular(c(time[1],time[m]),
                                                 nbasis=nbasis,
                                                 norder=min(nbasis, 4))
           },           
           "exponential"={                  
               bObj <- create.exponential.basis(c(time[1], time[m]), nbasis=nbasis)
           },
           "fourier"={
               nbasis <- fourierWarn(nbasis)
               bObj <- create.fourier.basis(c(time[1], time[m]), nbasis=nbasis)
           },
           "power"={
               bObj <- create.power.basis(c(time[1], time[m]), nbasis=nbasis)
           })
    phi <- eval.basis(time, bObj)
    return(list(bObj=bObj, phi=phi))
}


fourierWarn <- function(dimBase){
    if(dimBase%%2==0){
        dimBase <- dimBase+1
        warning(paste("basis nr", dimBase,
                      "was used because it has to be odd for fourier !!!"))
    }
    return(dimBase)
}


calcTimeNr <- function(data, reg){
    if(reg){
        res <- data.frame(table(apply(data,2,length)))
    }else{
        res <- data.frame(table(data[,1]))
    }
    colnames(res) <- c("obs","nr of timepoints")
    return(res)
}


##outList Functions *************************************
accordance <- function(cls, ctrs=NULL, relabel=FALSE){

    if(relabel){
        if(is.null(ctrs))
            stop("A list of cluster centers is needed to adapt cluster labels to similar clusters.")
        rel <- relabelMethods(cls=cls, ctrs=ctrs)
        cls <- rel$allClusters
    }
    
    nrMod <- dim(cls)[2]
    nrObs <- dim(cls)[1]
    k <- length(table(cls))
    votedCluster <- accordance <- rep(0, nrObs)
    
    for(i in 1:nrObs){
        temp <- summary(as.factor(cls[i,]))
        indx <- which(temp==max(temp))[1]
        votedCluster[i] <- names(temp)[indx]
        accordance[i] <- temp[indx]
    }
    accordance <- round(accordance/nrMod,3)
    votedCluster <- as.numeric(votedCluster)
    
    return(list(votedCluster=votedCluster, accordance=accordance))
}


rIMethods <- function(methodNames=NULL, cls, trueCluster=NULL){
    nrMethods <- dim(cls)[2]
    rI <- matrix(0, nrow=nrMethods, ncol=nrMethods)
    if(!is.null(methodNames))
        colnames(rI) <- rownames(rI) <-  methodNames
    else
        colnames(rI) <- rownames(rI) <-  1:nrMethods
    counter1 <- counter2 <- 1
    for(i in 1:nrMethods){
        for(j in 1:nrMethods){
            a <- cls[,i]
            b <- cls[,j]
            if(counter1 == counter2 & !is.null(trueCluster))
                rI[counter1,counter2] <-
                    rI[counter1,counter2]+randIndex(table(a,trueCluster))
            else if(counter1 == counter2 & is.null(trueCluster))
                rI[counter1,counter2] <- NA
            else
                rI[counter1, counter2] <- rI[counter1,counter2]+randIndex(table(a,b))
            counter2 <- counter2+1
        }
        counter2 <- 1
        counter1 <- counter1+1
    }
    rI <- round(rI, 4)
    return(rI)
}


relabelMethods <- function(methodNames=NULL, cls, ctrs){
    nrMethods <- dim(cls)[2]
    
    rI <- rIMethods(methodNames=methodNames, cls=cls)
    
    rIOld <- rI
    if(is.null(methodNames))
        methodNames <- 1:nrMethods
    
    names <-  methodNames
    
    if(nrMethods>2){
        counter <- 1
        couple <- matrix(0, nrow=2, ncol=nrMethods-1)
        indx <- which(rI==max(rI, na.rm=T), arr.ind=T)[1,]
        growIndx <- sort(indx)
        growNames <- names[sort(indx)]
        rI <- rI[indx,-indx]
        couple[,1] <- names[indx]
        counter <- 2
        
        while(!is.null(dim(rI))){
            indx <- which(rI==max(rI, na.rm=T), arr.ind=T)
            growIndx <- c(growIndx,match(colnames(rI)[indx[1,2]],names))
            growNames <- c(growNames,colnames(rI)[indx[1,2]])
            couple[,counter] <- c(colnames(rI)[indx[1,2]],rownames(rI)[indx[1,1]])
            rI <- rIOld[growIndx,-growIndx]
            counter <- counter+1
        }
        indx <- which(rI==max(rI, na.rm=T), arr.ind=T)[1]
        couple[,counter] <- c(names[!names%in%growNames],names(rI)[indx])
        mat <- t(apply(couple, 1, function(x) match(x,names)))
        firstMeth <- mat[1,]
        secMeth <- mat[2,]
        rownames(couple) <- c("from","to")
    }else{
        couple <- NULL
        firstMeth <- 1
        secMeth <- 2
    }
    
    counter <- 1
    if(nrMethods>1){
        for(pioneer in secMeth){
            follower <- firstMeth[counter]
            if(length(table(cls[,pioneer]))<=length(table(cls[,follower]))){
                res <- relabel(cls[,pioneer], cls[,follower], ctrs[[pioneer]],
                               ctrs[[follower]])
                cls[,follower] <- res$clusters
                ctrs[[follower]] <- res$centers
            }
            counter <- counter+1
        }
        
        return(list(allClusters=cls, allCenters=ctrs, fromTo=couple))
    }
}


squareGrid <- function(x, round=FALSE){
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    n <- length(factors)
    if(n==2 & round & x!=2)
        squareGrid(x+1)
    else{
        if(n%%2==0)
            indx1 <- n/2
        else
            indx1 <- floor(n/2)+1
        n1 <- factors[indx1]
        n2 <- x/n1
        return(c(n1,n2))
    }
}


longCtrY <- function(Yin, Tin, time, isobs, mean){
    indx_tin <- makeIndex(Tin, time, isobs)
    Mufcn <-t(apply(indx_tin, 1, function(x) mean[x]))
    ctrY <- (Yin-Mufcn)
    return(ctrY)
}


makeIndex <- function(Tin, time, isobs){
    if(is.null(dim(Tin)))
        indx_tin <- match(Tin, time)
    else
        indx_tin <- t(apply(Tin, 1, function(x) match(x,time)))
    if(sum(!isobs)> 0)
        indx_tin[!as.matrix(isobs)] <- NA
    return(indx_tin)
}


longCov <- function(ctrY, Tin, isobs, time){
    N <- rowSums(isobs)
    cov.raw <-numeric(sum(N*(N-1)/2))     
    cov.timeIn<-matrix(0,nrow=sum(N*(N-1)/2),ncol=2)
    var.raw <- rep(0, sum(N))
    var.timeIn <- rep(0, sum(N))
    count <- count2 <-1
    for (i in 1:length(N)){
        if(N[i]>1){
            yi <- ctrY[i,1:N[i]]
            ti <- Tin[i,1:N[i]]
            for (j in 1: N[i]){
                var.raw[count] <-(yi[j])^2
                var.timeIn[count] <- ti[j]
                count <- count+1
                if(j==N[i])
                    next
                for (l in (j+1): N[i]){
                    cov.raw[count2] <- (yi[j])*(yi[l])
                    cov.timeIn[count2,1] <- ti[j] 
                    cov.timeIn[count2,2] <- ti[l]
                    count2 <- count2+1
                }
            }
        }
    }
    return(list(cov.raw=cov.raw, cov.timeIn=cov.timeIn, var.raw=var.raw, var.timeIn=var.timeIn))
}

label2lowerk <- function(cluster){
    z <- as.factor(cluster)
    z <- factor(z, levels = levels(z), labels =
                    1:length(unique(z)))
    z <- as.numeric(z)
    k <- max(z)
    return(list(cluster=z, k=k))
}


getUniCl <- function(id, clusters, reduce=TRUE){
    if(reduce){
        pos <- which(!duplicated(id))
        cl <- clusters[pos]
    }else{
        cl <- rep(clusters, table(id))
    }
    return(cl)
}


makeClMat <- function(dist2centers, n=2){
    
    r <- t(matrix(apply(dist2centers, 1, rank, ties.method = "random"), 
                  nrow = ncol(dist2centers)))
    z <- list()
    for (k in 1:n) z[[k]] <- apply(r, 1, function(x) which(x == k))
    
    cldist <- cbind(dist2centers[cbind(1:nrow(dist2centers), z[[1]])],
                    dist2centers[cbind(1:nrow(dist2centers), z[[2]])])
    return(cldist)
}




