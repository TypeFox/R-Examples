#  Original  Copyright (C) 2008  Jie Peng and Hans-Georg Mueller
#  Modiefied Code: Copyright (C) 2011-2015 Christina Yassouridis
#  
#


##added function by Yassouridis
distclust <- function(data, k, reg, regTime, simpleCtrl,
                     fpcCtrl, method="pam")
    
    {
        ##control parameteres
        dimBase <- simpleCtrl@dimBase
        baseType <- simpleCtrl@baseType
        
        ##fpcCtrl parameters
        if(!is.null(fpcCtrl)){
            average <- fpcCtrl@average
            smoothFct1D <- match.fun(fpcCtrl@sm1Dim)
        }else{
            average <- TRUE
            smoothFct1D <- match.fun("sm.regression")
        }
                
        ##Reformat data
        res <- formatFuncy(data=data, format="Format3")
        
        Yin <- res$Yin; Tin <- res$Tin; N <- res$N; isobs <-
            res$isobs; t_all <- res$t_all

        ##Some Infos
        t_unique <- sort(unique(t_all))
        nt <- length(t_unique)
        nc <- dim(Yin)[1]
        
        ##Get basis for initial clustering
        if(baseType=="eigenbasis"){
            res <- fpc(Yin=Yin, Tin=Tin, isobs=isobs, dimBase=dimBase,
                      Tout=t_unique, fpcCtrl=fpcCtrl, reg=reg)
            mean <- res$meanfcn
            evar <- res$evar/diff(range(t_unique))
            base <- res$base
            plotParams <- res$plotParams
            Gamma  <- diag(res$eigval[1:dimBase])
        }else{
            res1 <- fitfclust(data=data, dimBase=dimBase, reg=reg, regTime=regTime, k=1,
                              h=1, hard=TRUE, baseType=baseType)
            base <- res1$base
            prms <- res1$parameters
            evar <- prms$sigma^2
            dimBase <- dim(base)[2]
            mean <- base%*%(prms$lambda.zero + prms$Lambda%*%t(prms$alpha))
            gprod <- res1$vars$gprod
            ind <- do.call(rbind,lapply(1:nc, function(x) diag(dimBase)))
            gsvd <- svd(gprod %*% ind/nc)
            Gamma <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u)
            plotParams <- NULL
        }

        ##*********************************************************
        ##Part 3: Clustering with distance
        cdist <- PCADist(Yin=Yin, Tin=Tin, isobs=isobs, time=t_unique, Gamma=Gamma, base=base,
                         mean=mean, evar=evar, K=dimBase)
        
        if(method=="hclust"){
            hc <- hclust(as.dist(cdist), "ave")
            cls <- cutree(hc, k=k)
            index <- rep(0,k)
            for(i in 1:k){
                ind1 <- which(cls==i)
                partM <- cdist[ind1,ind1]
                if(length(ind1)==1) l <- 1
                else l <- which(rowSums(partM)==min(rowSums(partM)))[1]
                index[i] <- ind1[l]
            }
            allTime <- rowSums(isobs[index,])
            
            dist <- cdist[,index]
        }else if(method=="pam"){
            hc <- NULL
            cl <- pam(cdist, k = k, diss=TRUE)
            cls <- cl$clustering
            dist <- cdist[,cl$medoid]
            allTime <- rowSums(isobs[cl$medoids,])
        }
        
        if(average){
            res <-makeSparse(Yin=Yin, Tin=Tin, isobs=isobs, time=t_unique)
            Yin <- res$longYin; Tin <- res$longTin; isobs <- res$longIsobs
        }
        ctrs <- smoothCenters(cluster=cls, Yin=Yin, Tin=Tin, isobs=isobs, time=t_unique,
                             average=average, smoothFct1D=smoothFct1D)
        props <- as.numeric(table(cls)/length(cls))
        dist2centers <- dist/allTime
      
        ##Output arguments
        res <- list(data=data, reg=reg, k=k, dimBase=dimBase,
                    cluster=cls, centers=ctrs, props=props,
                    dist2centers=dist2centers, plotParams=plotParams, grid=t_unique)
        
       return(res)
    }



PCADist<-function(Yin, Tin, isobs, time, Gamma, base,
                  mean, evar, K){
    nc <- dim(Tin)[1]
    m <- dim(Tin)[2]
    nt <- length(time)
    N <- rowSums(isobs)
    h1 <- as.vector(t(isobs))
    y <- Yin[h1==1] 
    t <- Tin[h1==1]
    indx_tin <- makeIndex(Tin, time, isobs)
    Mufcn <-t(apply(indx_tin, 1, function(x) mean[x]))
    result <- matrix(0,nc,nc)
  
    ##Pairwise conditional L2 distance 
    for (i in 1:(nc-1)){
        for (j in (i+1):nc){
            Ti <- Tin[i,1:N[i]]
            X <- Yin[i,1:N[i]]
            mux <- Mufcn[i,1:N[i]]
            if(sum(isobs[i,]==1)==1)
                Phix <- t(base[(indx_tin[i,][isobs[i,]==1]),][1:K])
            else
                Phix <- base[(indx_tin[i,][isobs[i,]==1]),][,1:K]
            Tj <- Tin[j,1:N[j]]
            Y <- Yin[j,1:N[j]]
            muy <- Mufcn[j,1:N[j]]
       
            if(sum(isobs[j,]==1)==1)
                Phiy <- t(base[(indx_tin[j,][isobs[j,]==1]),][1:K])
            else
                Phiy <- base[(indx_tin[j,][isobs[j,]==1]),][,1:K]
            result[i,j] <- CDist(Gamma, evar, Phix, Phiy, mux,
                                 muy, X, Y)[[1]]
            result[j,i] <- result[i,j]
            
        }
    }
    return(result)
}

CDist <- function(Gamma, evar, Phix, Phiy, mux, muy, X, Y){          
    ##conditional L2 distance for two subjects given mu,phi evaluated
    ##at observed points;distance based on conditional second moments and
    ##conditional PC scores of subjects X and Y
    Nx <- length(X)
    Ny <- length(Y)
       
    ##conditional mean
    cmeanx <- Gamma%*%t(Phix)%*%solve(Phix%*%Gamma%*%t(Phix)+evar*diag(Nx))%*%(X-mux)
    cmeany <- Gamma%*%t(Phiy)%*%solve(Phiy%*%Gamma%*%t(Phiy)+evar*diag(Ny))%*%(Y-muy)
    
    ##conditional covariance
    ccovx <-
        Gamma-Gamma%*%t(Phix)%*%solve(Phix%*%Gamma%*%t(Phix)+evar*diag(1,Nx))%*%Phix%*%Gamma
    ccovy <-
        Gamma-Gamma%*%t(Phiy)%*%solve(Phiy%*%Gamma%*%t(Phiy)+evar*diag(1,Ny))%*%Phiy%*%Gamma
    
    ##conditional second moments 
    cm2x <- diag(ccovx)+cmeanx^2
    cm2y <- diag(ccovy)+cmeany^2
    cm2.distance <- sum(cm2x+cm2y-2*cmeanx*cmeany)
    return(list(cm2.distance,cmeanx,cmeany))
}


smoothCenters <- function(cluster, Yin, Tin, isobs, time,
                        average, smoothFct1D){
    if(average)
        meanFct <- match.fun("avMean")
    else
        meanFct <- match.fun("longMean")
    
    cls <- split(1:length(cluster),cluster)
    nrCl <- lapply(cls,  length)
    Yin <- as.matrix(Yin); Tin <- as.matrix(Tin); isobs <- as.matrix(isobs)
    centers <- lapply(1:max(cluster), function(x){
        rawMean <- meanFct(Yin[cls[[x]],],
                           Tin[cls[[x]],],
                           isobs[cls[[x]],])
        hmu.cv <- h.select(rawMean$time,rawMean$mu, method="cv")
        smMean <- smoothFct1D(x=rawMean$time,
                              y=rawMean$mu,
                              h=hmu.cv,
                              eval.points=time,
                              weights=rep(1,length(rawMean$time)),
                              poly.index=1,
                              display="none")$estimate
        if(sum(is.na(smMean))==0)
            return(smMean)
        else if(sum(is.na(smMean))>0){
            hmu <- mean(diff(rawMean$time))
            smMean <- smoothFct1D(x=rawMean$time, y=rawMean$mu,
                                  h=hmu,eval.points=time,
                                  display="none")$estimate
        }
    }
                      )
    centers <- do.call(cbind, centers)
    return(centers)
}
