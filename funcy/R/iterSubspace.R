#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

iterSubspace <- function(data, k, reg, regTime,
                         funcyCtrlMbc, fpcCtrl, simplif){

    ##control parameters
    dimBase <- funcyCtrlMbc@dimBase
    baseType <- funcyCtrlMbc@baseType
    thd <- funcyCtrlMbc@thd
    maxit <- funcyCtrlMbc@maxit
    seed <- funcyCtrlMbc@seed
    init <- funcyCtrlMbc@init
    nrep <- funcyCtrlMbc@nrep
    flexDim <- funcyCtrlMbc@flexDim
    
    ##Reformat data
    res <- formatFuncy(data=data, format="Format3", regTime=regTime)
    if(reg)
        coeffFct <- integrateCoeffs
    
    Yin <- res$Yin; Tin <- res$Tin; N <- res$N; isobs <- res$isobs;
    t_all <- res$t_all

    oldK <- k
    
    ##Projection function
    if(baseType=="eigenbasis"){
        nrMaxTime <- fpcCtrl@nrMaxTime
        projFct <- function(Yin=Yin, Tin=Tin, Tout=t_unique, isobs=isobs,
                            dimBase=dimBase, fpcCtrl=fpcCtrl,
                            baseType=baseType, reg=reg){
            fpc(Yin=Yin, Tin=Tin, Tout=t_unique, isobs=isobs,
                dimBase=dimBase, fpcCtrl=fpcCtrl, reg=reg)
        }
        ##fpcCtrl parameters
        if(fpcCtrl@coeffsCalc=="estimate")
            coeffFct <- estimateCoeffs
        else if(fpcCtrl@coeffsCalc=="integrate")
            coeffFct <- integrateCoeffs
    }else
        projFct <- anyBasis
    
    ##Some infos
    t_unique <- sort(unique(t_all))
    nt <- length(t_unique)
    nc <- dim(Yin)[1] 
    m <- dim(Yin)[2] 
    id <- as.numeric(1:nc)

    ##Keep initial data
    YinOld <- Yin; TinOld <- Tin; isobsOld <- isobs
    indx_tinOld <- t(apply(TinOld, 1, function(x) match(x,t_unique)))
    
    ##Coefficients for initial clustering
    res <- projFct(Yin=Yin, Tin=Tin, Tout=t_unique, isobs=isobs,
                   dimBase=dimBase, fpcCtrl=fpcCtrl,
                   baseType=baseType, reg=reg)

    plotParams <- res$plotParams
    Yin <- res$Yin;
    Tin <- res$Tin;
    isobs <- res$isobs;
    coeffs <- res$coeffs
    ctrYin <- res$ctrY
    average <- res$average
    Delta <- res$Delta

    ##Check which center function should be used. 
    if(average)
        ctrFct <- match.fun("avCtrY")
    else
        ctrFct <- match.fun("longCtrY")

    ##INITIAL CLUSTERING --------------------------------------------
    res <- initClust(data=coeffs[,1:dimBase], k=k, init=init,
                     seed=seed, nrep=nrep)
    class <- initClass  <-  res$clusters
    
    reCluster <- matrix(0, nc, maxit)
    Err <- groupDimBase <- matrix(NA, maxit+1, k)
    groupMeansAll <- list()

    ##START
    ##ITERATION-------------------------------------------------------
    for(kk in 1:maxit){
        lastCl <- class
        clusts <- split(1:nc, lastCl, drop=TRUE)
        k <- length(clusts)
        
        res <- buildGroups(nt=nt,
                      dimBase=dimBase,
                      k=k,
                      clusts=clusts,
                      thd=thd,
                      isobs=isobsOld,
                      Tin=TinOld,
                      Yin=YinOld,
                      t_unique=t_unique,
                      average=average,
                      flexDim=flexDim,
                      fpcCtrl=fpcCtrl,
                      projFct=projFct,
                      baseType=baseType,
                      reg=reg)
        
        groupDimBase[kk, 1:k] <- res$groupDimBase
        groupMeans <-  res$groupMeans
        groupMeansAll[[kk]] <- groupMeans
        groupBase <- res$groupBase
        groupCoeffs <- res$groupCoeffs
        groupEigval <- res$groupEigval
        groupDelta <- res$groupDelta
        groupevar <- res$groupevar
        Err[kk, 1:k] <- res$groupErr
        errs <- 10000000 * matrix(1, nc, k)
        k <- res$k
        
        ##Project all curves (if !simplif except the ones of current
        ##cluster) into each space and assign each curve where error
        ##is minimal
        for(clusterRunner in 1:k){
            if(simplif)
                idOut <- id
            else
                idOut = id[-clusts[[clusterRunner]]]
            
            currentBase <- groupBase[[clusterRunner]]
            currentMean <- groupMeans[,clusterRunner]
            if(is.null(currentBase)){
                warning("The cluster number was reduced!")
                next
            }
            tempIsobs <- isobs[idOut,]
            tempTin <- Tin[idOut,]
            tempIndx_tin <- makeIndex(Tin=TinOld[idOut,],
                                      time=t_unique, isobs=isobsOld[idOut,])
            tempYin <- Yin[idOut,]
            tempCtrY <- ctrFct(Yin=tempYin, Tin=tempTin, time=t_unique, isobs=tempIsobs,
                               mean=currentMean)
            tempDelta <- Delta[idOut,]
            npcs <- groupDimBase[kk,clusterRunner]
            
            ##The functions to calculate the coefficients differ
            ##depending on basis in use
            if(baseType=="eigenbasis"){
                coeffs <- coeffFct(ctrY=tempCtrY,
                                   dimBase=npcs,
                                   base=currentBase,
                                   Delta=tempDelta,
                                   indx_tin=tempIndx_tin,
                                   eigval=groupEigval[[clusterRunner]],
                                   average=average,
                                   evar=groupevar[clusterRunner])$coeffs
            }else{
                if(!reg)
                    dataCtrY <- formatFuncy(data=list(Yin=tempCtrY,
                                                isobs=tempIsobs,
                                                Tin=tempTin),
                                            format="Format1"
                                           )
                else
                    dataCtrY <- tempCtrY
                coeffs <- makeCoeffs(data=dataCtrY,
                                     base=currentBase,
                                     reg=reg,
                                     dimBase=npcs,
                                     grid=t_unique,
                                     pert=0.05,
                                     baseType=baseType)$coeffs

            }
            
            err <- projError(coeffs=coeffs,
                             ctrY=tempCtrY,
                             base=currentBase,
                             dimBase=npcs,
                             isobs=tempIsobs,
                             indx_tin=tempIndx_tin,
                             average=average)
            
            errs[idOut, clusterRunner] <- sqrt(err)

            ##project each curve in current cluster to other
            ##spaces-------------------------
            if(!simplif){
                for(curveRunner in clusts[[clusterRunner]]){
                    coeffIndx <- which(curveRunner==clusts[[clusterRunner]])
                    inds <- clusts[[clusterRunner]][-coeffIndx]
                    ncurve <- length(inds)
                    if(ncurve == 0)
                        break
                    tempisobs <- isobsOld[inds,]
                    tempTin <- TinOld[inds,]
                    tempYin <- YinOld[inds,]

                    res <- projFct(Yin=tempYin,
                                   Tin=tempTin,
                                   Tout=t_unique,
                                   isobs=tempisobs,
                                   dimBase=npcs,
                                   fpcCtrl=fpcCtrl,
                                   baseType=funcyCtrlMbc@baseType,
                                   reg=reg)
                    Ypred <- auxmufcn <- t(res$meanfcn)
                    tempBase <- res$base 
                    isobs_i <- t(isobs[curveRunner,])
                    index_i <- t(indx_tinOld[curveRunner, 1:N[curveRunner]])
                    base_i <- t(t(tempBase))
                    Y_i <- Yin[curveRunner,]
                    ctrCurve_i <- t(ctrYin[curveRunner,])
                    coeffs_i <-
                        t(groupCoeffs[[clusterRunner]][coeffIndx,])
                    
                    a <-
                        (ctrCurve_i[isobs_i==1]-(base_i%*%t(coeffs_i))[index_i,])^2
                    
                    
                    err <- projError(coeffs=coeffs_i,
                                     ctrY=ctrCurve_i,
                                     base=base_i,
                                     dimBase=npcs,
                                     isobs=isobs_i,
                                     indx_tin=index_i,
                                     average=average)
                    errs[curveRunner, clusterRunner] <- sqrt(err)
                    
                }
            }
        }#end cluster runner
        
        class <- apply(errs, 1, function(x) which(x==min(x)))
        if(!is.null(dim(class))){
            class <- class[1,]
            warning("Class probability is equal")
        }
        reCluster[,kk] <- class   
        iditerSubspace <- class
        true <- 0
        ##Break up conditions ---------------------------------------
        ##1. If class assignment is the same as for starting class true=1
        if(all(initClass == iditerSubspace)){
            true <- 1
            ##2. If class assignment repeats at any iteration true=2
        }else{
            if(kk>2){
                for(j in 1:(kk-2)){
                    if(all(reCluster[,j]==iditerSubspace)){
                        true <- 2
                        break
                    }
                }
            }
        }
        
        if(true == 1){
            idselect <- 1
            break
        }else if(true == 2){
            sumSSE <- apply(Err, 1, function(x) sum(na.omit(x)))
            idselect <- which(sumSSE[1:(kk-1)]==min(sumSSE[1:(kk-1)]))
            idselect <- idselect[1]
            if(idselect==1){
                iditerSubspace <- initClass
            }else{
                iditerSubspace <- reCluster[,idselect-1]
            }
            break
        }
        ##3. If two class assignments in a row stay the same
        idmove <- (class == lastCl)
        nmove <- length(which(idmove == 0))
        if(nmove==0){
            idselect <- kk
            break
        }
        
        ##4. If no converging after maxit iterations
        if(kk == maxit){
            idselect <- kk
            print(paste('Warning- iterSubspace no covergence after ', maxit, ' iterations.'))
        }
    }
    ##-----------------------------------------------------
    groupDimBase <- (groupDimBase)[idselect,]
    centers <- groupMeansAll[[idselect]]
    groupErr <- Err[idselect,]
    nrIter <- idselect
    
    
    return(list(dimBase=dimBase,
                groupDimBase=groupDimBase, groupBase=groupBase, groupMeans=groupMeans,
                groupErr=groupErr, initClass=as.numeric(initClass),
                cls=iditerSubspace, ctrs=centers, plotParams=plotParams,
                grid=t_unique, nrIter=nrIter))
    
    
}
    

buildGroups <-
    function(nt, dimBase, k, clusts, thd, isobs, Tin, Yin, t_unique,
             average, flexDim, fpcCtrl, projFct,
             baseType=baseType, reg=reg){
        if(!is.null(fpcCtrl))
            nrMaxTime <- fpcCtrl@nrMaxTime
        else
            nrMaxTime <- 50
        
        nt <- min(nt, nrMaxTime)
        nt.new <- length(t_unique)
        N <- rowSums(isobs);
        m <- dim(Tin)[2]
        nc <- dim(Yin)[1]
        
        groupDimBase <- rep(NA, k)
        groupMeans <- matrix(NA, nt.new, k)
        groupBase <- list()
        groupEigval <- list()
        groupCoeffs <- list()
        ctrY <-  list()
        groupDelta <- list()
        groupevar <- matrix(NA, k, 1)
        groupErr <- rep(NA, k)
        
        for(l in 1:k){
            id <- clusts[[l]]
            ncurve <- length(id)
            if(ncurve <= 1){
                groupErr[l] <- -99
                groupDimBase[l] <- -99
                k <- k-1
            }else if(ncurve > 1){
                tempisobs <- isobs[id,];
                tempTin <- Tin[id,]
                tempYin <- Yin[id,];
                indx_tin <- makeIndex(tempTin, t_unique, tempisobs)
                
                res1 <- projFct(Yin=tempYin,
                                Tin=tempTin,
                                Tout=t_unique,
                                isobs=tempisobs,
                                dimBase=dimBase,
                                fpcCtrl=fpcCtrl,
                                baseType=baseType,
                                reg=reg)
                
                if(!is.null(res1$eigval))
                    groupEigval[[l]] <- res1$eigval 
                groupCoeffs[[l]] <- res1$coeffs
                groupBase[[l]] <- res1$base
                groupMeans[,l] <- res1$meanfcn
                groupevar[l] <- res1$evar
                ctrY[[l]] <- res1$ctrY
                groupDelta[[l]] <- res1$Delta
                tempisobs <- res1$isobs
                
                if(flexDim){
                    res <- findBasisNr(coeffs=groupCoeffs[[l]],
                                       ctrY=ctrY[[l]],
                                       base=groupBase[[l]],
                                       Delta=groupDelta[[l]],
                                       isobs=tempisobs,
                                       average=average,
                                       indx_tin=indx_tin,
                                       thd=thd,
                                       eigval=groupEigval[[l]],
                                       evar=groupevar[l])
                    groupDimBase[l] <- res$nrBase
                    groupErr[l] <- res$Err
                    groupCoeffs[[l]] <- res1$coeffs[,1:groupDimBase[l]]
                    groupBase[[l]] <- res1$base[,1:groupDimBase[l]]
                    if(groupDimBase[l]==1){
                        groupBase[[l]] <-t(t(groupBase[[l]]))
                        groupCoeffs[[l]] <- t(t(groupCoeffs[[l]]))
                    }
                }else{
                    groupDimBase[l] <- dimBase
                    groupErr[l] <-
                        sum(projError(coeffs=groupCoeffs[[l]],
                                      ctrY=ctrY[[l]],
                                      base=groupBase[[l]],
                                      dimBase=dimBase,
                                      isobs=tempisobs,
                                      indx_tin=indx_tin,
                                      average=average))
                    
                }
                
            }
        }
        return(list(groupDimBase=groupDimBase, groupCoeffs=
                        groupCoeffs, groupMeans=groupMeans,
                    groupBase=groupBase, groupevar=groupevar,
                    groupErr=groupErr, groupEigval=groupEigval,
                    groupDelta=groupDelta, k=k, grid=t_unique))
    }


projError <- function(coeffs, ctrY, base, dimBase, isobs, indx_tin, average){
    nc <- dim(ctrY)[1]
    m <- dim(base)[2]
    nt <- dim(base)[1]
    N <- rowSums(isobs)
    err <- rep(0,nc)

    if(dimBase==0)
        proj <-  matrix(mean, nrow=nt, ncol=nc)
    else if(dimBase==1){
        mat <- matrix(base[,1], nrow=nt, ncol=nc)
        proj <- t(mat)*coeffs[,1]
    }else{
        proj <- coeffs[,1:dimBase]%*%t(base[,1:dimBase])
    }
    if(average){
        err <- rowSums((ctrY-proj*isobs)^2)
    }else{
        for(i in 1:nc)
            err[i] <-
                sum((ctrY[i,isobs[i,]==1]-na.omit(proj[i,indx_tin[i,]]))^2)
        
    }
    err <- err/N
    return(err)
}


findBasisNr <- function(coeffs, ctrY, base, Delta, isobs, average,
                        indx_tin, thd, eigval=NULL, evar=NULL){
    nt <- dim(ctrY)[2]
    maxDim <- min(nt,dim(base)[2])
    tempsse <- Ratio_sse <- tempsse2 <- rep(0, nt+1)
    for(nk in 1:nt){
        tempsse[nk] <- sum(projError(coeffs=coeffs,
                                     ctrY=ctrY,
                                     base=base,
                                     dimBase=nk,
                                     isobs=isobs,
                                     indx_tin=indx_tin,
                                     average=average))
        if(nk == 1){
            Ratio_sse[nk] <- 1
        }else{
            Ratio_sse[nk] <- (tempsse[nk-1]-tempsse[nk])/tempsse[1]
        }
        if(Ratio_sse[nk] >= thd & nk != maxDim){
            nrBase <- nk
            Err <- tempsse[nk]
        }else if(Ratio_sse[nk] < thd | nk == maxDim){
            return(list(nrBase=nrBase, Err=Err))
        }
        
    }
}

anyBasis <- function(Yin, Tin, Tout=NULL, isobs, dimBase, plot=NULL,
                     fpcCtrl=NULL, baseType, reg=reg){

    t_unique <- sort(unique(c(Tout, Tin)))
    nt <- length(t_unique)
    nc <- dim(Yin)[1]
    tempDelta <- t(apply(Tin, 1, diff))
    Delta <- cbind(rep(mean(diff(t_unique)), nc),tempDelta)
    Delta <- cbind(rep(0, nc),tempDelta)
    indx_tin <- makeIndex(Tin, t_unique, isobs)
    ctrY <- Yin
    
    if(!reg){
        dat <- formatFuncy(data=list(Yin=Yin, isobs=isobs, Tin=Tin),
                           format="Format1")
        res <- makeCoeffs(data=dat, reg=FALSE, dimBase=dimBase,
                          grid=t_unique, pert=0.01, baseType=baseType)
    }else{
        res <- makeCoeffs(data=Yin, reg=TRUE, dimBase=dimBase,
                          grid=t_unique, pert=0.01,
                          baseType=baseType)
        
    }
    ##PART 2 - EIGENVALUES AND VECTORS------------------------------
    ##Calculate Eigenvalues and vectors-----------------------
    base <- fullBase <- res$base
    coeffs <- res$coeffs
    yReg <- fullBase%*%t(coeffs)
    mean <- apply(yReg, 1, mean)
    evar <- sum(projError(coeffs=coeffs, ctrY=ctrY, base=fullBase, dimBase=dimBase,
                          isobs=isobs, indx_tin=indx_tin, average=FALSE))
    eigval <- rep(1,dimBase)

    return(list(Yin=Yin, Tin=Tin, fullBase=fullBase, isobs=isobs, ctrY=ctrY,
                Delta=Delta, meanfcn=mean, covfcn=NULL,
                base=base, eigval=eigval, coeffs=coeffs,
                varprop=NULL, evar=evar, average=FALSE))
}
