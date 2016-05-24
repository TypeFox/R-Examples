#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  Original code in matlab, J.-M. Chiou and P.-L. Li J. R. Statist. Soc. B, Volume 69 (2007), 679-699
#

##does the same as fpc but is exported
fpca <- function(data=data, dimBase=4, fpcCtrl=NULL){
    
    chf <- checkFormat(data, reformat=FALSE)
    data <- chf$data
    reg <- chf$reg
    res <- formatFuncy(data, format="Format3")
    
    Yin <- res$Yin; Tin <- res$Tin; isobs=res$isobs
    
    fpcCtrl <- fpcCtrlCheck(fpcCtrl=fpcCtrl,
                            data=data,
                            reg=reg)
    
    res <- fpc(Yin=Yin, Tin=Tin, Tout=NULL, isobs=isobs,
               dimBase=dimBase, fpcCtrl=fpcCtrl, reg=reg)
    yreg <- res$base%*%t(res$coeffs)
    time <- res$plotParams$t_unique
    
    return(list(yreg=yreg, time=time, meanfcn=res$meanfcn, covfcn=res$covfcn,
                base=res$base, eigval=res$eigval, coeffs=res$coeffs,
                varprop=res$varprop))

}


fpc <- function(Yin=NULL, Tin=NULL, Tout=NULL, isobs=NULL, dimBase=dimBase, fpcCtrl=fpcCtrl, reg=reg){
    ##fpcCtrl parameters
    if(fpcCtrl@coeffsCalc=="estimate")
        coeffFct <- estimateCoeffs
    else if(fpcCtrl@coeffsCalc=="integrate")
        coeffFct <- integrateCoeffs
    hmu <- fpcCtrl@h1Dim
    hcov <- fpcCtrl@h2Dim
    nrMaxTime <- fpcCtrl@nrMaxTime
    average <- fpcCtrl@average
    smoothFct1D <- match.fun(fpcCtrl@sm1Dim)
    smoothFct2D <- match.fun(fpcCtrl@sm2Dim)
    
    ##Some Infos
    t_unique <- sort(unique(c(Tout, as.vector(Tin))))
    nt <- length(t_unique)
    nc <- dim(Yin)[1]
    
    tempDelta <- t(apply(Tin, 1, diff))
    Delta <- cbind(rep(mean(diff(t_unique)), nc),tempDelta)
    Delta <- cbind(rep(0, nc),tempDelta)
    indx_tin <- makeIndex(Tin=Tin, time=t_unique, isobs=isobs)

    if(average){
        meanFct <- avMean 
        covFct <- avCov 
        ctrFct <- avCtrY
        if(class(Yin)=="matrix"){
            sp <- makeSparse(Yin=Yin, Tin=Tin, isobs=isobs, time=t_unique)
            Yin <- sp$longYin
            isobs <-sp$longIsobs
            Tin <- sp$longTin
            Delta <- sp$longDelta
        }
        
    }else{
        meanFct <- match.fun("longMean")
        covFct <- match.fun("longCov")
        ctrFct <- match.fun("longCtrY")
    }   
    ## I. SMOOTH MU-FCT resp. RAW DATA ------------------------------
    ##raw mean
    meanRes <- meanFct(Yin=Yin, Tin=Tin, isobs=isobs)
    mean.raw <- meanRes$mu
    mean.timeIn <- meanRes$time

    ##smoothed mean
    isNa <- counter <- 0
    while(isNa==0){
        counter <- counter+1
        mu.smoothed <- smoothFct1D(x=mean.timeIn, y=mean.raw, h=hmu,
                                   eval.points=t_unique,
                                   weights=rep(1,length(mean.timeIn)),
                                   poly.index=1,
                                   display="none")$estimate

        if(sum(is.na(mu.smoothed))==0)
            isNa <- 1
        else{
            hmu <- hmu+counter*0.1*hmu
            warning("Another smoothing parameter", round(hmu,3),
                    " is used for mean smoothing.")
        }
        if(counter==10)
            stop("Smoothing problem for mean! Please choose another sm1Dim in fpcCtrl!")
    }
    
    ##II. SMOOTH COVARIANCE and VARIANCE-----------------------------
    ##Raw covariance
    ctrY <- ctrFct(Yin=Yin, Tin=Tin, time=t_unique, isobs=isobs, mean=mu.smoothed)
    covRes <- covFct(ctrY=ctrY, Tin=Tin, isobs=isobs, time=t_unique)
    
    ##Smoothing input
    cov.raw <- covRes$cov.raw
    cov.timeIn <- covRes$cov.timeIn
    cov.timeOut <- makeGrid(nt=nt, time=t_unique)
    var.timeIn <- covRes$var.timeIn
    var.raw <- covRes$var.raw
    t.red <- t_unique
    ntold <- nt
    
    ##Reduce evaluation points if unique time points is too long
    if(nt > nrMaxTime){
        t.red <- seq(min(t_unique), max(t_unique), length.out=nrMaxTime)
        cov.timeOut <- makeGrid(nt=nrMaxTime, time=t.red)
        ntold <- nt
        nt <- nrMaxTime
    }
    ##smoothed
    ##covariance----------------------------------------------
    hcovNew <- min(diff(t_unique))
    counter <- 0
    isNa <- 0
    while(isNa==0){
        counter <- counter+1
        cov.smoothed <- suppressWarnings(smoothFct2D(x=cov.timeIn,
                                                     y=cov.raw,
                                                     h=hcov,
                                                     eval.points=cov.timeOut,
                                                     weights=rep(1,length(cov.raw)),
                                                     poly.index=c(1,2),
                                                     eval.grid=FALSE,
                                                     display="none")$estimate)

        if(sum(is.na(cov.smoothed))==0)
            isNa <- 1
        else{
            hcovNew <- hcovNew+counter*0.1*hcovNew
            hcov <- rep(hcovNew,2)
            warning("Another smoothing parameter", round(hcov,3)[1],
                    " is used for covariance smoothing.")
        }
        if(counter==10)
            stop("Smoothing problem! Please choose another sm2Dim in fpcCtrl.")
    }
    
    ##Convert covariance to matrix
    covest <- makeCovMat(cov=cov.smoothed, nt=nt)
    
    ##Calculate Error variance ChiouLi - pc2.m--------------------------
    trmat <- matrix(c(1/sqrt(2), 1/sqrt(2), -1/sqrt(2),  1/sqrt(2)),
                    nrow=2, byrow=T)
    trxin <- trmat%*%t(cov.timeIn);
    trxin <- t(trxin)
    trxou <- matrix(0, ntold, 2);
    trxou[,1] <- trxou[,2] <- sqrt(2)*t_unique
    
    ##Smooth raw covariance on transformed time points as in and
    ##output
    hcovNew <- min(diff(trxin[,1]))
    counter <- 0
    isNa <- 0
    
    while(isNa==0){
        tildeG <- suppressWarnings(smoothFct2D(x=trxin, y=cov.raw,
                                               h=hcov, eval.points=trxou,
                                               weights=rep(1,length(cov.raw)),
                                               poly.index = c(1,2),
                                               eval.grid=FALSE,
                                               display="none")$estimate)
        if(sum(is.na(tildeG))==0)
            isNa <- 1
        else{
            hcovNew <- hcovNew+counter*0.1*hcovNew
            hcov <- rep(hcovNew,2)
            warning("Another smoothing parameter", round(hcov,3)[1],
                    " is used for covariance smoothing on transformed points.")
        }
        if(counter==10)
            stop("Smoothing problem!")
    }
    
    ##smooth raw variance on usual time points
    hatV <- smoothFct1D(x=var.timeIn, y=var.raw, h=hmu,
                        eval.points=t_unique,
                        weights=rep(1,length(var.raw)), poly.index=1,
                        eval.grid=FALSE, display="none")$estimate
    
    i1 <- diff(range(t_unique))/4
    i2 <- 3*i1  
    id <- which(t_unique >= i1 & t_unique <= i2)
    delta <- c(0,diff(t_unique))
    
    ##V(t)-C(t,t)
    diff <- (hatV[id]-tildeG[id])*delta[id]   
    evar <- abs(2*sum(diff)/diff(range(t_unique)))
    
    ##PART 2 - EIGENVALUES AND VECTORS------------------------------
    ##Calculate Eigenvalues and vectors-----------------------
    res <- eigCov(t.red=t.red, time=t_unique, cov=covest, hmu=hmu, dimBase=dimBase, smoothFct1D=smoothFct1D)
    eigval <- res$eigval
    base  <- res$base
    if(is.null(dim(base)))
        base <- t(t(base))
    if(reg)
        fullBase <- base
    else
        fullBase <- base[na.omit(as.vector(t(indx_tin))),]
    varprop <- res$varprop
    coeffs <- coeffFct(ctrY=ctrY,
                       base=base,
                       dimBase=dimBase,
                       Delta=Delta,
                       indx_tin=indx_tin,
                       eigval=eigval,
                       average=average,
                       evar=evar)$coeffs

    plotParams <- list(t_unique=t_unique, mean.timeIn=mean.timeIn,
                       mu.smoothed=mu.smoothed, mean.raw=mean.raw,
                       cov.timeIn=cov.timeIn, cov.raw=cov.raw,
                       cov.timeOut=cov.timeOut,
                       cov.smoothed=cov.smoothed, tildeG=tildeG,
                       hatV=hatV, base=base, nt=nt)
    
    return(list(Yin=Yin, Tin=Tin, isobs=isobs, fullBase=fullBase, ctrY=ctrY,
                Delta=Delta, meanfcn=mu.smoothed, covfcn=covest,
                base=base, eigval=eigval, coeffs=coeffs,
                varprop=varprop, evar=evar, average=average, plotParams=plotParams))
}


estimateCoeffs <- function(ctrY, dimBase, base, Delta, indx_tin, eigval,
                           average, evar){
    if(dimBase>1){
        Gamma  <- diag(eigval[1:dimBase])
        base <- base[,1:dimBase]
    }else{
        Gamma <- eigval[1:dimBase]
        base <- t(t(base))
    }
    
    nt <- dim(ctrY)[2]
    if(is.null(dim(ctrY))){
        nc <- 1
        ctrY <- t(ctrY)
    }else
        nc <- dim(ctrY)[1]

    if(average){
        coeffs <- as.matrix(t(Gamma%*%t(na.omit(base))%*%solve(na.omit(base)%*%Gamma%*%t(na.omit(base))+evar*diag(nt))%*%t(ctrY)))         
    }else{
        if(dimBase==1){
            coeffs <- t(sapply(1:nc, function(x)
                t(na.omit(base[indx_tin[x,],]))%*%solve(na.omit(base[indx_tin[x,],])%*%t(na.omit(base[indx_tin[x,],]))+evar/Gamma*diag(length(na.omit(ctrY[x,]))))%*%na.omit(ctrY[x,])))
            
        }else{
            coeffs <- t(sapply(1:nc, function(x)
                Gamma%*%t(na.omit(base[indx_tin[x,],]))%*%solve(na.omit(base[indx_tin[x,],])%*%Gamma%*%t(na.omit(base[indx_tin[x,],]))+evar*diag(length(na.omit(ctrY[x,]))))%*%na.omit(ctrY[x,])))
        }
    }
    if(dimBase==1)
        coeffs <- t(coeffs)
    
    return(list(coeffs=coeffs, ctrY=ctrY))
}


integrateCoeffs <- function(ctrY, dimBase, base, Delta, indx_tin,
                            eigval=NULL, average=NULL, evar=NULL){
    if(is.null(dim(ctrY))){
        n <- 1
        ctrY <- t(ctrY)
    }
    n <- dim(ctrY)[1]
    coeffs <- matrix(0,nrow=n, ncol=dimBase)

    if(average){
        coeffs <- as.matrix((ctrY*Delta)%*%base)
        return(list(coeffs=coeffs, ctrY=ctrY))
    }else{
        for(i in 1:n){
            if(is.null(dim(base))){
                basei <- base[na.omit(indx_tin[i,])]
            }else{
                basei <-  base[na.omit(indx_tin[i,]),1:dimBase]
            }
            yi <- na.omit(ctrY[i,])*na.omit(Delta[i,])
            if(is.null(dim(yi))){
                yi <- t(yi)
            }
            coeffs[i,] <- yi%*%basei
        }
        return(list(coeffs=coeffs, ctrY=ctrY))
    }
}


avCtrY <- function(Yin, Tin, time, isobs, mean){
    ctrY <- t(t(Yin)-mean)*isobs
    return(ctrY)
}


makeGrid <- function(nt, time){
    mm <- nt*(nt+1)/2
    tempid <- matrix(1,nt,nt);
    tempid <- upper.tri(tempid, diag=TRUE)
    src <- which(t(tempid) == 1, arr.ind=TRUE)
    idr <- src[,1]
    idc <- src[,2] 
    xout <- cbind(time[idc],time[idr])
    return(xout)
}


makeFullGrid <- function(time){
    grid <- expand.grid(time, time)
    indx <- grid[,1]==grid[,2]
    res <- apply(grid,2,function(x) x[!indx])[,c(2,1)]
    return(res)
}


makeCovMat <- function(cov, nt){
    covest<-matrix(0, nt, nt)  
    count <-0
    for (i in 1:nt){
        for (j in i:nt){
            count<-count+1
            covest[i,j]<-cov[count]
            covest[j,i]<-covest[i,j]
        }
        if (covest[i,i] < 0) covest[i,i] <- 0
    }
    return(covest)
}


eigCov <- function(t.red=NULL, time, cov, hmu, dimBase, smoothFct1D){
    if(is.null(t.red)){
        t.red <- time
        quadwt <- diag(c(mean(diff(t.red)),diff(t.red)))
        tempA <- sqrt(quadwt)%*%cov%*%sqrt(quadwt)
    }else{
        tempA <- cov
    }
    eig <- eigen(tempA, symmetric=TRUE)
    ##Eigval
    eigval <- eig$values[1:dimBase]
    eigval[eigval<0] <- 0
    
    ##Varprop
    varprop <- (eigval/sum(eigval))[1:dimBase]
    
    ##Eigfcn
    vecs <- apply(eig$vectors, 2, function(x)
        smoothFct1D(x=t.red, y=x, h=hmu,
                    eval.points=time,
                    poly.index=2,
                    weights=rep(1,length(t.red)),
                    display="none")$estimate)
    
    if(sum(is.na(vecs))>0)
        stop("Please choose a larger number of maximal time points in fpcCtrl (not enought to smooth).")    
    temp <- apply(vecs, 2, function(x)
        1/sqrt(trapz(time,x^2)))
    invW <- diag(temp)
    base <- (vecs%*%invW)[,1:dimBase]
    
    return(list(base=base, eigval=eigval, varprop=varprop))
}


avCov <- function(ctrY, Tin, isobs, time){
    nt <- length(time)
    Wsum <- as.matrix(t(isobs)%*%isobs)
    cov.raw <- as.matrix(t(ctrY)%*%ctrY/Wsum)
    id <- which(diag(1,nt)==1)
    var.raw <- cov.raw[id]
    cov.raw <- cov.raw[-id]
    isobs <- Wsum[-id]!=0
    cov.raw <- cov.raw[isobs==1]
    cov.timeIn <- makeFullGrid(time)
    cov.timeIn <- cov.timeIn[isobs==1,]
    var.timeIn <- time
    
    return(list(cov.timeIn=cov.timeIn, cov.raw=cov.raw, var.raw=var.raw, var.timeIn=var.timeIn))
}


longCov <- function(ctrY, Tin, isobs, time){
    N <- rowSums(isobs)
    cov.raw <-numeric(sum(N*(N-1)/2))     
    cov.timeIn <- matrix(0,nrow=sum(N*(N-1)/2),ncol=2)
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


avMean <- function(Yin, Tin, isobs){
    time <- sort(unique(Tin[which(isobs==1,
                                  arr.ind=T)]))
    if(is.null(dim(Yin))){
        Yin <- (t(Yin)); Tin <- (t(Tin)); isobs <- (t(isobs))
    }
    obsPerTime <- colSums(isobs)
    mean <- colSums(Yin)[obsPerTime!=0]/obsPerTime[obsPerTime!=0]
    return(list(mu=mean, time=time))
}


longMean <- function(Yin, Tin, isobs){
    isobsVec <- as.vector(t(isobs))
    yVec <- as.vector(t(Yin))
    mean <- yVec[isobsVec==1] 
    timeVec <- as.vector(t(Tin))
    time <- timeVec[isobsVec==1]
    return(list(mu=mean, time=time))
}


