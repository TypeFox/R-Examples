
hiddenICP <-
function(X,Y,ExpInd,alpha=0.1,mode="asymptotic",intercept=FALSE){
    if(is.vector(X) & is.numeric(X)) X <- matrix(X,ncol=1)
    if(!is.matrix(X) & !is.data.frame(X)) stop("'X' must be a matrix or data frame")
    if(!is.vector(Y)) stop("'Y' must be a vector")
    if(nrow(X) <= ncol(X)) stop(paste( "hiddenICP not suitable for high-dimensional data (at the moment) \n -- need nrow(X)>ncol(X) but have nrow(X)=",nrow(X)," and ncol(X)=",ncol(X)))

    
    if(!is.list(ExpInd)){
        if(length(ExpInd)!=length(Y)) stop("if `ExpInd' is a vector, it needs to have the same length as `Y'")
        uni <- unique(ExpInd)
        if(length(uni)==1) stop(paste("there is just one environment ('ExpInd'=",uni[1]," for all observations) and the method needs at least two distinct environments",sep=""))
        if(min(table(ExpInd))<=2){
            cat("\n out put of 'table(ExpInd)':\n ")
            print(table(ExpInd))
            stop("one environment has just one or two observations (as supllied by 'ExpInd'); there need to be at least 3 (and ideally dozens) of observations in each environment; the out put of 'table(ExpInd)' is given below to show the number of observations in each unique environment as supplied by 'ExpInd'")
        }
        
        K <- length(uni)
        ExpIndNEW <- list()
        for (uc in 1:K){
            ExpIndNEW[[uc]] <- which(ExpInd==uni[uc])
            attr(ExpIndNEW[[uc]],"value") <- uni[uc]
        }
        ExpInd <- ExpIndNEW
        rm(ExpIndNEW)
    }else{
        ran <- range(unlist(ExpInd))
        if(ran[1]<1) stop(paste("if `ExpInd' is a list with indicies of observations, \n minimal entry has to be at least 1 but is",ran[1]))
        if(ran[2]>length(Y)) stop(paste("if `ExpInd' is a list with indicies of observations, \n maximal entry has to be at most equal \n to the length",length(Y),"of the observations but is",ran[2]))
    }
    X <- as.matrix(X)
    if(length(ucol <- unique(colnames(X)))< min(3, ncol(X))) colnames(X) <- paste("Variable",1:ncol(X),sep="_")
    colX <- colnames(X)

    
    if(intercept){
        X <- cbind(rep(1,nrow(X)),X)
    }

    K <- length(ExpInd)
    p <- ncol(X)
    n <- nrow(X)

    kc <- 1
    KC <- if(K>2) K else 1
    ConfInt <- matrix(NA, nrow=2,ncol=p)
    pvalues <- rep(1,p)
    for (kc in 1: KC){
        ins <- ExpInd[[kc]]
        out <- (1:n)[-ins]
        DS <- (t(X[ins,])%*%X[ins,])/length(ins) - (t(X[out,])%*%X[out,])/length(out)
        Drho <- (t(X[ins,])%*%Y[ins])/length(ins) - (t(X[out,])%*%Y[out])/length(out)
        DSI <- solve(DS)
    
        
        betahat <- as.numeric(solve( DS,Drho))
        if(kc==1) betahatall <- betahat else betahatall <- betahatall+betahat
        Zin <- matrix(nrow=length(ins),ncol=p)
        Zout <- matrix(nrow=length(out),ncol=p)
        for (i in 1:length(ins)){
            tmp <- DSI %*% X[ins[i],]
            Zin[i,] <- as.numeric(- tmp * sum(tmp*Drho) + Y[ins[i]] *tmp)
        }
        for (i in 1:length(out)){
            tmp <- DSI %*% X[out[i],]
            Zout[i,] <- as.numeric(- tmp * sum(tmp*Drho) + Y[out[i]] *tmp)
        }
        sigmavec <- sqrt(diag((cov(Zin)/length(ins)+cov(Zout)/length(out))))

        pvalues <- pmin(pvalues, 2*K* (1-pt( abs(betahat)/pmax(10^(-10),sigmavec),df=n-1)))
        
        addvar <- qnorm(max(0.5,1-alpha/(2*K))) * sigmavec
        maximineffectsN <- sign(betahat) * pmax( 0, abs(betahat) - addvar)
        ConfInt[1,] <- pmax(ConfInt[1,], betahat - addvar,na.rm=TRUE) 
        ConfInt[2,] <- pmin(ConfInt[2,], betahat + addvar,na.rm=TRUE) 
        if(kc==1){
            maximineffects <- maximineffectsN
        }else{
            for (varc in 1:p){
                if(abs(maximineffectsN[varc]) > abs(maximineffects[varc])) maximineffects[varc] <- maximineffectsN[varc]
            }
        }
    }
    betahat <- betahatall/KC
    maximinCoefficients <- maximineffects
    if(intercept){
        betahat <- betahat[-1]
        maximinCoefficients <- maximinCoefficients[-1]
        ConfInt <- ConfInt[,-1]
        pvalues <- pvalues[-1]
    }
    ConfInt <- apply(ConfInt,2,sort,decreasing=TRUE)
    retobj <- list(betahat=betahat,maximinCoefficients=maximinCoefficients,ConfInt=ConfInt,pvalues=pvalues,colnames=colX,alpha=alpha)
    class(retobj) <- "hiddenInvariantCausalPrediction"
    return(retobj)
                       
}



