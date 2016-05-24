calcMSCritMCC <-
function(workDir, myLabel="model choice for ...", H0=3, whatToDoList=c("approxMCL", "approxML", "postMode") ) {

prevwd <- getwd()

setwd(workDir)

outFileNamesMCC <- list.files(pattern="[0123456789].RData" )   

xval <- numeric(length(outFileNamesMCC)) 

Hmax <- Hnr <- length(xval) 

logLike <- logXiPrior <- logEtaPrior <- logClassLike <- xi.m <- K <- Prior <- eta.m <- N <- Njk.i <- Mcmc <- entropy <- c0 <- Initial <- NULL

toDoIsPostMode <- toDoIsApproxML <- toDoIsApproxMCL <- FALSE

MSCritList <- list()

for ( whatToDo in whatToDoList ) {

    # MCC results

    print("MCC")
    print(whatToDo)
    
    toDoIsPostMode  <- if ( whatToDo=="postMode" )  TRUE else FALSE
    toDoIsApproxML  <- if ( whatToDo=="approxML" )  TRUE else FALSE
    toDoIsApproxMCL <- if ( whatToDo=="approxMCL" ) TRUE else FALSE

    stopifnot( setequal( c(toDoIsPostMode, toDoIsApproxML, toDoIsApproxMCL), c(FALSE, TRUE, FALSE) ) )

    BicMCC <- numeric(Hmax)
    adjBicMCC <- numeric(Hmax)
    AicMCC <- numeric(Hmax)
    AweMCC <- numeric(Hmax)
    IclBicMCC <- numeric(Hmax)
    ClcMCC <- numeric(Hmax)
    Dic2MCC <- numeric(Hmax)
    Dic4MCC <- numeric(Hmax)
    logICLMCC <- numeric(Hmax)
    corrClassMCC <- numeric(Hmax)
    mMaxMCC <- numeric(Hmax)
    
    if ( toDoIsPostMode )  maxLogPostDens <- numeric(Hmax) 
    if ( toDoIsApproxML )  maxLogLike <- numeric(Hmax)
    if ( toDoIsApproxMCL ) maxLogClassLike <- numeric(Hmax)

    for ( HH in 1:Hnr ) { 
        results <- load( outFileNamesMCC[HH] )
    
        logPostDens <- logLike + logXiPrior + logEtaPrior

        if ( toDoIsPostMode )  { mMax <- which.max(logPostDens);  maxLogPostDens[HH] <- max(logPostDens) }
        if ( toDoIsApproxML )  { mMax <- which.max(logLike);      maxLogLike[HH] <- max(logLike) }
        if ( toDoIsApproxMCL ) { mMax <- which.max(logClassLike); maxLogClassLike[HH] <- max(logClassLike) }
    
        mMaxMCC[HH] <- mMax
    
        xiPostMode <- array(xi.m[mMax,,,], c(K+1, K+1, Prior$H))
        etaPostMode <- eta.m[ , mMax ]
    
        estGroupSizeMat <- matrix(etaPostMode, N, Prior$H, byrow=TRUE)
        tempLK <- 
            estGroupSizeMat * matrix(mapply(function(i,h) prod(xiPostMode[,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE )
        LK <- sum(log(rowSums(tempLK)))
        
        Nobs <- N  

        dh <- (K + 1)*K * Prior$H + (Prior$H - 1) 
        BicMCC[HH] <- bick <- -2*LK + dh*log(Nobs)
    
        adjBicMCC[HH] <- adjbick <- bick - 2*Prior$H*log(H0) + 2*lgamma(Prior$H + 1) + 2*H0

        AicMCC[HH] <- aick <- -2*LK + 2*dh
    
        sProbs <- tempLK
        sProbsNorm  <-  sProbs/rowSums(sProbs)    
        class <- max.col(sProbsNorm)  
    
        # classification log likelihood
        CLK <- sum( log( ( matrix( mapply(function(i,h) 
            prod(xiPostMode[,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE ) )[cbind(1:N,class)]*etaPostMode[class]))

        AweMCC[HH] <- awek <- -2*CLK + 2*dh*( 3/2 + log(Nobs) )
    
        tik <- sProbsNorm
    
        stopifnot( min(tik) > 1e-320 )
    
        EK <- -sum(tik*log(tik))
    
        IclBicMCC[HH] <- Iclbic <- -2*LK +  log(Nobs)*dh + 2*EK
    
        ClcMCC[HH] <- clc <- -2*LK + 2*EK
    
        Dic2MCC[HH] <- dic2 <- -4*mean(logLike[Mcmc$M0:Mcmc$M]) + 2* LK
    
        Dic4MCC[HH] <- dic4 <- dic2 + 2*mean(entropy[Mcmc$M0:Mcmc$M])
    
        xval[HH] <- Prior$H
        
        N_h <- table( class )
        e0jkh <- c0
        alpha0h <- rep( Prior$e0, Prior$H )
        Njkh <- array(0, dim(e0jkh))
        for ( h in 1:Prior$H ) { Njkh[,,h] <- apply( Njk.i[,,class==h] , c(1, 2), sum) }
        logICLMCC[HH] <- logicl <- -2*(
        sum(
        colSums( lgamma ( apply(e0jkh, c(3), colSums) ) )
         - colSums( apply( lgamma(e0jkh), c(3), colSums) )
         + colSums( apply( lgamma(e0jkh + Njkh), c(3), colSums) )
         - colSums( lgamma ( apply(e0jkh + Njkh, c(3), colSums) ) )
        ) + lgamma( sum( alpha0h ) ) - sum ( lgamma( alpha0h ) ) + sum ( lgamma( alpha0h + N_h ) ) - lgamma( sum( alpha0h + N_h ) ) )
    
        corrClassMCC[HH]<-corrClass<-if ( Prior$H > 1 ) if (is.null(Initial$S.i.start)) NA else suppressWarnings(round(e1071::classAgreement(table(Initial$S.i.start, class))$crand*100,5))     
     
        print(paste(HH, ". MCC: H =", Prior$H, "; BIC =", round(bick,2), "; adjBIC =", round(adjbick,2), "; AIC =", round(aick,2), 
             "; AWE =", round(awek,2), "; CLC =", round(clc,2), "; IclBic =", round(Iclbic,2), 
             "; DIC2 =", round(dic2,2), "; DIC4a =", round(dic4,2), "; ICL =", round(logicl,2), 
             "; adj Rand:", corrClass, "%" ))

        flush.console()

        rm( list = setdiff( ls(all=TRUE), c("outFileNamesMCC", "Hmax", "whatToDo", "BicMCC", "adjBicMCC", "AicMCC", "AweMCC", "H0",
                                        "IclBicMCC", "ClcMCC", "Dic2MCC", "Dic4MCC", "corrClassMCC", "mMaxMCC", "HH", "Hnr", "xval", "myLabel", "logICLMCC", "MSCritList",  
                                        "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd",                                        
                                        if ( toDoIsPostMode ) "maxLogPostDens",
                                        if ( toDoIsApproxML ) "maxLogLike",
                                        if ( toDoIsApproxMCL ) "maxLogClassLike"                                        
                                        ) ) )  
        
    }   # end for HH

    indi <- 1:Hnr
 
    postscript( paste("MCCModSelCrits_",whatToDo,".eps",sep=""), width = 10, height = 7)
    par(mfrow=c(2,4), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1)) # c(bottom, left, top, right)
    plot(xval, BicMCC[indi], main="BIC", xlab="H", ylab="")
    plot(xval, adjBicMCC[indi], main="adjBIC", xlab="H", ylab="")
    plot(xval, Dic2MCC[indi], main="DIC2", xlab="H", ylab="")
    plot(xval, AweMCC[indi], main="AWE", xlab="H", ylab="")
    plot(xval, ClcMCC[indi], main="CLC", xlab="H", ylab="")
    plot(xval, IclBicMCC[indi], main="ICL-BIC", xlab="H", ylab="")
    plot(xval, logICLMCC[indi], main="ICL", xlab="H", ylab="")
    plot(xval, Dic4MCC[indi], main="DIC4a", xlab="H", ylab="")
    mtext(paste(myLabel, "MCC", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
    dev.off()

    pdf( paste("MCCModSelCrits_",whatToDo,".pdf",sep=""), width = 10, height = 7)
    par(mfrow=c(2,4), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1)) # c(bottom, left, top, right)
    plot(xval, BicMCC[indi], main="BIC", xlab="H", ylab="")
    plot(xval, adjBicMCC[indi], main="adjBIC", xlab="H", ylab="")
    plot(xval, Dic2MCC[indi], main="DIC2", xlab="H", ylab="")
    plot(xval, AweMCC[indi], main="AWE", xlab="H", ylab="")
    plot(xval, ClcMCC[indi], main="CLC", xlab="H", ylab="")
    plot(xval, IclBicMCC[indi], main="ICL-BIC", xlab="H", ylab="")
    plot(xval, logICLMCC[indi], main="ICL", xlab="H", ylab="")
    plot(xval, Dic4MCC[indi], main="DIC4a", xlab="H", ylab="")
    mtext(paste(myLabel, "MCC", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
    dev.off()

    print(getwd())
    print(whatToDo)
    
    MSCritTable <- cbind(H=xval, mMax=mMaxMCC[indi], 
          if ( toDoIsPostMode ) maxLPD=maxLogPostDens[indi] , 
          if ( toDoIsApproxML ) maxLL=maxLogLike[indi] , 
          if ( toDoIsApproxMCL ) maxLCL=maxLogClassLike[indi] ,
          BIC=BicMCC[indi], adjBIC=adjBicMCC[indi], AIC=AicMCC[indi], AWE=AweMCC[indi], CLC=ClcMCC[indi], 
          IclBic=IclBicMCC[indi], DIC2=Dic2MCC[indi], DIC4a=Dic4MCC[indi], logICL=logICLMCC[indi], adjRand=corrClassMCC[indi])
    
    print(MSCritTable)
        
    cat("\n\n")
        
    if ( toDoIsPostMode  ) MSCritList$postMode  <- MSCritTable
    if ( toDoIsApproxML  ) MSCritList$approxML  <- MSCritTable
    if ( toDoIsApproxMCL ) MSCritList$approxMCL <- MSCritTable
    
    # save.image( paste("MCCModSelCrits_",whatToDo,".RData",sep="") )
    save(list = c("outFileNamesMCC", "Hmax", "whatToDo", "BicMCC", "adjBicMCC", "AicMCC", "AweMCC", "H0",
                                        "IclBicMCC", "ClcMCC", "Dic2MCC", "Dic4MCC", "corrClassMCC", "mMaxMCC", "HH", "Hnr", "xval", "myLabel", "logICLMCC", "MSCritList",  
                                        "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd", "indi",                                       
                                        if ( toDoIsPostMode ) "maxLogPostDens",
                                        if ( toDoIsApproxML ) "maxLogLike",
                                        if ( toDoIsApproxMCL ) "maxLogClassLike"                                        
                                        ), file = paste("MCCModSelCrits_",whatToDo,".RData",sep="") )


} # end for whatToDo

print(cbind(outFileNames=outFileNamesMCC), quote=FALSE)
    
MSCritList$outFileNames <- outFileNamesMCC

setwd(prevwd)

return( invisible( MSCritList ) )

}
