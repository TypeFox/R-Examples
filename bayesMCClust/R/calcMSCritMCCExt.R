calcMSCritMCCExt <-
function(workDir, NN, myLabel="model choice for ...", ISdraws=3, H0=3, whatToDoList=c("approxMCL", "approxML", "postMode") ) {

prevwd <- getwd()

setwd(workDir)

outFileNamesMCC <- list.files(pattern="[0123456789].RData")    

xval <- numeric(length(outFileNamesMCC)) 

Hmax <- Hnr <- length(xval) 

N <- logLike <- logXiPrior <- Prior <- logBetaPrior <- logClassLike <- xi.m <- K <- Beta.m <- Data <- Njk.i <- Mcmc <- entropy <- c0 <- Initial <- NULL

toDoIsPostMode <- toDoIsApproxML <- toDoIsApproxMCL <- FALSE

MSCritList <- list()

for ( whatToDo in whatToDoList ) { 

# MCC results

print("MCC Logit")
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
classMCC <- matrix(0, NN, Hmax)
corrClassMCC <- numeric(Hmax)
mMaxMCC <- numeric(Hmax)

if ( toDoIsPostMode ) maxLogPostDens <- numeric(Hmax) 
if ( toDoIsApproxML ) maxLogLike <- numeric(Hmax)
if ( toDoIsApproxMCL )maxLogClassLike <- numeric(Hmax)

for ( HH in 1:Hnr ) { 
    
    results <- load( outFileNamesMCC[HH] )
        
    if ( exists("N") )  { if ( dim(classMCC)[1] != N ) next } else { break }
        
    logPostDens <- logLike + logXiPrior + if ( Prior$betaPrior!="uninformative" ) logBetaPrior else 0 

    if ( toDoIsPostMode )  { mMax <- which.max(logPostDens);  maxLogPostDens[HH] <- max(logPostDens) }
    if ( toDoIsApproxML )  { mMax <- which.max(logLike);      maxLogLike[HH] <- max(logLike) }
    if ( toDoIsApproxMCL ) { mMax <- which.max(logClassLike); maxLogClassLike[HH] <- max(logClassLike) }
    
    mMaxMCC[HH] <- mMax
    
    xiPostMode <- array(xi.m[mMax,,,], c(K+1, K+1, Prior$H)) 
    betaPostMode <- Beta.m[,,mMax] 
        
    XBetak <- crossprod(t(Data$X), apply( betaPostMode, c(1,2), mean) )
    logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) 
    
    # observed log likelihood
    tempLK <-logit.temp * matrix(mapply(function(i,h) prod(xiPostMode[,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE )
    LK <- sum(log(rowSums(tempLK)))

    dh <- (K + 1)*K * Prior$H + ncol(Data$X) * (Prior$H-1)

    BicMCC[HH] <- bick <- -2*LK + dh*log(N)
    
    adjBicMCC[HH] <- adjbick <- bick - 2*Prior$H*log(H0) + 2*lgamma(Prior$H + 1) + 2*H0

    AicMCC[HH] <- aick <- -2*LK + 2*dh
    
    sProbs <- tempLK
    sProbsNorm  <-  sProbs/rowSums(sProbs)    
    classMCC[,HH] <- class <- max.col(sProbsNorm) 
    
    # classification log likelihood
    CLK <- sum( log( ( matrix( mapply( function(i,h) 
       prod( xiPostMode[,,h]^Njk.i[,,i] ), rep(1:N, each=Prior$H), 1:Prior$H ), N,Prior$H, byrow=TRUE))[cbind(1:N,class)]*logit.temp[cbind(1:N,class)])) 

    AweMCC[HH] <- awek <- -2*CLK + 2*dh*( 3/2 + log(N) )
    
    tik <- sProbsNorm
    
    stopifnot( min(tik) > 1e-320 )
    
    EK <- -sum(tik*log(tik))
    
    IclBicMCC[HH] <- Iclbic <- -2*LK +  log(N)*dh + 2*EK
    
    ClcMCC[HH] <- clc <- -2*LK + 2*EK
    
    Dic2MCC[HH] <- dic2 <- -4*mean(logLike[Mcmc$M0:Mcmc$M]) + 2* LK
    
    Dic4MCC[HH] <- dic4 <- dic2 + 2*mean(entropy[Mcmc$M0:Mcmc$M])
    
    xval[HH] <- Prior$H
    
    e0jkh <- c0
    
    Njkh <- array(0, dim(e0jkh))
    for ( h in 1:Prior$H ) { Njkh[,,h] <- apply( Njk.i[,,class==h] , c(1, 2), sum) }
    
    betaMeanPropos <- apply(Beta.m[,,Mcmc$M0:Mcmc$M], c(1,2), mean) # dim = . x H
    betaCovarPropos <- array(0, c(dim(betaMeanPropos)[1], dim(betaMeanPropos)[1], Prior$H) ) # dim = . x . x H
    for ( h in 1:Prior$H ) {
        betaCovarPropos[,,h] <- cov( t( Beta.m[,h,Mcmc$M0:Mcmc$M] ) )    
    }
    
    pS_temp <- numeric(ISdraws)
    
    for ( im in 1:length(pS_temp) ) {
        # propose betas:
        betaProposed <- matrix(0, dim(betaMeanPropos)[1], dim(betaMeanPropos)[2] ) # dim = . x H
        for ( h in 2:Prior$H ) {
            betaProposed[, h] <- rmnorm(n = 1, mean = betaMeanPropos[,h], varcov = betaCovarPropos[,,h])         
        }
        XBetak <- crossprod(t(Data$X), betaProposed) # dim = N x Prior$H
        logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) # dim = N x Prior$H
        logMNLike <- sum( log( logit.temp[ cbind( 1:N, class ) ] ) ) 
        logPriorDens <- numeric(Prior$H)
        for ( h in 2:Prior$H ) logPriorDens[h] <- dmnorm( betaProposed[, h], mean = rep( Prior$betaPriorMean, dim(betaMeanPropos)[1]), 
                                                   varcov = diag( rep( Prior$betaPriorVar, dim(betaMeanPropos)[1] ) ), log = TRUE )
        logPrior <- sum( logPriorDens[2:Prior$H] )
        logProposDens <- numeric(Prior$H)
        for ( h in 2:Prior$H ) {
            logProposDens[h] <- dmnorm( betaProposed[, h], mean = betaMeanPropos[,h], varcov = betaCovarPropos[,,h], log = TRUE )
        }
        logPropos <- sum( logProposDens[2:Prior$H] )
        pS_temp[im] <- logMNLike + logPrior - logPropos
    }
    pS <- log( mean( exp( pS_temp - max( pS_temp ) ) ) ) + max( pS_temp ) 
    
    logICLMCC[HH] <- logicl <- -2*(
    sum(
    colSums( lgamma ( apply(e0jkh, c(3), colSums) ) )
    - colSums( apply( lgamma(e0jkh), c(3), colSums) )
    + colSums( apply( lgamma(e0jkh + Njkh), c(3), colSums) )
    - colSums( lgamma ( apply(e0jkh + Njkh, c(3), colSums) ) )
    ) + pS )

    corrClassMCC[HH]<-corrClass<-if ( Prior$H > 1 ) if (is.null(Initial$S.i.start)) NA else suppressWarnings(round(e1071::classAgreement(table(Initial$S.i.start, class))$crand*100,5))     
    
    print(paste(HH, ". MCC Logit: H =", Prior$H, "; BIC =", round(bick,2), "; adjBIC =", round(adjbick,2), "; AIC =", round(aick,2), 
             "; AWE =", round(awek,2), "; CLC =", round(clc,2), "; IclBic =", round(Iclbic,2), 
             "; DIC2 =", round(dic2,2), "; DIC4a =", round(dic4,2), "; ICL =", round(logicl,2),
             "; adj Rand:", corrClass, " " ))
                
    flush.console()
    
    rm( list = setdiff( ls(all=TRUE), c("outFileNamesMCC", "Hmax", "whatToDo", "BicMCC", "adjBicMCC", "AicMCC", "AweMCC", "classMCC", "NN",
                                        "IclBicMCC", "ClcMCC", "Dic2MCC", "Dic4MCC", "corrClassMCC", "mMaxMCC", "HH", "Hnr", "xval", "myLabel", "H0",
                                        "ISdraws", "logICLMCC", "MSCritList", "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd",
                                        if ( toDoIsPostMode ) "maxLogPostDens",
                                        if ( toDoIsApproxML ) "maxLogLike",
                                        if ( toDoIsApproxMCL ) "maxLogClassLike" 
                                        ) ) )  
        
} # end for HH

indi <- 1:Hnr

postscript( paste("MCCLogitModSelCrits_",whatToDo,".eps",sep=""), width = 10, height = 7)
par(mfrow=c(2,4), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1)) # c(bottom, left, top, right)
plot(xval, BicMCC[indi], main="BIC", xlab="H", ylab="")
plot(xval, adjBicMCC[indi], main="adjBIC", xlab="H", ylab="")
plot(xval, Dic2MCC[indi], main="DIC2", xlab="H", ylab="")
plot(xval, AweMCC[indi], main="AWE", xlab="H", ylab="")
plot(xval, ClcMCC[indi], main="CLC", xlab="H", ylab="")
plot(xval, IclBicMCC[indi], main="ICL-BIC", xlab="H", ylab="")
plot(xval, logICLMCC[indi], main="ICL", xlab="H", ylab="")
plot(xval, Dic4MCC[indi], main="DIC4a", xlab="H", ylab="")
mtext(paste(myLabel, "MCC Logit", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
dev.off()

pdf( paste("MCCLogitModSelCrits_",whatToDo,".pdf",sep=""), width = 10, height = 7)
par(mfrow=c(2,4), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1)) # c(bottom, left, top, right)
plot(xval, BicMCC[indi], main="BIC", xlab="H", ylab="")
plot(xval, adjBicMCC[indi], main="adjBIC", xlab="H", ylab="")
plot(xval, Dic2MCC[indi], main="DIC2", xlab="H", ylab="")
plot(xval, AweMCC[indi], main="AWE", xlab="H", ylab="")
plot(xval, ClcMCC[indi], main="CLC", xlab="H", ylab="")
plot(xval, IclBicMCC[indi], main="ICL-BIC", xlab="H", ylab="")
plot(xval, logICLMCC[indi], main="ICL", xlab="H", ylab="")
plot(xval, Dic4MCC[indi], main="DIC4a", xlab="H", ylab="")
mtext(paste(myLabel, "MCC Logit", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
dev.off()

print(getwd())
print(whatToDo)

MSCritTable <- cbind(H=xval, mMax=mMaxMCC[indi], 
          if ( toDoIsPostMode ) maxLPD=maxLogPostDens[indi] , 
          if ( toDoIsApproxML ) maxLL=maxLogLike[indi] , 
          if ( toDoIsApproxMCL ) maxLCL=maxLogClassLike[indi] ,
          BIC=BicMCC[indi], adjBIC=adjBicMCC[indi], AIC=AicMCC[indi], AWE=AweMCC[indi], CLC=ClcMCC[indi], IclBic=IclBicMCC[indi], 
          DIC2=Dic2MCC[indi], DIC4a=Dic4MCC[indi], logICL=logICLMCC[indi], adjRand=corrClassMCC[indi])

print(MSCritTable)
cat("\n\n")

if ( toDoIsPostMode  ) MSCritList$postMode  <- MSCritTable
if ( toDoIsApproxML  ) MSCritList$approxML  <- MSCritTable
if ( toDoIsApproxMCL ) MSCritList$approxMCL <- MSCritTable

#save.image( paste("MCCLogitModSelCrits_",whatToDo,".RData",sep="") )
save(list = c("outFileNamesMCC", "Hmax", "whatToDo", "BicMCC", "adjBicMCC", "AicMCC", "AweMCC", "classMCC", "NN",
              "IclBicMCC", "ClcMCC", "Dic2MCC", "Dic4MCC", "corrClassMCC", "mMaxMCC", "HH", "Hnr", "xval", "myLabel", "H0",
              "ISdraws", "logICLMCC", "MSCritList", "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd", "indi",
              if ( toDoIsPostMode ) "maxLogPostDens",
              if ( toDoIsApproxML ) "maxLogLike",
              if ( toDoIsApproxMCL ) "maxLogClassLike" 
              ), file = paste("MCCLogitModSelCrits_",whatToDo,".RData",sep="") )

} # end for whatToDo

cat("No of ISdraws =", ISdraws, "\n\n")

MSCritList$ISdraws <- ISdraws

print(cbind(outFileNames=outFileNamesMCC), quote=FALSE)
    
MSCritList$outFileNames <- outFileNamesMCC

setwd(prevwd)

return( invisible( MSCritList ) )

}
