calcMSCritDMC <-
function(workDir, myLabel="model choice for ...", myN0="N0 = ...", whatToDoList=c("approxMCL", "approxML", "postMode") ) {

prevwd <- getwd()

setwd(workDir)

outFileNamesDMC <-  list.files(pattern="[0123456789].RData" )   

xval <- numeric(length(outFileNamesDMC)) 

Hmax <- Hnr <- length(xval)

logLike <- logEtaPrior <- logEPrior <- logClassLike <- e_h_m <- K <- Prior <- eta_m <- N <- Prior <- Njk.i <- Initial <- NULL

toDoIsPostMode <- toDoIsApproxML <- toDoIsApproxMCL <- FALSE

MSCritList <- list()

for ( whatToDo in whatToDoList ) {

# DMC results

print("DMC")
print(whatToDo)

toDoIsPostMode  <- if ( whatToDo=="postMode" )  TRUE else FALSE
toDoIsApproxML  <- if ( whatToDo=="approxML" )  TRUE else FALSE
toDoIsApproxMCL <- if ( whatToDo=="approxMCL" ) TRUE else FALSE

stopifnot( setequal( c(toDoIsPostMode, toDoIsApproxML, toDoIsApproxMCL), c(FALSE, TRUE, FALSE) ) )

margDataLL_0 <- function(i, h, Njk.i=Njk.i, e_h=e_h, AA=AA, BB=BB) {
    CC <- sum( lgamma( Ne <- (Njk.i[,,i] + e_h[,,h]) ) )
    DD <- sum( lgamma( rowSums( Ne ) ) )
    return( exp(AA[h] - BB[h] + CC - DD) )
}


BicDMC <- numeric(Hmax)
AicDMC <- numeric(Hmax)
AweDMC <- numeric(Hmax)
IclBicDMC <- numeric(Hmax)
ClcDMC <- numeric(Hmax)
corrClassDMC <- numeric(Hmax)
mMaxDMC <- numeric(Hmax)

    if ( toDoIsPostMode ) maxLogPostDens <- numeric(Hmax) 
    if ( toDoIsApproxML ) maxLogLike <- numeric(Hmax)
    if ( toDoIsApproxMCL )maxLogClassLike <- numeric(Hmax)

for ( HH in 1:Hnr ) { 
    results <- load( outFileNamesDMC[HH] )
    
    logPostDens <- logLike + logEtaPrior + logEPrior
    
    if ( toDoIsPostMode )  { mMax <- which.max(logPostDens);  maxLogPostDens[HH] <- max(logPostDens) }
    if ( toDoIsApproxML )  { mMax <- which.max(logLike);      maxLogLike[HH] <- max(logLike) }
    if ( toDoIsApproxMCL ) { mMax <- which.max(logClassLike); maxLogClassLike[HH] <- max(logClassLike) }
   
    mMaxDMC[HH] <- mMax
       
    ePostMode <- array(e_h_m[,,,mMax], c(K+1, K+1, Prior$H)) 
    etaPostMode <- eta_m[ mMax ,  ]

    AA <- apply( lgamma(apply(ePostMode, 3, rowSums)), 2, sum )   
    BB <- apply( lgamma(ePostMode), 3, sum)   
    LK <- sum(sapply(1:N, function(i) log(sum( sapply(1:Prior$H, function(h) 
                                            margDataLL_0(i=i, h=h, Njk.i=Njk.i, e_h=ePostMode, AA=AA, BB=BB))*etaPostMode ))))
    Nobs <- sum(Njk.i) + N 

    dh <- (K + 1)^2 * Prior$H + (Prior$H - 1)  
    BicDMC[HH] <- bick <- -2*LK + dh*log(Nobs)    
    
    AicDMC[HH] <- aick <- -2*LK + 2*dh


    margDataLikeli <- matrix(0, N, Prior$H)
    for ( i in 1:N) {
        for ( h in 1:Prior$H ) {  
            margDataLikeli[i, h] <- margDataLL_0(i=i, h=h, Njk.i=Njk.i, e_h=ePostMode, AA=AA, BB=BB)
        }
    }  
    sProbs  <- margDataLikeli*matrix(etaPostMode, N, Prior$H, byrow=TRUE) 

    sProbsNorm  <-  sProbs/rowSums(sProbs)    

    class <- max.col(sProbsNorm)  

    CLK <- sum(sapply(1:N, function(i) log( ( sapply(1:Prior$H, function(h) 
                                            margDataLL_0(i=i, h=h, Njk.i=Njk.i, e_h=ePostMode, AA=AA, BB=BB))*etaPostMode )[class[i]]  )))

    AweDMC[HH] <- awek <- -2*CLK + 2*dh*( 3/2 + log(Nobs) )
    
    tik <- sProbsNorm
    
    stopifnot( min(tik) > 1e-320 )
    
    EK <- -sum(tik*log(tik))
    
    IclBicDMC[HH] <- Iclbic <- -2*LK +  log(Nobs)*dh + 2*EK
    
    ClcDMC[HH] <- clc <- -2*LK + 2*EK
    
    xval[HH] <- Prior$H
    
    corrClassDMC[HH]<-corrClass <- if ( Prior$H > 1 ) if (is.null(Initial$S.i.start)) NA else suppressWarnings(round(e1071::classAgreement(table(Initial$S.i.start, class))$crand*100,5))     
    
    print(paste(HH, ". DMC: H =", Prior$H, "; BIC =", round(bick,2), "; AIC =", round(aick,2), "; AWE =", round(awek,2), 
                "; CLC =", round(clc,2), "; IclBic =", round(Iclbic,2), "; adj Rand:", corrClass, "%" ))
    flush.console()
    
    rm( list = setdiff( ls(all=TRUE), c("outFileNamesDMC", "Hmax", "whatToDo", "margDataLL_0", "BicDMC", "AicDMC", "AweDMC", 
                                        "IclBicDMC", "ClcDMC", "corrClassDMC", "mMaxDMC", "HH", "Hnr", "xval", "myLabel", "MSCritList",
                                        "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd",
                                        if ( toDoIsPostMode ) "maxLogPostDens",
                                        if ( toDoIsApproxML ) "maxLogLike",
                                        if ( toDoIsApproxMCL ) "maxLogClassLike"
                                        ) ) )    
    
} # end for HH

indi <- 1:Hnr

postscript( paste("DMCModSelCrits_",whatToDo,".eps",sep=""), width = 10, height = 7)
par(mfrow=c(2,3), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1))
plot(xval, AicDMC[indi], main="AIC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, BicDMC[indi], main="BIC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, AweDMC[indi], main="AWE", xlab="H", ylab="", bg="black", pch=23)
plot(xval, ClcDMC[indi], main="CLC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, IclBicDMC[indi], main="ICL-BIC", xlab="H", ylab="", bg="black", pch=23)
mtext(paste(myLabel, "DMC", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
dev.off()

pdf( paste("DMCModSelCrits_",whatToDo,".pdf",sep=""), width = 10, height = 7)
par(mfrow=c(2,3), omi=c(0.0,0.0,0.25,0.0), mai=c(0.6, 0.3, 0.5, 0.1))
plot(xval, AicDMC[indi], main="AIC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, BicDMC[indi], main="BIC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, AweDMC[indi], main="AWE", xlab="H", ylab="", bg="black", pch=23)
plot(xval, ClcDMC[indi], main="CLC", xlab="H", ylab="", bg="black", pch=23)
plot(xval, IclBicDMC[indi], main="ICL-BIC", xlab="H", ylab="", bg="black", pch=23)
mtext(paste(myLabel, "DMC", whatToDo, format(Sys.time(), "%d %b %Y"), sep="   ***   "), outer=TRUE, cex=1.2)
dev.off()

print(getwd())
print(whatToDo)

MSCritTable <- cbind(H=xval, mMax=mMaxDMC[indi], if ( toDoIsPostMode ) maxLPD=maxLogPostDens[indi] , 
          if ( toDoIsApproxML ) maxLL=maxLogLike[indi] , 
          if ( toDoIsApproxMCL ) maxLCL=maxLogClassLike[indi] , BIC=BicDMC[indi], AIC=AicDMC[indi], AWE=AweDMC[indi], CLC=ClcDMC[indi], IclBic=IclBicDMC[indi], 
          adjRand=corrClassDMC[indi])

print(MSCritTable)
cat("\n\n")

if ( toDoIsPostMode  ) MSCritList$postMode  <- MSCritTable
if ( toDoIsApproxML  ) MSCritList$approxML  <- MSCritTable
if ( toDoIsApproxMCL ) MSCritList$approxMCL <- MSCritTable

# save.image( paste("DMCModSelCrits_",whatToDo,".RData",sep="") )
save(list = c("outFileNamesDMC", "Hmax", "whatToDo", "margDataLL_0", "BicDMC", "AicDMC", "AweDMC", 
                                        "IclBicDMC", "ClcDMC", "corrClassDMC", "mMaxDMC", "HH", "Hnr", "xval", "myLabel", "MSCritList",
                                        "toDoIsPostMode", "toDoIsApproxML", "toDoIsApproxMCL", "prevwd", "indi",
                                        if ( toDoIsPostMode ) "maxLogPostDens",
                                        if ( toDoIsApproxML ) "maxLogLike",
                                        if ( toDoIsApproxMCL ) "maxLogClassLike"
                                        ), file = paste("DMCModSelCrits_",whatToDo,".RData",sep="") )

} # end for whatToDo

print(cbind(outFileNames=outFileNamesDMC), quote=FALSE)
    
MSCritList$outFileNames <- outFileNamesDMC

setwd(prevwd)

return( invisible( MSCritList ) )

}
