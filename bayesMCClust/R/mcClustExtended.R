
mcClustExtended <- function( Data = list( dataFile = stop("'dataFile' (=> Njk.i) must be specified: either 'filename' (path) or data"),
                                          storeDir = "try01", # will be created if not existing (in current working directory!)
                                          priorFile = NULL,   # 'priorFile' (=> mcc) must be specified: either "filename" (path) or data
                                          X = stop("X (matrix of covariates) must be specified")), # 
                 Prior = list( H = 4, 
                               c = 1, 
                               cOff = 1, 
                               usePriorFile = FALSE, 
                               xiPooled = FALSE, 
                               N0 = 5,
                               betaPrior = "informative",  # 'uninformative' (improper) prior pars for beta (betaPriorVar = infty)
                               betaPriorMean = 0,
                               betaPriorVar = 1),   # 'informative' prior
                 Initial = list( xi.start.ind = 3, 
                                 pers = 0.7,
                                 S.i.start = rep(1:H, N),
                                 Beta.start = NULL ),  #    or     sample(1:H, N, replace=TRUE)
                 Mcmc = list( M = 50, 
                              M0 = 20, 
                              mOut = 5, 
                              mSave = 10, 
                              seed = 12345)
               ) {

# store start date and time
startDate <- format( Sys.time(), "%Y%m%d_%H%M%S" )  # eg. "20080603_105216"
startTime <- proc.time()[3]

# store parameter values
H <-     Prior$H         #   number of groups                              
M <-     Mcmc$M          #   number of MCMC draws                          
M0 <-    Mcmc$M0         #   number of discarded draws (burn-in)           
set.seed(Mcmc$seed)      #   set random seed                               
mOut <-  Mcmc$mOut       #   report after every <mOut> draws               
mSave <- Mcmc$mSave      #   save Workspace after every <mSave> draws      
xi.start.ind <- Initial$xi.start.ind   #   set 'mode' for initial values   
pers <-  Initial$pers    #   only needed if  xi.start.ind == 3             
c <-     Prior$c         #   prior for xi                                  
cOff <-  Prior$cOff      #   prior for xi                                  
dirName  <- Data$storeDir 

# Define covariates
X <- Data$X     # set covariates (incl intercept)
XNo <- ncol(X)  # determine number of covariates (incl intercept)
dimnamesX <- dimnames(X)[[2]]

# create directory 'storeDir' if necessary
if ( !file.exists(dirName) ) dir.create(dirName)

# define file name where the results are saved
fileName <- paste("MCC_Logit_newAux_H", H, "_M", M, "_", startDate, ".RData", sep="")  # eg. "MCC_H4_M10000_20080603_105408.RData"
dimnamesH <- paste("h", 1:H, sep="")

# set prior pars for beta (bk0, Bk0inv)
if (Prior$betaPrior == "uninformative") {
    bk0 <-   matrix(0, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"0",sep="") ) )
    Bk0inv <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"0inv",sep="")))
} else {
    bk0 <-   matrix(Prior$betaPriorMean, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"0",sep="") ) )
    Bk0inv <- array( chol2inv( chol( diag(Prior$betaPriorVar, XNo) ) ), c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"0inv",sep="")))
}

# define path to results file (in current working directory)
workspaceFile <- file.path(dirName, fileName)  # eg. "try02\\MCC_H4_M10000_20080603_105408.RData"

# normal mixture approximation for the logistic distribution 

# 6 components
wr  <- c(0.3347, 5.8701, 5.8726, 22.4268, 28.7401, 36.7557)/100
sr2 <- c(15.9234, 8.9109, 0.8468, 5.0772, 1.6100, 2.8904)
sr  <- c(3.9904, 2.9851, 0.9202, 2.2533, 1.2689, 1.7001)

logwr <- log(wr)
logsr <- log(sr)

# store logLike in each iteration
logLike <- numeric(M)
logClassLike <- numeric(M)
entropy <- numeric(M)
logBetaPrior <- numeric(M) 
logXiPrior <- numeric(M)

# define several function for later use
freqMat  <- function(i) S.i.counts[i, S.i.temp[i] ] <<- S.i.counts[i, S.i.temp[i] ] + 1

countMat  <- function(h) Njk.h.temp[,,h] <<- apply( Njk.i.ind[,,S.i.temp.ind==h] , c(1, 2), sum)   # to replace with tapply!?

dirDensFun <- function(h, j) log( ddirichlet( xi.m[m, j, , h ], c0[j, , h] ) )  # for calculation of prior densitiy of xi

# load data depending on format of dataFile
if ( is.character(Data$dataFile) ) {
    tp <- load(Data$dataFile)
    stopifnot( tp == "Njk.i" ) 
} else {
    Njk.i <- Data$dataFile
}

# calculate N and K from data
N <- dim(Njk.i)[3]
K <- dim(Njk.i)[1] - 1   

# load prior data (if required) depending on format of priorFile
if (Prior$usePriorFile | Initial$xi.start.ind == 4) { 
    if ( is.character(Data$priorFile) ) {
        tp2 <- load(Data$priorFile)
        stopifnot( tp2 == "mcc" ) 
    } else {
        mcc <- Data$priorFile
    } 
} else { 
    mcc <- NULL 
} 

# sum up all transitions of all individuals 
Njk <- matrix(0, K + 1, K + 1, dimnames = dimnames(Njk.i)[1:2] )
Njk <- apply(Njk.i, c(1, 2), sum) 

# classical estimation of transition matrix (-> relative frequencies or empirical transition matrix)
xi.hat <- Njk/apply(Njk, 1, sum) # scaled to row sum = 1 

# matrix-matching (pattern recognition)
pasteNjk <- apply(Njk.i, c(3), paste, collapse=" ") # by column!
names(pasteNjk) <- 1:N 
freq <- as.numeric( table(pasteNjk) ) # frequencies of different (!) transition matrices (in ascending order)
indizes <- as.numeric( names( sort( pasteNjk[!duplicated(pasteNjk)] ) ) ) # indices of different (!) transition matrices
R <- length(freq) # number of different (!) transition matrices

Njk.i.ind <- array(0, c(K + 1, K + 1, R), dimnames = dimnames(Njk.i) )
for (asd in 1:R) Njk.i.ind[,,asd] <- Njk.i[,,indizes[asd]]*freq[asd]

# prepare prior information
c0 <- 
if (Prior$usePriorFile) {    
    if (Prior$xiPooled) {            
            array( ceiling(Prior$N0*mcc[[1]]$xi), c(K + 1, K + 1, H) )            
            } else {    
            array( ceiling(Prior$N0*mcc[[H]]$xi), c(K + 1, K + 1, H) )            
            }    
    } else {    
        array( diag(c, K + 1) + cOff, c(K + 1, K + 1, H) ) 
}

# prepare initial values

if (xi.start.ind == 1) {
    xi.h <- array( 1/(K + 1), c(K + 1, K + 1, H)) # 'uniform' transition matrix in group h
}
if (xi.start.ind == 2) {
    xi.h <- array(0, c(K + 1, K + 1, H))
    invisible(sapply(1:H, function(h) xi.h[,,h] <<- xi.hat ))      # empirical transition matrix
}
if (xi.start.ind == 3) {
    trans <- (1 - pers)/K
    xi.start <- (pers - trans)*diag(1, K + 1, K + 1) + trans*matrix(1, K + 1, K + 1)
    xi.h <- array(0, c(K + 1, K + 1, H))
    invisible(sapply(1:H, function(h) xi.h[,,h] <<- xi.start ))      #  diagonal transition matrix  
}
if (xi.start.ind == 4) {
    xi.h <- array( mcc[[H]]$xi, c(K + 1, K + 1, H))    #  use prior data
}

# start value for Beta
Beta.start <- if ( is.null( Initial$Beta.start ) ) matrix( 0, XNo, H, dimnames=list( dimnamesX , dimnamesH )) else Initial$Beta.start
stopifnot( dim(Beta.start) == c(XNo, H) )

#   a)          
# ======        

# to store frequencies of (individual) classifications
S.i.counts <- matrix(0, N, H)

# initialise posterior pars for beta
bkN <- matrix(0, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"N",sep="") ) )
BkNinv <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"Ninv",sep="")))
BkN    <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"N",sep="")))

# initialise (unknown) parameter matrix
Beta.m <- array( 0, c(XNo, H, M), dimnames=list( dimnamesX , dimnamesH, NULL ))

# first iteration of MCMC sampler                                  
ptm <- proc.time()[3]

###########
  m <- 1  #
###########

# indicator matrix for group membership (Si != h) !!!
Ind <- matrix(1, N, H, dimnames=list( NULL , dimnamesH ) )

# set indicator matrix  ( 10000 executions take approx 600 seconds = 10 minutes !! )
for (i in 1:N) Ind[ i, Initial$S.i.start[i] ] <- 0  
IndCut <- Ind[,2:H]

# new aux mix: sample uniform random numbers -- used to sample utilities
Unif <- matrix( runif(N*(H-1)), N, H-1, dimnames=list( NULL , dimnamesH[-1] ) )

# calculate XBeta and Lambda
XBeta  <- crossprod(t(X), Beta.start ) 
Lambda <- exp( XBeta )

# new aux mix: 
LambdaTemp <- matrix(rowSums(Lambda), nrow=N, ncol=H) - Lambda
LambdaTempCut <- LambdaTemp[,2:H]
LambdaStar <- Lambda[,2:H]/LambdaTempCut  

# NEW aux mix: calculate/draw utilities (
Util <- log( LambdaStar * Unif + ( 1 - IndCut ) ) - log( 1 - Unif + LambdaStar * IndCut  )

# calculate some other auxiliary variables
bnm <- Util - XBeta[,2:H]  + log( LambdaTempCut ) # of dim N x (H-1)

#---------------------------

Rih <- matrix(0,N,H-1)

tTrickmat <- t(outer(seq_len(length(wr)),seq_len(length(wr)), "<="))

if (H > 2) {
    logrProbs <- apply( bnm , c(1,2), function(x) logwr - logsr - 1/2*(x/sr)^2 ) 
    rProbs <- exp( logrProbs )
    for (h in 1:(H-1)) {
        vertfkt <- tTrickmat %*% rProbs[,,h]    
        Rih[,h] <- rowSums( vertfkt[length(wr),]*runif(N) >  t(vertfkt) )+1  
    } 
} else { 
    logrProbs <- sapply( bnm , function(x) logwr - logsr - 1/2*(x/sr)^2 ) 
    rProbs <- exp( logrProbs )
    vertfkt <- tTrickmat %*% rProbs
    Rih[,1] <- rowSums( vertfkt[length(wr),]*runif(N) >  t(vertfkt) )+1
}

#---------------------------

SRmat  <- matrix( sr[Rih], N, H - 1, dimnames=list(NULL, dimnamesH[-1]) )
SR2mat <- matrix(sr2[Rih], N, H - 1, dimnames=list(NULL, dimnamesH[-1]) )

# another auxiliary variable
UtilStd <- (Util  + log( LambdaTempCut ) )/SR2mat  

#---------------------------

tooBig <- ( N > 50000 )

if ( !tooBig ) { # faster but needs much more memory space 
    myXtX <- matrix(0, N, XNo*XNo) 
    for ( i in 1:N ) {    
        myXtX[i,] <- as.vector(tcrossprod( X[i,] ))    
    }
}

for ( h in 2:H ) {    
    
    if ( tooBig ) { # slower but even feasible   
        XStd <- X/sqrt(SR2mat[,h-1])  
        mySum1 <- matrix(0, XNo, XNo)  
        for ( i in 1:N ) { 
            mySum1 <- mySum1 + tcrossprod( XStd[i,] )     
        }            
    } else { # faster but needs much more memory space    
        mySum1 <- matrix( colSums( myXtX / SR2mat[,h-1] ), nrow=XNo, ncol=XNo )            
    }  

    myXtW <- X * UtilStd[,h-1]  
    mySum2 <- colSums( myXtW )
    BkNinv[,,h] <- Bk0inv[,,h] + mySum1     
    BkN[,,h] <- chol2inv( chol( BkNinv[,,h] ) )        
    bkN[,h] <- crossprod( BkN[,,h], mySum2 + crossprod( Bk0inv[,,h], bk0[,h] )    )    
    Beta.m[,h,m] <- mvrnorm( n=1, mu=bkN[,h], Sigma=BkN[,,h] )
}

#------------------------

## Pr(Si=h|betah,h=1,...,H)
XBetak <- crossprod(t(X), Beta.m[,,m])   
logit.start <- exp( XBetak ) / rowSums( exp( XBetak ) ) # N x H -matrix

# Bayes' classification for each individual i 
# =========================================== 

trickmat <- outer(seq_len(H),seq_len(H), "<=")

sProbs <-  matrix( mapply( function(i, h) prod( xi.h[,,h]^Njk.i[,,i] ), rep(1:N, each=H), 1:H ), N, H, byrow=TRUE ) * logit.start

if (H > 1) {    
            vertfkt <- sProbs %*% trickmat
            S.i.temp <- rowSums( vertfkt[,H]*runif(N) > vertfkt )+1        
} else {        
            S.i.temp <- rep(1,N) # individual classification 
}       

if (M0 <= 1) invisible( sapply(1:N, function(i) S.i.counts[i, S.i.temp[i] ] <<- S.i.counts[i, S.i.temp[i] ] + 1) )   # TO BE IMPROVED!?!?!

#   b)          
# ======        
S.i.temp.ind <- S.i.temp[indizes]

#   c)          
# ======        
xi.m <- array(0, c(M, K + 1, K + 1, H))
Njk.h.temp <- array(0, c(K + 1, K + 1, H) )

# calculate absolute transition frequencies in each of H groups
invisible( sapply( 1:H, function(h) Njk.h.temp[,,h] <<- apply( Njk.i.ind[,,S.i.temp.ind==h] , c(1, 2), sum) ) )

# first draw of transition matrices
hj.grid <- as.matrix(expand.grid(h=1:H, j=1:(K + 1)))
transMat <- function(h, j) xi.m[m, j, , h] <<- gtools::rdirichlet(n= 1, alpha= Njk.h.temp[ j,,h] + c0[j,,h] )
invisible( mapply( transMat, hj.grid[,1], hj.grid[,2] ) )

# store log likelihood: 
# ===================== 

# observed log likelihood
tempMix <- matrix(mapply(function(i,h) prod(xi.m[m,,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE )
tempLK  <- logit.start * tempMix

tik <- tempLK/rowSums(tempLK)
sProbsMax <- max.col(tempLK) 

logLike[m] <- sum( log( rowSums( tempLK ) ) )   # "observed" log likelihood!!

logClassLike[m] <- sum( log( ( tempLK )[ cbind( 1:N, sProbsMax ) ] ) )   # classification log likelihood!!

entropy[m] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )

logBetaPrior[m] <- sum( dnorm( Beta.m[ , -1, m], mean = Prior$betaPriorMean, sd = sqrt(Prior$betaPriorVar), log=TRUE ) )

# calculate prior densities for xi and eta
logXiPrior[m] <- sum( mapply( dirDensFun , hj.grid[,1], hj.grid[,2]  ) )

# report something...
cat("workspaceFile: ", workspaceFile, "  (within current working directory!) \n")
cat("m =", m, "; duration of iter proc so far: ", round( proc.time()[3] - ptm, 2 ), "sec. \n")
flush.console()

# MCMC Sampler
for (m in 2:M) { 
    
    # indicator matrix for group membership (Si != h) !!!
    Ind <- matrix(1, N, H, dimnames=list( NULL , dimnamesH ) )
    
    # set indicator matrix  ( 10000 executions take approx 600 seconds = 10 minutes !! )
    for (i in 1:N) Ind[ i, S.i.temp[i] ] <- 0  
    IndCut <- Ind[,2:H]

    # new aux mix: sample uniform random numbers -- used to sample utilities
    Unif <- matrix( runif(N*(H-1)), N, H-1, dimnames=list( NULL , dimnamesH[-1] ) )

    # calculate XBeta and Lambda
    XBeta  <- XBetak    #    XBetak from the previous iteration    
    Lambda <- exp( XBeta )
    
    # new aux mix: 
    LambdaTemp <- matrix(rowSums(Lambda), nrow=N, ncol=H) - Lambda
    LambdaTempCut <- LambdaTemp[,2:H]
    LambdaStar <- Lambda[,2:H]/LambdaTempCut  

    # NEW aux mix: calculate/draw utilities (
    Util <- log( LambdaStar * Unif + ( 1 - IndCut ) ) - log( 1 - Unif + LambdaStar * IndCut  )

    # calculate some other auxiliary variables
    bnm <- Util - XBeta[,2:H] + log( LambdaTempCut )  # of dim N x (H-1) 
    
    if (H > 2) {
        logrProbs <- apply( bnm , c(1,2), function(x) logwr - logsr - 1/2*(x/sr)^2 ) 
        rProbs <- exp( logrProbs )
        for (h in 1:(H-1)) {
            vertfkt <- tTrickmat %*% rProbs[,,h]    
            Rih[,h] <- rowSums( vertfkt[length(wr),]*runif(N) >  t(vertfkt) )+1  
        }
    } else { 
        logrProbs <- sapply( bnm , function(x) logwr - logsr - 1/2*(x/sr)^2 ) 
        rProbs <- exp( logrProbs )
        vertfkt <- tTrickmat %*% rProbs   
        Rih[,1] <- rowSums( vertfkt[length(wr),]*runif(N) >  t(vertfkt) )+1          
    }

    SRmat  <- matrix( sr[Rih], N, H - 1, dimnames=list(NULL, dimnamesH[-1]) )
    SR2mat <- matrix(sr2[Rih], N, H - 1, dimnames=list(NULL, dimnamesH[-1]) )

    # another auxiliary variable
    UtilStd <- (Util + log( LambdaTempCut ) )/SR2mat   
    
    for ( h in 2:H ) {    
        
        if ( tooBig ) {    
            XStd <- X/sqrt(SR2mat[,h-1])  
            mySum1 <- matrix(0, XNo, XNo)  
            for ( i in 1:N ) { 
                mySum1 <- mySum1 + tcrossprod( XStd[i,] )     
            }            
        } else {    
            mySum1 <- matrix( colSums( myXtX / SR2mat[,h-1] ), nrow=XNo, ncol=XNo )          
        }         
        
        myXtW <- X * UtilStd[,h-1]   
        mySum2 <- colSums( myXtW )
        BkNinv[,,h] <- Bk0inv[,,h] + mySum1 
        BkN[,,h] <- chol2inv( chol( BkNinv[,,h] ) )      
        bkN[,h] <- crossprod( BkN[,,h], mySum2 + crossprod( Bk0inv[,,h], bk0[,h] )    )    
        Beta.m[,h,m] <- mvrnorm( n=1, mu=bkN[,h], Sigma=BkN[,,h] )
    }
    
    ## Pr(Si=h|betah,h=1,...,H)
    XBetak <- crossprod(t(X), Beta.m[,,m])   
    logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) # N x H -matrix
    
    # Bayes' classification for each individual i 
    # =========================================== 
    sProbs <- tempMix * logit.temp   # tempMix from the previous iteration # logit.temp from the current iteration !!!!    
    
    if (H > 1) {    
            vertfkt <- sProbs %*% trickmat
            S.i.temp <- rowSums( vertfkt[,H]*runif(N) > vertfkt )+1        
        } else {        
            S.i.temp <- rep(1,N) # individual classification 
    }       
    
    # aggregate (for technical (storage) reasons) classifications from M0 on 
    if (m > M0) sapply(1:N, freqMat ) 
    
    # sample mixing proportions (group sizes)
    S.i.temp.ind <- S.i.temp[indizes]
    
    # transition matrix 
    Njk.h.temp <- array(0, c(K + 1, K + 1, H) )  
    sapply(1:H, countMat )
    mapply( transMat, hj.grid[,1], hj.grid[,2] )    

    # store log likelihood: 
    # ===================== 

    # observed log likelihood
    tempMix <- matrix(mapply(function(i,h) prod(xi.m[m,,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE )               
    tempLK  <- logit.temp * tempMix
    
    tik <- tempLK/rowSums(tempLK)
    sProbsMax <- max.col(tempLK)    
    
    logLike[m] <- sum( log( rowSums( tempLK ) ) )   # "observed" log likelihood!!
    
    logClassLike[m] <- sum( log( ( tempLK )[ cbind( 1:N, sProbsMax ) ] ) )   # classification log likelihood!!
    
    entropy[m] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )
    
    logBetaPrior[m] <- sum( dnorm( Beta.m[ , -1, m], mean = Prior$betaPriorMean, sd = sqrt(Prior$betaPriorVar), log=TRUE ) )

    # calculate prior densities for xi and eta
    logXiPrior[m] <- sum( mapply( dirDensFun , hj.grid[,1], hj.grid[,2]  ) )
 
    # report something...
    if ( identical(all.equal(m %% mOut, 0), TRUE) || is.element(m, c(1:5, 10, 20, 50, 100, 200, 500) ) ) {
            cat("m =", m, "; duration of iter proc so far:", 
                            round( diff <- proc.time()[3] - ptm, 2 ), "sec.,  exp time to end:", round( (diff/(m-1)*M - diff)/60, 2 ), " min. \n")
            flush.console()
            }
    # intermediate storage
    if (identical(all.equal(m %% mSave, 0), TRUE))
                save( list = c("Data", "Prior", "Initial", "Mcmc", "Beta.m", "bk0", "Bk0inv", "c0", "fileName",
                    "freq", "indizes", "K", "mcc", "N", "Njk.i", "Njk.i.ind", "R", "S.i.counts", "workspaceFile",
                    "xi.hat", "xi.m", if (exists("xi.start")) "xi.start", "xi.start.ind", "bkN", "BkN", 
                    "logLike", "logBetaPrior", "logXiPrior", "logClassLike", "entropy" ), 
                    file = file.path(dirName, paste(startDate, "_temp.RData", sep="")) )
}

# store ergodic average of xi_h for all groups
estTransProb <- if (H > 1) apply( xi.m[M0:M,,,], c(2, 3, 4), mean) else array( apply( xi.m[M0:M,,,], c(2, 3), mean), c(K+1,K+1,H))
dimnames(estTransProb) <- dimnames(Njk.i)
dimnames(estTransProb)[[3]] <- paste("Group", 1:H)

logPostDens <- logLike + logXiPrior + if ( Prior$betaPrior!="uninformative" ) logBetaPrior else 0 
mMax <- which.max(logPostDens)

# save total time
totalTime <- proc.time()[3] - startTime

# save results in 'workspaceFile'
save( list = c("Data", "Prior", "Initial", "Mcmc", "Beta.m", "bk0", "Bk0inv", "c0", "estTransProb", "fileName", 
               "freq", "indizes", "K", "mcc", "N", "Njk.i", "Njk.i.ind", "R", "S.i.counts", "totalTime", "workspaceFile",
               "xi.hat", "xi.m", if (exists("xi.start")) "xi.start", "xi.start.ind", "bkN", "BkN", 
               "logLike", "logBetaPrior", "logXiPrior", "logPostDens", "mMax", "logClassLike", "entropy" ), 
               file = workspaceFile ) 

# delete intermediate storage
unlink( file.path(dirName, paste(startDate, "_temp.RData", sep="")) )

# report total time
cat("Total time:", totalTime%/%3600, "hours", (totalTime%%3600)%/%60, "min \n")

# build list to be returned
resultsList <- list( workspaceFile=workspaceFile, Data=Data, Prior=Prior, Initial=Initial, Mcmc=Mcmc, Beta.m=Beta.m, bk0=bk0, Bk0inv=Bk0inv, c0=c0, estTransProb=estTransProb, 
               fileName=fileName, freq=freq, indizes=indizes, K=K, mcc=mcc, N=N, Njk.i=Njk.i, Njk.i.ind=Njk.i.ind, R=R, S.i.counts=S.i.counts, totalTime=totalTime,
               xi.hat=xi.hat, xi.m=xi.m, xi.start= if (exists("xi.start")) xi.start else NULL, xi.start.ind=xi.start.ind, bkN=bkN, BkN=BkN, logLike=logLike, logBetaPrior=logBetaPrior, 
               logXiPrior=logXiPrior, logPostDens=logPostDens, mMax=mMax, logClassLike=logClassLike, entropy=entropy )

# return results file name
return( invisible( resultsList ) )

}
