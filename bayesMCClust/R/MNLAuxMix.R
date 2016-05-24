
MNLAuxMix <- function( Data = list( storeDir = "try01", # will be created if not existing (in current working directory!)
                                    X = stop("X (matrix of covariates) must be specified")), # individual specific variables
                 Prior = list( H = 4, # number of alternatives 1,...,H
                               betaPrior = "informative",  # 'uninformative' (improper) prior pars for beta (betaPriorVar = infty)
                               betaPriorMean = 0,
                               betaPriorVar = 1),   # 'informative' prior
                 Initial = list( S.i.start = rep(1:H, N),   # have to be the natural numbers (starting with 1)
                                 Beta.start = NULL ),  #    vector of multinomial outcomes / choice made
                 Mcmc = list( M = 50, 
                              M0 = 20, 
                              mOut = 5, 
                              mSave = 10, 
                              seed = 12345)
               ) {

stopifnot( min(Initial$S.i.start)==1 )

# store start date and time
startDate <- format( Sys.time(), "%Y%m%d_%H%M%S" )  # eg. "20080603_105216"
startTime <- proc.time()[3]

# store parameter values
H <-     Prior$H         #   number of alternatives (groups)               
M <-     Mcmc$M          #   number of MCMC draws                          
M0 <-    Mcmc$M0         #   number of discarded draws (burn-in)           
set.seed(Mcmc$seed)      #   set random seed                               
mOut <-  Mcmc$mOut       #   report after every <mOut> draws               
mSave <- Mcmc$mSave      #   save Workspace after every <mSave> draws      
dirName  <- Data$storeDir 

# Define covariates
X <- Data$X     # set covariates (incl intercept)
XNo <- ncol(X)  # determine number of covariates (incl intercept)
dimnamesX <- dimnames(X)[[2]]

# create directory 'storeDir' if necessary
if ( !file.exists(dirName) ) dir.create(dirName)

# define file name where the results are saved
fileName <- paste("mnLogit_newAux_H", H, "_M", M, "_", startDate, ".RData", sep="")  # eg. "mnLogit_H4_M10000_20080603_105408.RData"
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
wr  <- c(0.3347, 5.8701, 5.8726, 22.4268, 28.7401, 36.7557)/100
sr2 <- c(15.9234, 8.9109, 0.8468, 5.0772, 1.6100, 2.8904)
sr  <- c(3.9904, 2.9851, 0.9202, 2.2533, 1.2689, 1.7001)
logwr <- log(wr)
logsr <- log(sr)

# store logLike in each iteration
logLike <- numeric(M)

# calculate N from data
N <- dim(X)[1] 

# start value for Beta
Beta.start <- if ( is.null( Initial$Beta.start ) ) matrix( 0, XNo, H, dimnames=list( dimnamesX , dimnamesH )) else Initial$Beta.start
stopifnot( dim(Beta.start) == c(XNo, H) )

#   a)          
# ======        

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

    if ( tooBig ) {  # slower but even feasible     
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

    BkNinv[,,h] <- Bk0inv[,,h] + mySum1 # 
    
    BkN[,,h] <- chol2inv( chol( BkNinv[,,h] ) )    
    
    bkN[,h] <- crossprod( BkN[,,h], mySum2 + crossprod( Bk0inv[,,h], bk0[,h] )    )    

    Beta.m[,h,m] <- mvrnorm( n=1, mu=bkN[,h], Sigma=BkN[,,h] )
}

#------------------------

# store log likelihood: 
# ===================== 
XBetak <- crossprod(t(X), Beta.m[,,m])
logit.start <- exp( XBetak ) / rowSums( exp( XBetak ) ) 

logLike[m] <- sum( log( logit.start[cbind(1:N, Initial$S.i.start)]))   # "observed" log likelihood!!

# report something...
cat("workspaceFile: ", workspaceFile, "  (within current working directory!) \n")
cat("m =", m, "; duration of iter proc so far: ", round( proc.time()[3] - ptm, 2 ), "sec. \n")
flush.console()

# MCMC Sampler
for (m in 2:M) {
    
    # new aux mix: sample uniform random numbers -- used to sample utilities
    Unif <- matrix( runif(N*(H-1)), N, H-1, dimnames=list( NULL , dimnamesH[-1] ) )

    # calculate XBeta and Lambda
    XBeta  <- XBetak   #    XBetak from the previous iteration #    
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

    # store log likelihood: 
    # ===================== 
    XBetak <- crossprod(t(X), Beta.m[,,m])
    logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) 

    logLike[m] <- sum( log( logit.temp[cbind(1:N, Initial$S.i.start)]))   # "observed" log likelihood!!
    
    # report something...
    if ( identical(all.equal(m %% mOut, 0), TRUE) || is.element(m, c(1:5, 10, 20, 50, 100, 200, 500) ) ) {
            cat("m =", m, "; duration of iter proc so far:", 
                            round( diff <- proc.time()[3] - ptm, 2 ), "sec.,  exp time to end:", round( (diff/(m-1)*M - diff)/60, 2 ), " min. \n")
            flush.console()
            }
    # intermediate storage
    if (identical(all.equal(m %% mSave, 0), TRUE))
                save( list = c("Data", "Prior", "Initial", "Mcmc", "Beta.m", "bk0", "Bk0inv", "fileName" ), 
                      file = file.path(dirName, paste(startDate, "_temp.RData", sep="")) )                     
                      

}

# save total time
totalTime <- proc.time()[3] - startTime

# save results in 'workspaceFile'
save( list = c("Data", "Prior", "Initial", "Mcmc", "Beta.m", "bk0", "Bk0inv", "fileName", 
               "N", "totalTime", "workspaceFile", "bkN", "BkN", "logLike" ), file = workspaceFile ) 

# delete intermediate storage
unlink( file.path(dirName, paste(startDate, "_temp.RData", sep="")) )

# report total time
cat("Total time:", totalTime%/%3600, "hours", (totalTime%%3600)%/%60, "min \n")

# build list to be returned
resultsList <- list( workspaceFile=workspaceFile, Data=Data, Prior=Prior, Initial=Initial, Mcmc=Mcmc, Beta.m=Beta.m, bk0=bk0, Bk0inv=Bk0inv, fileName=fileName, N=N, totalTime=totalTime, 
                     bkN=bkN, BkN=BkN, logLike=logLike )

# return results file name
return( invisible( resultsList ) )

}
