
mcClust <- function( Data = list( dataFile = stop("'dataFile' (=> Njk.i) must be specified: either 'filename' (path) or data"),
                                  storeDir = "try01", # will be created if not existing (in current working directory!)
                                  priorFile = NULL),  # 'priorFile' (=> mcc) must be specified: either "filename" (path) or data
                 Prior = list( H = 4, 
                               e0 = 4, 
                               c = 1, 
                               cOff = 1, 
                               usePriorFile = FALSE, 
                               xiPooled = FALSE, 
                               N0 = 5), 
                 Initial = list( xi.start.ind = 3, 
                                 pers = 0.7, 
                                 S.i.start = NULL ),
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
e0 <-    Prior$e0        #   prior for eta                                 
c <-     Prior$c         #   prior for xi                                  
cOff <-  Prior$cOff      #   prior for xi                                  
dirName  <- Data$storeDir 

# create directory 'storeDir' if necessary
if ( !file.exists(dirName) ) dir.create(dirName)

# define file name where the results are saved
fileName <- paste("MCC_H", H, "_M", M, "_", startDate, ".RData", sep="")  # eg. "MCC_H4_M10000_20080603_105408.RData"

# define path to results file (in current working directory)
workspaceFile <- file.path(dirName, fileName)  # eg. "try02\\MCC_H4_M10000_20080603_105408.RData"

# store logLike and priors in each iteration
logLike <- numeric(M)
logClassLike <- numeric(M)
entropy <- numeric(M)
logXiPrior <- numeric(M)
logEtaPrior <- numeric(M) 

# define several function for later use
NhSfunPat <- function(S, H, freq) { 
    sapply( 1:H ,  function(h) sum( (S==h)*freq ) )
}

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
            array( (Prior$N0*mcc[[1]]$xi), c(K + 1, K + 1, H) )            
            } else {    
            array( (Prior$N0*mcc[[H]]$xi), c(K + 1, K + 1, H) )            
            }    
    } else {    
        array( diag(c, K + 1) + cOff, c(K + 1, K + 1, H) ) 
}

# prepare initial values
eta.start <- 1/H  # relative group size

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
    eta.start <- mcc[[H]]$eta
}

eta.h <- array(eta.start, H)

#   a)          
# ======        

# to store frequencies of (individual) classifications
S.i.counts <- matrix(0, N, H)

# first iteration of MCMC sampler                                  
ptm <- proc.time()[3]

###########
  m <- 1  #
###########

# Bayes' classification for each individual i 
# =========================================== 

trickmat <- outer(seq_len(H),seq_len(H), "<=")
sProbs <-  matrix( mapply( function(i, h) prod( xi.h[,,h]^Njk.i[,,i] ), rep(1:N, each=H), 1:H ), N, H, byrow=TRUE ) * matrix(eta.h, N, H, byrow=TRUE) # dim = N x H

if ( is.null( Initial$S.i.start ) ) {
    if (H > 1) {    
            vertfkt <- sProbs %*% trickmat
            S.i.temp <- rowSums( vertfkt[,H]*runif(N) > vertfkt )+1        
        } else {        
            S.i.temp <- rep(1,N) # individual classification 
    }       
} else {
    S.i.temp <- if (H > 1) Initial$S.i.start else rep(1,N)
}

if (M0 <= 1) invisible( sapply(1:N, function(i) S.i.counts[i, S.i.temp[i] ] <<- S.i.counts[i, S.i.temp[i] ] + 1) )   

#   b)          
# ======        
eta.m <- matrix(0, H, M )
S.i.temp.ind <- S.i.temp[indizes]
eta.m[, m] <- if (H > 1) gtools::rdirichlet(n= 1, alpha= NhSfunPat(S.i.temp.ind, H, freq) + e0 ) else 1  # sample mixing proportions

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

# calculate prior densities for xi and eta
logXiPrior[m] <- sum( mapply( dirDensFun , hj.grid[,1], hj.grid[,2]  ) )
logEtaPrior[m] <- log( ddirichlet( eta.m[,m], rep(e0 ,H) ) )

# report something...
cat("workspaceFile: ", workspaceFile, "  (within current working directory!) \n")
cat("m =", m, "; duration of iter proc so far: ", round( proc.time()[3] - ptm, 2 ), "sec. \n")
flush.console()

# MCMC Sampler
for (m in 2:M) {

    # Bayes' classification for each individual i 
    # =========================================== 
    
    sProbsPart1 <- matrix( mapply( function(i, h) prod( xi.m[m - 1,,,h]^Njk.i[,,i] ), rep(1:N, each=H), 1:H ), N, H, byrow=TRUE )
    
    sProbs <- sProbsPart1 * matrix(eta.m[, m - 1], N, H, byrow=TRUE)                                                                                       
    tik <- sProbs/rowSums(sProbs)
    sProbsMax <- max.col(sProbs)
                                                                                                                         
    # individual classification 
    if (H > 1) {    
        vertfkt <- sProbs %*% trickmat
        S.i.temp <- rowSums( vertfkt[,H]*runif(N) > vertfkt )+1        
    } else {        
        S.i.temp <- rep(1,N) # individual classification 
    }       
    
    # store log likelihood: 
    # ===================== 

    # observed and classification log likelihood
    logLike[m-1]      <- sum( log( rowSums( sProbs ) ) )   # "observed" log likelihood!!
    logClassLike[m-1] <- sum( log( ( sProbsPart1 )[cbind(1:N,sProbsMax)]* eta.m[ sProbsMax , m - 1 ] ))    
    
    entropy[m-1] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )
    
    # aggregate (for technical (storage) reasons) classifications from M0 on 
    if (m > M0) sapply(1:N, freqMat ) 
    
    # sample mixing proportions (group sizes)
    S.i.temp.ind <- S.i.temp[indizes]
    eta.m[, m] <- if (H > 1) gtools::rdirichlet(n= 1, alpha= NhSfunPat(S.i.temp.ind, H, freq) + e0 ) else 1
    
    # transition matrix 
    Njk.h.temp <- array(0, c(K + 1, K + 1, H) )  
    sapply(1:H, countMat )
    mapply( transMat, hj.grid[,1], hj.grid[,2] )
    
    # calculate prior densities for xi and eta
    logXiPrior[m] <- sum( mapply( dirDensFun , hj.grid[,1], hj.grid[,2]  ) )
    logEtaPrior[m] <- log( ddirichlet( eta.m[,m], rep(e0 ,H) ) )
    
    # report something...
    if ( identical(all.equal(m %% mOut, 0), TRUE) || is.element(m, c(1:5, 10, 20, 50, 100, 200, 500) ) ) {
            cat("m =", m, "; duration of iter proc so far:", 
                            round( diff <- proc.time()[3] - ptm, 2 ), "sec.,  exp time to end:", round( (diff/(m-1)*M - diff)/60, 2 ), " min. \n")
            flush.console()
            }
    # intermediate storage
    if (identical(all.equal(m %% mSave, 0), TRUE))
                save( list = c("Data", "Prior", "Initial", "Mcmc", "c0", "eta.start", "eta.m", "fileName",
                    "freq", "indizes", "K", "mcc", "N", "Njk.i", "Njk.i.ind", "R", "S.i.counts", "workspaceFile",
                    "xi.hat", "xi.m", if (exists("xi.start")) "xi.start", "xi.start.ind", 
                    "logLike", "logClassLike", "entropy", "logXiPrior", "logEtaPrior"),
                          file = file.path(dirName, paste(startDate, "_temp.RData", sep="")) )

}

# store log likelihood for m=M: 
# ============================= 

# observed and classification log likelihood

estGroupSizeMat <- matrix(eta.m[, m], N, Prior$H, byrow=TRUE)
sProbsPart1 <-     matrix(mapply(function(i,h) prod(xi.m[m,,,h]^Njk.i[,,i]),rep(1:N,each=Prior$H),1:Prior$H ), N, Prior$H, byrow=TRUE )
tempLK <- estGroupSizeMat * sProbsPart1
tik <- tempLK/rowSums(tempLK)
sProbsMax <- max.col(tempLK) 

logLike[m]      <- sum( log( rowSums( tempLK ) ) )   # "observed" log likelihood!!
logClassLike[m] <- sum( log( ( sProbsPart1 )[cbind(1:N,sProbsMax)]* eta.m[ sProbsMax , m ] ))

entropy[m] <-  if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )

# store ergodic average of xi_h for all groups
estTransProb <- if (H > 1) apply( xi.m[M0:M,,,], c(2, 3, 4), mean) else array( apply( xi.m[M0:M,,,], c(2, 3), mean), c(K+1,K+1,H))
dimnames(estTransProb) <- dimnames(Njk.i)
dimnames(estTransProb)[[3]] <- paste("Group", 1:H)

# store ergodic average of eta_h for all groups
estGroupSize <- if (H > 1) apply( eta.m[,M0:M], 1, mean ) else mean( eta.m[,M0:M] )

logPostDens <- logLike + logXiPrior + logEtaPrior
mMax <- which.max(logPostDens)

# save total time
totalTime <- proc.time()[3] - startTime

# save results in 'workspaceFile'
save( list = c("Data", "Prior", "Initial", "Mcmc", "c0", "eta.start", "estGroupSize", "estTransProb", "eta.m", "fileName",
               "freq", "indizes", "K", "mcc", "N", "Njk.i", "Njk.i.ind", "R", "S.i.counts", "totalTime", "workspaceFile",
               "xi.hat", "xi.m", if (exists("xi.start")) "xi.start", "xi.start.ind", "logLike", "logClassLike", "entropy",
               "logXiPrior", "logEtaPrior", "logPostDens", "mMax" ), 
               file = workspaceFile ) 
               
# delete intermediate storage
unlink( file.path(dirName, paste(startDate, "_temp.RData", sep="")) )

# report total time
cat("Total time:", totalTime%/%3600, "hours", (totalTime%%3600)%/%60, "min \n")

# build list to be returned
resultsList <- list( workspaceFile=workspaceFile, Data=Data, Prior=Prior, Initial=Initial, Mcmc=Mcmc, c0=c0, eta.start=eta.start, estGroupSize=estGroupSize, estTransProb=estTransProb, 
               eta.m=eta.m, fileName=fileName, freq=freq, indizes=indizes, K=K, mcc=mcc, N=N, Njk.i=Njk.i, Njk.i.ind=Njk.i.ind, R=R, S.i.counts=S.i.counts, totalTime=totalTime,
               xi.hat=xi.hat, xi.m=xi.m, xi.start= if (exists("xi.start")) xi.start else NULL, xi.start.ind=xi.start.ind, logLike=logLike, logClassLike=logClassLike, entropy=entropy,
               logXiPrior=logXiPrior, logEtaPrior=logEtaPrior, logPostDens=logPostDens, mMax=mMax )

# return results file name
return( invisible( resultsList ) )

}
