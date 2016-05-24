
dmClustExtended <- function( Data = list(dataFile = stop("'dataFile' must be specified: filename or data"),
                                 storeDir = "try01", 
                                 mccFile = "mcc.RData",
                                 X = stop("X (matrix of covariates) must be specified")), 
                     Prior = list(H = 4, 
                                  a0 = 1, 
                                  alpha = 1, 
                                  N0 = 10,     
                                  isPriorNegBin = FALSE, 
                                  mccAsPrior = FALSE, 
                                  xiPooled = TRUE, 
                                  persPrior = 7/10,
                                  betaPrior = "informative",  # 'uninformative' (improper) prior pars for beta (betaPriorVar = infty)
                                  betaPriorMean = 0,
                                  betaPriorVar = 1),   # 'informative' prior
                     Initial = list(mccUse = FALSE, 
                                    pers = 1/6,
                                    S.i.start = rep(1:H, N),
                                    Beta.start = NULL ),    #    or     sample(1:H, N, replace=TRUE)....
                     Mcmc = list(kNo = 2, 
                                 M = 50, 
                                 M0 = 20, 
                                 mOut = 5, 
                                 mSave = 10, 
                                 showAcc = TRUE, 
                                 monitor = FALSE, 
                                 seed = 12345)                
) {

# start point to measure (total) used time 
pts <- proc.time()[3]

# save date of function call 
startDate <- format( Sys.time(), "%Y%m%d_%H%M%S" )

# store parameter values
H       <- Prior$H      # number of groups 
kNo     <- Mcmc$kNo     # how many updates per row? 
M       <- Mcmc$M       # how many iterations?
M0      <- Mcmc$M0      # burn-in 
mOut    <- Mcmc$mOut    # report after every mOut iterations plus first five
mSave   <- Mcmc$mSave   # save workspace after every mSave iterations
showAcc <- Mcmc$showAcc # report acceptance rate during iteration process
monitor <- Mcmc$monitor # show figures of draws during iteration process
seed    <- Mcmc$seed    # set seed
set.seed(seed)

storeDir <- Data$storeDir # paste(Data$storeDir, "/", sep="")

# Define covariates
X <- Data$X     # set covariates (incl intercept)
XNo <- ncol(X)  # determine number of covariates (incl intercept)
dimnamesX <- dimnames(X)[[2]]

isPriorNegBin <- Prior$isPriorNegBin  # prior information for group parameters e_h (NegMultinom)
a0            <- Prior$a0             # a0, alpha, N0 is for Negative Multinomial-prior
alpha         <- Prior$alpha
N0            <- Prior$N0

mccAsPrior <- Prior$mccAsPrior  # prior parameter xi 
xiPooled   <- Prior$xiPooled
persPrior  <- Prior$persPrior

mccUse     <- Initial$mccUse   # should preliminary results be used for initial values?
if (!mccUse) pers <- Initial$pers 

# create directory 'storeDir' if necessary
if (!file.exists(Data$storeDir)) dir.create(Data$storeDir)

# define file name for workspace to be saved
fileName    <- paste("DMC_Logit_newAux_H", H, "_M", M, "_", startDate, ".RData", sep="")
fileNameTxt <- paste("DMC_Logit_newAux_H", H, "_M", M, "_", startDate, "_logfile.txt", sep="")
workspaceFile <- file.path( storeDir, fileName )

# define log file name
logFileNameTemp <- paste(startDate, "_logfile.txt", sep="")
logFileName <- file.path(storeDir, logFileNameTemp) # paste(storeDir, startDate, "_logfile.txt", sep="")

# store logLike in each iteration
logLike <- numeric(M)
logClassLike <- numeric(M)
entropy <- numeric(M)
logBetaPrior <- numeric(M) 
logEPrior <- numeric(M) 

# normal mixture approximation for the logistic distribution 
wr  <- c(0.3347, 5.8701, 5.8726, 22.4268, 28.7401, 36.7557)/100
sr2 <- c(15.9234, 8.9109, 0.8468, 5.0772, 1.6100, 2.8904)
sr  <- c(3.9904, 2.9851, 0.9202, 2.2533, 1.2689, 1.7001)
logwr <- log(wr)
logsr <- log(sr)

# define several function for later use:

# exp( log( p( y_i | e_h_0 ) ) ) 
margDataLL_0 <- function(i, h, Njk.i=Njk.i, e_h=e_h, AA=AA, BB=BB) {
    CC <- sum( lgamma( Ne <- (Njk.i[,,i] + e_h[,,h]) ) )
    DD <- sum( lgamma( rowSums( Ne ) ) )
    return( exp(AA[h] - BB[h] + CC - DD) )
}

# exp( log( p( y_i | e_h_m ) ) ) 
margDataLL_m <- function(i, h, m, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB) {
    CC <- sum( lgamma( Ne <- (Njk.i[,,i] + e_h_m[,,h,m])   ) )
    DD <- sum( lgamma( rowSums( Ne ) ) )
    return( exp(AA[h] - BB[h] + CC - DD) )
}

# p( e_{h,j.}^{prop} | y, S )  
# ---------------------------  
# p( e_{h,j.}^{old}  | y, S )  
condDistrRowRatio <- function(h, j, e_hj_prop, e_hj_old, m=m, isPriorNegBin=isPriorNegBin, # e_hj_prop=e_hj_new, e_hj_old=e_hj_prev, 
                              N_h_temp=N_h_temp, S_i_h_ind=S_i_h_ind, Njk.i=Njk.i, a0=a0, N0=N0, xi_prior=xi_prior, alpha=alpha ) {    
    logBB <- 0
    if ( N_h_temp >= 1 ) { 
        for ( i in S_i_h_ind ) {
            logBB <- logBB + sum( lgamma( Nep <- (Njk.i[j,,i] + e_hj_prop) ) - lgamma( Neo <- (Njk.i[j,,i] + e_hj_old) ) ) + lgamma( sum( Neo ) ) - 
            lgamma( sum( Nep ) ) 
        } 
    } 
    return(        
    exp( log( priorDens( xVec=e_hj_prop - 1, k=a0, alpha=alpha, lambdaVec=N0*xi_prior[j,,h], isPriorNegBin=isPriorNegBin ) ) - 
         log( priorDens( xVec=e_hj_old  - 1, k=a0, alpha=alpha, lambdaVec=N0*xi_prior[j,,h], isPriorNegBin=isPriorNegBin ) ) +
         N_h_temp*( lgamma( sum( e_hj_prop ) ) - lgamma( sum( e_hj_old ) ) + sum( lgamma( e_hj_old ) - lgamma( e_hj_prop ) )  ) + 
         logBB) )
}

q_k <- function(ehjk_cond) { 
    return( 1/( min(ehjk_cond, 1+1) + 1) )  
}

# update k elements randomly ( H*(K+1) times per iteration !!) 
ehjUpdate_k <- function(ehj, h, j, kNo=kNo, m, isPriorNegBin=isPriorNegBin, K=K, N_h_temp=N_h_temp, S_i_h_ind=S_i_h_ind, Njk.i=Njk.i, 
                        a0=a0, N0=N0, xi_prior=xi_prior, alpha=alpha) {   
    
    e_hj_prev <- ehj    
    kSet <- sample(1:(K+1), size=kNo)     
    updateVector <- rep(0, K+1)    
    for (k in kSet) {
        updateVector[k] <- sample( (min(e_hj_prev[k], 1+1)*(-1)+1):1, 1 ) 
    }    
    e_hj_new  <- e_hj_prev + updateVector      
    r_1 <- condDistrRowRatio(h=h, j=j, e_hj_prop=e_hj_new, e_hj_old=e_hj_prev, m=m, isPriorNegBin=isPriorNegBin, N_h_temp=N_h_temp, 
                             S_i_h_ind=S_i_h_ind, Njk.i=Njk.i, a0=a0, N0=N0, xi_prior=xi_prior, alpha=alpha ) 
    r_2 <- 1
    for ( k in kSet ) { r_2 <- r_2*q_k(ehjk_cond=e_hj_new[k])/q_k(ehjk_cond=e_hj_prev[k]) }     
    accept <- c( min(1, r_1*r_2), 0 )    
    accept[2] <- rbinom(n=1, size=1, prob=accept[1])    
    if (accept[2]) return( list( e_hj=e_hj_new, accept=accept ) ) else return( list( e_hj=e_hj_prev, accept=accept ) )
}

priorDens <- function( xVec, k, alpha=k, lambdaVec=N0*xi_prior[j,,h], isPriorNegBin=isPriorNegBin ) {
    if (isPriorNegBin) { 
        prod( dnbinom( x=xVec, size=k, prob=k/(k + lambdaVec) ) )
    } else {
        dNegMultinom(xVec=xVec, pVec=lambdaVec/(alpha + sum(lambdaVec)), k=k)
    }
}

dNegMultinom <- function(xVec, pVec, k) {   
    exp( lgamma( k + sum(xVec) ) - lgamma(k) - sum( lgamma(xVec + 1) ) + k*log( 1 - sum(pVec) ) + sum( xVec*log(pVec) ) )    
}

monitor_e_h <- function(j, k, e_h_m=e_h_m, ylim=c(0, max(100, max(e_h_m) + 3)), from=1, m=m, M=M, H=H ) { 
    plot(e_h_m[j, k, 1, from:m], type="l", col=2, ylim=ylim, xlim=c(from, min( (m %/% 50 + 1)*50, M) ), xlab="", ylab="", 
                                 main=paste("MCMC draws for", expression(e[h]), "[",j,k,"]") )
    for (h in 2:H) lines( e_h_m[j, k, h, from:m ], type="l", col=h + 1)
}

monitor_xi_h <- function(j, k, xi_h_m=xi_h_m, ylim=c(0, 1), from=1, m=m, M=M, H=H ) { 
    plot(xi_h_m[j, k, 1, from:m], type="l", col=2, ylim=ylim, xlim=c(from, min( (m %/% 50 + 1)*50, M) ), xlab="", ylab="", 
                                 main=paste("MCMC draws for", expression(xi[h]), "[",j,k,"]") )
    for (h in 2:H) lines( xi_h_m[j, k, h, from:m ], type="l", col=h + 1)
}

# load data depending on format of dataFile
if ( is.character(Data$dataFile) ) {
    tp <- load(Data$dataFile)
    stopifnot( tp == "Njk.i" ) 
} else {
    Njk.i <- Data$dataFile
}

# load prior data (if required) depending on format of mccFile
if (Prior$mccAsPrior | Initial$mccUse) { 
    if ( is.character(Data$mccFile) ) {
        tp2 <- load(Data$mccFile)
        stopifnot( tp2 == "mcc" ) 
    } else {
        mcc <- Data$mccFile
    } 
} else { 
    mcc <- NULL 
} 

# report something...
cat(textTemp <- paste("workspaceFile: ", workspaceFile, "  (within current working directory!) \n") )
write(textTemp, file = logFileName, append = TRUE)
cat("Data loaded!", fill=TRUE)
flush.console()

# calculate N and K from data
N <- dim(Njk.i)[3]
K <- dim(Njk.i)[1] - 1

# Define dimnames...
dnj <- paste("j",1:(K+1),sep="")
dnk <- paste("k",1:(K+1),sep="")
dnh <- paste("h",1:H,sep="")
dnm <- paste("m",1:M,sep="")
dnn <- paste("n",1:N,sep="") 
dnhj<- paste("hj",1:(H*(K+1)),sep="")

# report something...
cat(textTemp <- paste("Data Information: Datafile =", if ( is.character(Data$dataFile) ) Data$dataFile else "no file name", 
                            ", N =", N, ", K =", K ), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
cat(textTemp <- paste("Manual Settings: No of groups H =", H, ", kNo =", kNo), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
cat(textTemp <- paste("MCMC Parameters: M =", M, ", M0 =", M0, ", mOut =", mOut, ", mSave =", mSave, ", seed =", seed, ", showAcc =", showAcc), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
textTemp <- if (isPriorNegBin) {
                    paste("Prior Parameters for e_h (Neg Binom): a0 =", a0, ", N0 =", N0, ", xi_prior (see below)")
            } else {
                    paste("Prior Parameters for e_h (Neg Multinom): a0 =", a0, ", alpha =", alpha, ", N0 =", N0, ", xi_prior (see below)")
            }
cat(textTemp, fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
flush.console()

# prepare prior parameter for xi
xi_prior <- 
    if (mccAsPrior) {    
        if (xiPooled) {            
            xiPriorText <- "estimated transition matrix (Bayesian) from overall pooling (mccAsPrior = TRUE & xiPooled = TRUE)"            
            array( mcc[[1]]$xi, c(K + 1, K + 1, H) )  # eg. (Bayesian) estimated transition matrix from overall pooling  OR empirical trans mat  
        } else {    
            xiPriorText <- "estimated transition matrix (Bayesian) in each group (mccAsPrior = TRUE & xiPooled = FALSE)"            
            array( mcc[[H]]$xi, c(K + 1, K + 1, H) )  
        }    
    } else {    
        xiPriorText <- paste( "with persPrior = ", round(persPrior, 4), " created xi_prior (equal in each group, mccAsPrior = FALSE & xiPooled = FALSE)" )        
        transPr <- (1 - persPrior)/K         # off-diagonal    
        xiTemp <- (persPrior - transPr)*diag(1, K + 1, K + 1) + matrix(transPr, K + 1, K + 1, dimnames=list(dnj, dnk)) # prior par (for Neg Multinomial)
        array( xiTemp, c(K + 1, K + 1, H) ) 
    }

dimnames(xi_prior) <- list(dnj, dnk, dnh)

cat(textTemp <- paste("Information on xi_prior (for Neg Bin/Neg Multinom Prior): ", xiPriorText ), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
cat("Prior information and parameters set!", fill=TRUE)
flush.console()

# set prior pars for beta (bk0, Bk0inv)
if (Prior$betaPrior == "uninformative") {
    bk0 <-   matrix(0, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"0",sep="") ) )
    Bk0inv <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"0inv",sep="")))
} else {
    bk0 <-   matrix(Prior$betaPriorMean, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"0",sep="") ) )
    Bk0inv <- array( chol2inv( chol( diag(Prior$betaPriorVar, XNo) ) ), c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"0inv",sep="")))
}

# start value for Beta
Beta.start <- if ( is.null( Initial$Beta.start ) ) matrix( 0, XNo, H, dimnames=list( dimnamesX , dnh )) else Initial$Beta.start
stopifnot( dim(Beta.start) == c(XNo, H) )

# prepare initial values for e_h
if (mccUse) {     
    xi_h_start <- array( mcc[[H]]$xi, c(K + 1, K + 1, H) ) 
    dimnames(xi_h_start) <- list(dnj, dnk, dnh)
    e_start <- N0*xi_h_start   
} else { 
    trans <- (1 - pers)/K  # off-diagonal
    xi_h_start <- (pers - trans)*diag(1, K + 1, K + 1) + matrix(trans, K + 1, K + 1, dimnames=list(dnj, dnk))
    e_start <- N0*array(xi_h_start, c(K +1, K +1, H), dimnames=list(dnj, dnk, dnh)) 
}
e_h_0 <- ceiling(e_start) #  has to be discrete
dimnames(e_h_0) <- list(dnj, dnk, dnh)

# report something...
extra <- if (mccUse) "" else paste(", pers =", round(pers, 5) )
cat(textTemp <- paste("Inital Values Information: mccUse =", mccUse, extra), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
cat("Initial values set!", fill=TRUE)
flush.console()

# initialisations 
e_h_m    <- array(0, c(K + 1, K + 1, H, M), dimnames=list( dnj, dnk, dnh, dnm ) )
S_i_freq <- matrix(0, H, N, dimnames=list(dnh, dnn ))
accept   <- array(0, c(M, H*(K+1), 2), dimnames = list(dnm, dnhj, c("accProb", "accYesNo") ) )
xi_h_m   <- array(0, c(K + 1, K + 1, H, M), dimnames=list( dnj, dnk, dnh, dnm ) )

# initialise posterior pars for beta
bkN <- matrix(0, XNo, H, dimnames=list(dimnamesX, paste("b",1:H,"N",sep="") ) )
BkNinv <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"Ninv",sep="")))
BkN    <- array(0, c(XNo, XNo, H), dimnames=list(dimnamesX, dimnamesX, paste("B",1:H,"N",sep="")))

# initialise (unknown) parameter matrix
Beta.m <- array( 0, c(XNo, H, M), dimnames=list( dimnamesX , dnh, NULL ))

cat("Initialisations done!", fill=TRUE)
cat("MCMC Iteration...", fill=TRUE)
flush.console()

# MCMC-sampler 
# ============ 

# Bayes' classification for each subject i 
# -> sample S_i from a discrete probability distribution p(S_i | y_i, eta, e_1, ..., e_H )

# indicator matrix for group membership (Si != h) !!!
Ind <- matrix(1, N, H, dimnames=list( NULL , dnh ) )

# set indicator matrix  ( 10000 executions take approx 600 seconds = 10 minutes !! )
for (i in 1:N) Ind[ i, Initial$S.i.start[i] ] <- 0  
IndCut <- Ind[,2:H]

# new aux mix: sample uniform random numbers -- used to sample utilities
Unif <- matrix( runif(N*(H-1)), N, H-1, dimnames=list( NULL , dnh[-1] ) )

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

SRmat  <- matrix( sr[Rih], N, H - 1, dimnames=list(NULL, dnh[-1]) )
SR2mat <- matrix(sr2[Rih], N, H - 1, dimnames=list(NULL, dnh[-1]) )

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

    Beta.m[,h,1] <- mvrnorm( n=1, mu=bkN[,h], Sigma=BkN[,,h] )
}

#------------------------

## Pr(Si=h|betah,h=1,...,H)
XBetak <- crossprod(t(X), Beta.m[,,1]) 
logit.start <- exp( XBetak ) / rowSums( exp( XBetak ) ) # N x H -matrix
logit.temp <- logit.start

# calculate marginal likelihood
margDataLikeli <- matrix(0, N, H)
e_h <- e_h_0

if (H > 1) {
    AA <- apply( lgamma(apply(e_h, 3, rowSums)), 2, sum )   
    BB <- apply( lgamma(e_h), 3, sum)                       
    for ( i in 1:N) {
        for ( h in 1:H ) {
            margDataLikeli[i, h] <- margDataLL_0(i=i, h=h, Njk.i=Njk.i, e_h=e_h, AA=AA, BB=BB) # TO REVISE (WITH MAPPLY?)...
        }
    }    
    likeli <- margDataLikeli*logit.start  # (N x H)-matrix
    S_i_m_temp <- matrix(0, H, N, dimnames=list(dnh, dnn ))    
    for (i in 1:N) S_i_m_temp[,i] <- rmultinom(n=1, size=1, prob=likeli[i,] ) # is internally normalized to sum 1
} else {    
    S_i_m_temp <- matrix(1, H, N, dimnames=list(dnh, dnn ))     
}

# distribution of group sizes 
N_h <- rowSums(S_i_m_temp) 

# distribution of group parameters 
# first update
for ( h in 1:H ) {    
    S_i_h_ind <- which( as.logical( S_i_m_temp[h,] ) )    # TO REVISE...
    N_h_temp <- N_h[h]         
    for ( j in 1:(K+1) ) {        
        upDate <- ehjUpdate_k(ehj=e_h_0[j,,h], h=h, j=j, kNo=kNo, m=1, isPriorNegBin=isPriorNegBin, K=K, N_h_temp=N_h_temp, S_i_h_ind=S_i_h_ind, 
                              Njk.i=Njk.i, a0=a0, N0=N0, xi_prior=xi_prior, alpha=alpha)        
        e_h_m[j,,h,1] <- upDate$e_hj        
        accept[1, j + (h-1)*(K+1), ] <- upDate$accept            
    }
    xi_h_m[,,h,1] <- e_h_m[,,h,1]/rowSums(e_h_m[,,h,1])
}

# report something...
cat(textTemp <- paste("m = 1", if (showAcc) paste("; Acc Rate of first draws =", round( mean( accept[1, , 2] ), 2)  )), fill=TRUE)
write(textTemp, file = logFileName, append = TRUE)
flush.console()

# if necessary monitor draws...
if (monitor) {
    dev.new()
    if (K <= 2) par(mfrow=c(1, 3)) else par(mfrow=c(2, 3))
    e_dev <- dev.cur()
    dev.new()
    if (K <= 2) par(mfrow=c(1, 3)) else par(mfrow=c(2, 3))
    xi_dev <- dev.cur() 
}

# calculate prior densities
for ( h in 1:H) {   
    for ( j in 1:(K+1) ) {            
            logEPrior[1] <- logEPrior[1] +  
                            log( priorDens( xVec=e_h_m[j, ,h,1] - 1, k=a0, alpha=alpha, lambdaVec=N0*xi_prior[j,,h], 
                                            isPriorNegBin=isPriorNegBin ) )
    }            
}

ptm <- proc.time()[3]
for ( m in 2:M ) {

    # indicator matrix for group membership (Si != h) !!!
    Ind <- matrix(1, N, H, dimnames=list( NULL , dnh ) )

    # set indicator matrix  ( 10000 executions take approx 600 seconds = 10 minutes !! )
    S.i.temp <- apply(S_i_m_temp, 2, which.max)
    for (i in 1:N) Ind[ i, S.i.temp[i] ] <- 0  
    IndCut <- Ind[,2:H]    

    # new aux mix: sample uniform random numbers -- used to sample utilities
    Unif <- matrix( runif(N*(H-1)), N, H-1, dimnames=list( NULL , dnh[-1] ) )

    # calculate XBeta and Lambda
    XBeta  <- XBetak # from previous iteration   
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

    SRmat  <- matrix( sr[Rih], N, H - 1, dimnames=list(NULL, dnh[-1]) )
    SR2mat <- matrix(sr2[Rih], N, H - 1, dimnames=list(NULL, dnh[-1]) )

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
    logit.prev <- logit.temp  # from  'm-1'
    
    XBetak <- crossprod(t(X), Beta.m[,,m])   
    logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) # N x H -matrix


    # update classifications 
    if (H > 1) {  
        AA <- apply( lgamma(apply(e_h_m[,,,m-1], 3, rowSums)), 2, sum )  
        BB <- apply( lgamma(e_h_m[,,,m-1]), 3, sum )
        for ( i in 1:N) {
            for ( h in 1:H ) {
                margDataLikeli[i, h] <- margDataLL_m(i=i, m=m-1, h=h, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB)
            }
        }  
        likeli <- margDataLikeli*logit.temp # margDataLikeli from previous iteration & logit.temp from current iteration   #  updated (N x H)-Matrix   
        for (i in 1:N) S_i_m_temp[,i] <- rmultinom(n=1, size=1, prob=likeli[i,] ) # is internally normalized to sum 1
    } else {
        S_i_m_temp <- matrix(1, H, N)         
    }     
    if (m > M0) S_i_freq <- S_i_freq + S_i_m_temp  # M-M0 times

    # store log likelihood: 
    # ===================== 

    # observed log likelihood
    tempLK  <- margDataLikeli * logit.prev    
    tik <- tempLK/rowSums(tempLK)
    sProbsMax <- max.col(tempLK) 
    
    logLike[ m - 1 ] <- sum( log( rowSums( tempLK ) ) )   # "observed" log likelihood!!
    
    logClassLike[ m - 1 ] <- sum( log( ( tempLK )[ cbind( 1:N, sProbsMax ) ] ) )   # classification log likelihood!!

    entropy[ m - 1 ] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )
   
    logBetaPrior[ m - 1 ] <- sum( dnorm( Beta.m[ , -1, m - 1], mean = Prior$betaPriorMean, sd = sqrt(Prior$betaPriorVar), log=TRUE ) )

    # update group sizes 
    N_h_temp0 <- rowSums(S_i_m_temp)
        
    # update group parameters  
    for ( h in 1:H ) {        
        S_i_h_ind <- which( as.logical( S_i_m_temp[h,] ) ) # TO REVISE...
        N_h_temp <- N_h_temp0[h]         
        for ( j in 1:(K+1) ) {         
            upDate <- ehjUpdate_k(ehj=e_h_m[j,,h, m - 1 ], h=h, j=j, K=K, kNo=kNo, m=m, isPriorNegBin=isPriorNegBin, N_h_temp=N_h_temp, 
                                         S_i_h_ind=S_i_h_ind, Njk.i=Njk.i, a0=a0, N0=N0, xi_prior=xi_prior, alpha=alpha)            
            e_h_m[j,,h,m] <- upDate$e_hj             
            accept[m, j + (h-1)*(K+1), ] <- upDate$accept     
            
            # calculate prior densities
            logEPrior[m] <- logEPrior[m] +  
                            log( priorDens( xVec=e_h_m[j, ,h,m] - 1, k=a0, alpha=alpha, lambdaVec=N0*xi_prior[j,,h], 
                                            isPriorNegBin=isPriorNegBin ) )                    
        }    
        xi_h_m[,,h,m] <- e_h_m[,,h,m]/rowSums(e_h_m[,,h,m])  
    } 
    
    # report progress 
    if ( identical(all.equal(m %% mOut, 0), TRUE) || is.element(m, c(1:5, 10, 20, 50, 100, 200, 500) )  ) {  
            if (identical(all.equal(m %% mSave, 0), TRUE)) save(list = ls(), file = file.path(storeDir, paste(startDate, "_temp.RData", sep="")), compress = TRUE)                                                      
            if (monitor) {  
                dev.set(which = e_dev )
                sapply(1:(K+1), function(k) { monitor_e_h(j=k, k=k, e_h_m=e_h_m, from=1, m=m, M=M, H=H)
                                              #if (exists( "e_true" ) ) abline(h=e_true[k,k,], col=2:(H+1), lty=2) 
                                            } 
                )                
                dev.set(which = xi_dev )
                sapply(1:(K+1), function(k) { monitor_xi_h(j=k, k=k, xi_h_m=xi_h_m, from=1, m=m, M=M, H=H)
                                              #if (exists( "xi_true" ) ) abline(h=xi_true[k,k,], col=2:(H+1), lty=2) 
                                            } 
                )                
            }            
            cat(textTemp <- paste("m =", m, "; duration of iter proc so far:", round(diff <- proc.time()[3] - ptm,2), "sec. ,  exp time to end:", 
                round( (diff/(m-1)*M - diff)/60, 2 ), " min.",  
                if (showAcc) paste("; Acc Rate of last", min(m, mOut), "draws =", round( mean( accept[max(1, m - mOut):m, , 2] ), 2)  )),
                fill=TRUE)
            flush.console()
            write(textTemp, file = logFileName, append = TRUE)
            }     
}

# store log likelihood: 
# ===================== 

if (H > 1) {  
    AA <- apply( lgamma(apply(e_h_m[,,,m], 3, rowSums)), 2, sum )  
    BB <- apply( lgamma(e_h_m[,,,m]), 3, sum )
    for ( i in 1:N) {
        for ( h in 1:H ) {
            margDataLikeli[i, h] <- margDataLL_m(i=i, m=m, h=h, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB)
        }
    }  
} 

# observed log likelihood
tempLK  <- margDataLikeli * logit.temp  
tik <- tempLK/rowSums(tempLK)
sProbsMax <- max.col(tempLK) 

logLike[ m ] <- sum( log( rowSums( tempLK ) ) )   # "observed" log likelihood!!

logClassLike[ m ] <- sum( log( ( tempLK )[ cbind( 1:N, sProbsMax ) ] ) )    # classification log likelihood!!

entropy[ m ] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )

logBetaPrior[ m ] <- sum( dnorm( Beta.m[ , -1, m], mean = Prior$betaPriorMean, sd = sqrt(Prior$betaPriorVar), log=TRUE ) )

logPostDens <- logLike + logEPrior  + if ( Prior$betaPrior!="uninformative" ) logBetaPrior else 0 
mMax <- which.max(logPostDens)

# report and write out totel time
cat( totalTimeTemp <- paste("Total time: ", (totTime <- proc.time()[3] - pts)%/%3600, "hours", (totTime%%3600)%/%60, "min" ), "\n" )
write( totalTimeTemp , file = logFileName, append = TRUE)

# save selection of final workspace
save(list = c("accept", "Beta.m", "bk0", "Bk0inv", "Data", "e_h_0", "e_h_m", "fileName", "Initial", "K", "logFileName", "mcc", "Mcmc", 
              "N", "Njk.i", "Prior", "S_i_freq", "xi_h_m", "xi_prior", "bkN", "BkN", "workspaceFile", 
              "logLike", "logBetaPrior", "logEPrior", "logPostDens", "mMax", "logClassLike", "entropy"), 
              file = workspaceFile, compress = TRUE) # 

# delete temporary buffer store 'temp.RData'
unlink( file.path(storeDir, paste(startDate, "_temp.RData", sep="")) )   

# rename log file with relevant information
file.rename(logFileName, file.path( storeDir, fileNameTxt) )

# build list to be returned
resultsList <- list( workspaceFile=workspaceFile, accept=accept, Beta.m=Beta.m, bk0=bk0, Bk0inv=Bk0inv, Data=Data, e_h_0=e_h_0, e_h_m=e_h_m, fileName=fileName, Initial=Initial, K=K, 
                     logFileName=logFileName, mcc=mcc, Mcmc=Mcmc, N=N, Njk.i=Njk.i, Prior=Prior, S_i_freq=S_i_freq, xi_h_m=xi_h_m, xi_prior=xi_prior, bkN=bkN, BkN=BkN, logLike=logLike, 
                     logBetaPrior=logBetaPrior, logEPrior=logEPrior, logPostDens=logPostDens, mMax=mMax, logClassLike=logClassLike, entropy=entropy )

# return results file name
return( invisible( resultsList ) )

}
