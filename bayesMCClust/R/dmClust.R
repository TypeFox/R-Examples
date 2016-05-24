
dmClust <- function( Data = list(dataFile = stop("'dataFile' must be specified: filename or data"),
                                 storeDir = "try01", 
                                 mccFile = "mcc.RData"), 
                     Prior = list(H = 4, 
                                  alpha0 = 4, 
                                  a0 = 1, 
                                  alpha = 1, 
                                  N0 = 10,     
                                  isPriorNegBin = FALSE, 
                                  mccAsPrior = FALSE, 
                                  xiPooled = TRUE, 
                                  persPrior = 7/10),
                     Initial = list(mccUse = FALSE, 
                                    pers = 1/6, 
                                    S.i.start = NULL ),
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

alpha0        <- Prior$alpha0         # prior information for eta (group size)
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
fileName    <- paste("DMC_H", H, "_M", M, "_", startDate, ".RData", sep="")
fileNameTxt <- paste("DMC_H", H, "_M", M, "_", startDate, "_logfile.txt", sep="")
workspaceFile <- file.path( storeDir, fileName)

# define log file name
logFileNameTemp <- paste(startDate, "_logfile.txt", sep="")
logFileName <- file.path(storeDir, logFileNameTemp) # paste(storeDir, startDate, "_logfile.txt", sep="")

# store logLike and priors in each iteration
logLike <- numeric(M)
logClassLike <- numeric(M)
entropy <- numeric(M)
logEtaPrior <- numeric(M) 
logEPrior <- numeric(M) 

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
cat(textTemp <- paste("Prior Parameters for eta (Dirichlet): alpha0 =", alpha0), fill=TRUE)
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

# prepare initial values for eta
eta_start <- if (mccUse) {                     
                    mcc[[H]]$eta                   
             } else {             
                    rep(1, H)/H
             }

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
eta_m    <- matrix(0, M, H, dimnames=list(dnm, dnh ))
S_i_freq <- matrix(0, H, N, dimnames=list(dnh, dnn ))
accept   <- array(0, c(M, H*(K+1), 2), dimnames = list(dnm, dnhj, c("accProb", "accYesNo") ) )
xi_h_m   <- array(0, c(K + 1, K + 1, H, M), dimnames=list( dnj, dnk, dnh, dnm ) )

cat("Initialisations done!", fill=TRUE)
cat("MCMC Iteration...", fill=TRUE)
flush.console()

# MCMC-sampler 
# ============ 

# Bayes' classification for each subject i 
# -> sample S_i from a discrete probability distribution p(S_i | y_i, eta, e_1, ..., e_H )

margDataLikeli <- matrix(0, H, N)
e_h <- e_h_0

if ( is.null( Initial$S.i.start ) ) {

    if (H > 1) {
        AA <- apply( lgamma(apply(e_h, 3, rowSums)), 2, sum )   
        BB <- apply( lgamma(e_h), 3, sum)                       
        for ( i in 1:N) {
            for ( h in 1:H ) {
                margDataLikeli[h, i] <- margDataLL_0(i=i, h=h, Njk.i=Njk.i, e_h=e_h, AA=AA, BB=BB)
            }
        }    
        likeli <- margDataLikeli*eta_start  # (H x N)-matrix
        S_i_m_temp <- matrix(0, H, N, dimnames=list(dnh, dnn ))    
        for (i in 1:N) S_i_m_temp[,i] <- rmultinom(n=1, size=1, prob=likeli[,i] ) # is internally normalized to sum 1
    } else {    
        S_i_m_temp <- matrix(1, H, N, dimnames=list(dnh, dnn ))     
    }

} else {
    
    if (H > 1) {
    
        S_i_m_temp <- matrix(0, H, N, dimnames=list(dnh, dnn )) 
        S_i_m_temp[cbind(Initial$S.i.start, 1:N)] <- 1
    
    } else {    
    
        S_i_m_temp <- matrix(1, H, N, dimnames=list(dnh, dnn )) 
    
    }

}

# distribution of group sizes 
N_h <- rowSums(S_i_m_temp) 
a_h <- N_h + alpha0
eta_1 <- gtools::rdirichlet(n=1, alpha=a_h )  # produces row-vector (1 x H) !!
eta_m[1,] <- eta_1

# distribution of group parameters 
# first update
for ( h in 1:H ) {    
    S_i_h_ind <- which( as.logical( S_i_m_temp[h,] ) )    
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
cat(textTemp <- paste("workspaceFile: ", workspaceFile, "  (within current working directory!) \n") )
write(textTemp, file = logFileName, append = TRUE)
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
logEtaPrior[1] <- log( ddirichlet( eta_m[1,], rep(alpha0,H)) )
for ( h in 1:H) {   
    for ( j in 1:(K+1) ) {            
            logEPrior[1] <- logEPrior[1] +  
                            log( priorDens( xVec=e_h_m[j, ,h,1] - 1, k=a0, alpha=alpha, lambdaVec=N0*xi_prior[j,,h], 
                                            isPriorNegBin=isPriorNegBin ) )
    }            
}

ptm <- proc.time()[3]
for ( m in 2:M ) {  
    
    # update classifications 
    if (H > 1) {  
        AA <- apply( lgamma(apply(e_h_m[,,,m-1], 3, rowSums)), 2, sum )  
        BB <- apply( lgamma(e_h_m[,,,m-1]), 3, sum )
        for ( i in 1:N) {
            for ( h in 1:H ) {
                margDataLikeli[h, i] <- margDataLL_m(i=i, m=m-1, h=h, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB)
            }
        }  
        likeli <- margDataLikeli*eta_m[m-1,] # updated (H x N)-Matrix   
        for (i in 1:N) S_i_m_temp[,i] <- rmultinom(n=1, size=1, prob=likeli[,i] ) # is internally normalized to sum 1
    } else {
        S_i_m_temp <- matrix(1, H, N)         
    }     
    if (m > M0) S_i_freq <- S_i_freq + S_i_m_temp  # M-M0 times

    # store log likelihood: 
    # ===================== 

    # observed log likelihood
    if (H == 1) {      
        AA <- apply( lgamma(apply(array(e_h_m[,,,m-1], c(K+1, K+1, H)), 3, rowSums)), 2, sum )     # 
        BB <- apply( lgamma(array(e_h_m[,,,m-1], c(K+1, K+1, H))), 3, sum )
        for ( i in 1:N) {
            for ( h in 1:H ) {
                margDataLikeli[h, i] <- margDataLL_m(i=i, m=m-1, h=h, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB)
            }
        }  
        likeli <- margDataLikeli*eta_m[m-1,] # updated (H x N)-Matrix   
    }
    
    sProbs <- t(likeli)
    tik <- sProbs/rowSums(sProbs)
    sProbsMax <- max.col(sProbs) 
    
    # "observed" log likelihood!!
    logLike[ m-1 ] <- sum( log( colSums( likeli ) ) )
    
    logClassLike[ m-1 ] <- sum( log( likeli[cbind(sProbsMax, 1:N)] ) )          
            
    entropy[ m-1 ] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )            

    # update group sizes 
    N_h_temp0 <- rowSums(S_i_m_temp)
    a_h <- N_h_temp0 + alpha0 
    eta_m[m,] <- if (H > 1) gtools::rdirichlet(n=1, alpha=a_h ) else 1  # produces row-Vector (1 x H) !!
        
    # update group parameters  
    for ( h in 1:H ) {        
        S_i_h_ind <- which( as.logical( S_i_m_temp[h,] ) ) 
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
    
    # calculate prior densities
    logEtaPrior[m] <- log( ddirichlet( eta_m[m,], rep(alpha0,H)) )
    
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

# update classifications 
AA <- apply( lgamma(apply(array(e_h_m[,,,m], c(K+1, K+1, H)), 3, rowSums)), 2, sum )     # 
BB <- apply( lgamma(array(e_h_m[,,,m], c(K+1, K+1, H))), 3, sum )
for ( i in 1:N) {
    for ( h in 1:H ) {
        margDataLikeli[h, i] <- margDataLL_m(i=i, m=m, h=h, Njk.i=Njk.i, e_h_m=e_h_m, AA=AA, BB=BB)
    }
}  

likeli <- margDataLikeli*eta_m[m,] # updated (H x N)-Matrix   
sProbs <- t(likeli)  
tik <- sProbs/rowSums(sProbs)

sProbsMax <- max.col(sProbs)

# store log likelihood: 
# ===================== 

# observed log likelihood
logLike[ m ] <- sum( log( colSums( likeli ) ) ) # "observed" log likelihood!!

logClassLike[ m ] <- sum( log( likeli[cbind(sProbsMax, 1:N)] ) )  

entropy[ m ] <- if ( min(tik) > 1e-320 ) -sum( tik*log(tik) )         
        
logPostDens <- logLike + logEtaPrior + logEPrior
mMax <- which.max(logPostDens)

# report and write out totel time
cat( totalTimeTemp <- paste("Total time: ", (totTime <- proc.time()[3] - pts)%/%3600, "hours", (totTime%%3600)%/%60, "min" ), "\n" )
write( totalTimeTemp , file = logFileName, append = TRUE)

# save selection of final workspace
save(list = c("accept", "Data", "e_h_0", "e_h_m", "eta_m", "fileName", "Initial", "K", "logFileName", "mcc", "Mcmc", "N", "Njk.i", "Prior", 
              "S_i_freq", "xi_h_m", "xi_prior", "logLike", "logClassLike", "entropy", "logEtaPrior", "logEPrior", "logPostDens", "mMax", "workspaceFile"), 
              file = workspaceFile) # 

# delete temporary buffer store 'temp.RData'
unlink( file.path(storeDir, paste(startDate, "_temp.RData", sep="")) )   

# rename log file with relevant information
file.rename(logFileName, file.path( storeDir, fileNameTxt) )

# build list to be returned
resultsList <- list( workspaceFile=workspaceFile, accept=accept, Data=Data, e_h_0=e_h_0, e_h_m=e_h_m, eta_m=eta_m, fileName=fileName, Initial=Initial, K=K, logFileName=logFileName,
                     mcc=mcc, Mcmc=Mcmc, N=N, Njk.i=Njk.i, Prior=Prior, S_i_freq=S_i_freq, xi_h_m=xi_h_m, xi_prior=xi_prior, logLike=logLike, logClassLike=logClassLike, 
                     entropy=entropy, logEtaPrior=logEtaPrior, logEPrior=logEPrior, logPostDens=logPostDens, mMax=mMax )

# return results file name
return( invisible( resultsList ) )

}
