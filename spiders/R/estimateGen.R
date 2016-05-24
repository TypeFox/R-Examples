##' @title estimate a reparameterization of the  hypothesis c_st
##'
##' @description estimates parameters from hypothesis lambda != gamma, where lambda
##' is indepdent of gamma
##'
##' @details There are 2*S*T free parameters under this hypothesis.
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period    
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estGen <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    ## some numbers
    S <- ncol(Ydst); s <- seq_len(S)
    T <- nrow(Ydst); t <- seq_len(T)
    ST <- S*T; st <- seq_len(ST)

    if (EM) {
        
        ## estimate gamma
        gammaHat <- Ydst/I

        ## initialize some things
        em_iter <- 1
        lambdaHat <- lambdaHat_old <- gammaHat

        ## iterate EM
        while ( TRUE ) {
            ## update lambda
            elambda <- exp(-lambdaHat)
            EX <- lambdaHat/(1-elambda)
            lambdaHat <- Xdst*EX/J

            ## check convergence
            if ( converged(lambdaHat, lambdaHat_old) ) {
                break
            }
            lambdaHat_old <- lambdaHat
            em_iter <- em_iter+1
            
            ## limit iterations
            if ( em_iter > em_maxiter ) {
                stop(sprintf("estGen: max EM iterations, %d, reached. Please adjust accordingly.", em_maxiter))
            }
        }
        
        ## since estimated 0s don't follow from calculations above
        lambdaHat[which(Xdst == 0, arr.ind=T)] <- 0

        ## calc standard error with est params
        ## SE <- seEM(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)
        Info <- diag(2*ST)             # initialize information matrix
        expl <- exp(lambdaHat)
        
        ## fill in second derivatives
        diag(Info)[-st] <- unlist(Ydst/gammaHat^2) # gamma
        diag(Info)[st] <- unlist(J/(expl-1)^2)     # lambda
        
        ## calc loglik wit est params
        loglik <- llEM(Xdst, Ydst, lambdaHat, gammaHat, J, I)

        list(lambda=as.matrix(lambdaHat), gamma=as.matrix(gammaHat),
             em_iters=em_iter, ll=loglik, var=solve(Info))

    } else {
        ## some numbers
        ## not sure this is the right spot for these checks
        T <- nrow(Xdst)
        if ( length(J) != T ) {
            stop("J indexed oddly says estGen")
        }
        if ( length(I) != T ) {
            stop("I indexed oddly says estGen")
        }

        gammaHat <- Ydst / I                # row-wise division
        lambdaHat <- Xdst / J

        ## calc standard error with est params
        ## SE <- se(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)
        Info <- diag(2*ST)             # initialize information matrix

        ## fill Info with second derivatives
        diag(Info)[-st] <- unlist(Ydst/gammaHat^2) # gamma
        diag(Info)[st] <- unlist(Xdst/lambdaHat^2) # lambda

        ## calc loglik with est params
        loglik <- ll(Xdst, Ydst, lambdaHat, gammaHat, J, I)
        
        list(lambda=as.matrix(lambdaHat), gamma=as.matrix(gammaHat),
             ll=loglik, var=solve(Info))
    }

}
