##' @title estimate hypothesis c_st = 1
##'
##' @description estimates parameters from hypothesis lambda = gamma
##'
##' @details There are S*T free parameters under this hypothesis.
##'
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
est1 <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    JI <- J + I
    if (EM) {
        em_iter <- 1
        init <- est1(Xdst, Ydst, J, I, FALSE, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma

        ## iterate EM
        while (EM) {
            
            ## some numbers
            lambda <- gammaHat
            elambda <- exp(lambda)
            EX <- lambda*elambda / (elambda - 1)
            ZEX <- Xdst*EX

            ## update parameters
            gammaHat <- (ZEX + Ydst) / JI

            ## check convergence of EM
            if ( converged(gammaHat, gammaHat_old) ) {
                break
            }
            
            ## if not converged, store updated estimates
            gammaHat_old <- gammaHat

            ## limit iterations
            em_iter <- em_iter+1
            if ( em_iter > em_maxiter ) {
                stop(sprintf("est1: max EM iterations, %d, reached. Please adjust accordingly.", em_maxiter))
            }
        }

        ## calc standard error with est params
        egamma <- exp(gammaHat)
        Info <- diag(unlist(Ydst/gammaHat^2 + J*egamma/(egamma - 1)^2))
        tryCatch(var <- solve(Info),
                 error=function(e) {
                     print("Variances not calculated; system is singular.")
                     var <- NULL
                 })
        
        ## calc loglik with est params
        loglik <- llEM(Xdst, Ydst, gammaHat, gammaHat, J, I)
        
        list(lambda=as.matrix(gammaHat), gamma=as.matrix(gammaHat),
             ll=loglik, var=var)
        
    } else {
        
        XYdst <- Xdst + Ydst
        gammaHat <-  XYdst / JI

        ## calc standard error with est params
        Info <- diag(unlist(XYdst/gammaHat^2))

        ## calc loglik with est params
        loglik <- ll(Xdst, Ydst, gammaHat, gammaHat, J, I)

        list(lambda=as.matrix(gammaHat), gamma=as.matrix(gammaHat),
             ll=loglik, var=solve(Info))        
    }
}
