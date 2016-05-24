##' @title estimate hypothesis c_st = c_s
##'
##' @description estimates parameters from hypothesis lambda_s = c_s * gamma_s
##'
##' @details There are S*T + S free parameters under this hypothesis.
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estCs <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    ## some numbers
    S <- ncol(Xdst); s <- seq_len(S)
    T <- nrow(Xdst); t <- seq_len(T)
    ST <- S*T

    if (EM) {
        ## avoiding estimated zeros
        pen <- matrix(0, T, S)
        zos <- which(Xdst==0, arr.ind=TRUE)
        if (length(zos)>0) {
            pen[zos] <- 1
            pen[zos] <- pen[zos]/(J[zos[,1]]+1) # add imaginary observation
            pen <- 0.001*pen                    # weight it
        }
        
        ## initialize some values
        em_iter <- 1
        init <- est1(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma
        cHat <- cHat_old <- sumT(Xdst) / sumT(J*gammaHat) + sumT(pen)


        ## iterate EM
        while ( TRUE ) {
            
            ## expected value of Xjst
            lambda <- sapply(s, function(j) cHat[j]*gammaHat[,j]) # col-wise `*`
            elambda <- exp(-lambda)
            EX <- lambda/(1-elambda)

            ## convenience
            ZEX <- Xdst*EX + pen        # pull ZEX away from zero

            ## iterate only once simultaneous eqs
            ZEXY <- ZEX + Ydst
            gammaHat <- sapply(s, function(j) ZEXY[,j] / (cHat[j]*J + I)) # col-wise `/`
            cHat <- sumT(ZEX) / sumT(J*gammaHat)
            
            ## check convergence of EM
            if ( converged(cHat, cHat_old) &&
                converged(gammaHat, gammaHat_old) ) {
                break
            }
            
            ## if not converged, store updated estimates
            cHat_old <- cHat
            gammaHat_old <- gammaHat
            ## print(sprintf('EM iteration %d found values: c = %f', em_iter, cHat))
            em_iter <- em_iter+1
            
            ## limit iterations
            if ( em_iter > em_maxiter ) {
                stop(sprintf("estCs: max EM iterations, %d, reached. Please adjust accordingly.", em_maxiter))
            }
            
        }

        ## calc standard error with est params
        ## SE <- seEM(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
        Info <- diag(ST+S)             # initialize information matrix
        g2 <- gammaHat^2               # gamma^2
        l <- sapply(s, function(j) gammaHat[,j]*cHat[j]) # col-wise `*`
        expl <- exp(l); tmp <- J*expl/(expl-1)^2  # a common term

        ## fill Info with second derivatives
        d2g2 <- sapply(s, function(j) tmp[,j]*cHat[j]^2 + Ydst[,j]/g2[,j]) # col-wise `*`
        diag(Info)[-s] <- unlist(d2g2)     # gamma
        diag(Info)[s] <- sumT(tmp*g2)      # c
        for ( i in s ) {           # fill in off diags; upper tri only
            Info[i,-s][t+T*(i-1)] <- unlist(-J/(expl[,i] -1) + l[,i]*tmp[,i])
        }
        lowmat <- lower.tri(Info)
        Info[lowmat] <- t(Info)[lowmat]      # make symmetric from upper tri
        tryCatch(var <- solve(Info),
                 error=function(e) {
                     print("Variances not calculated; system is singular.")
                     var <- NULL
                 })
        
        ## calc log-lik with est params
        loglik <- llEM(Xdst, Ydst, NA, gammaHat, J, I, cHat)
        
        list(c=cHat, gamma=as.matrix(gammaHat), em_iters=em_iter,
             ll=loglik, var=var)
    } else {

        XYdst <- Xdst + Ydst
        stXdst <- sumST(Xdst)
        iter <- 1; maxiter <- 500

        ## not sure this is the right spot for these checks
        ## ensure J & I have dimension T or 1
        if ( length(J) != T ) stop("J indexed oddly says est0.")
        if ( length(I) != T ) stop("I indexed oddly says est0.")

        ## initialize some values
        ## cHat <- cHat_old <- rep(0.5, S) # runif(S)
        ## gammaHat <- gammaHat_old <- sapply(s, function(j) XYdst[,j] / (J*cHat[j] + I))
        gammaHat <- gammaHat_old <- sapply(s, function(j) XYdst[,j] / (J + I))

        ## iteratively update; relies on concavity of log-lik
        while ( TRUE ) {

            ## update parameters
            cHat <- sumT(Xdst) / (sumT(J*gammaHat))
            gammaHat <- sapply(s, function(j) XYdst[,j] / (J*cHat[j] + I)) # col-wise `/`

            ## check convergence
            if ( converged(gammaHat, gammaHat_old) &&
                converged(cHat, cHat_old) ) {
                break
            }

            ## if not converged, update estimates for next iteration
            gammaHat_old <- gammaHat
            cHat_old <- cHat
            iter <- iter+1
            
            if ( iter > maxiter ) {
                stop(sprintf("estCt: %d not sufficient iterations for simultaneous equations.", maxiter))
            }
            
        }

        ## calc standard error with est params
        ## SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
        Info <- diag(ST+S)              # initialize information matrix

        ## fill Info with second derivatives
        diag(Info)[-s] <- unlist(XYdst/gammaHat^2) # gamma
        diag(Info)[s] <- sumT(Xdst)/cHat^2         # c
        for ( i in s ) {
            Info[i,-s][t+T*(i-1)] <- J        # fill off diags; upper tri only
        }
        lowmat <- lower.tri(Info)
        Info[lowmat] <- t(Info)[lowmat]      # make symmetric from upper tri
        

        ## calc log-lik with est params
        loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
            
        list(gamma=as.matrix(gammaHat), c=cHat, iters=iter,
             ll=loglik, var=solve(Info))
    }
}
