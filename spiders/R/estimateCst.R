##' @title estimate hypothesis c_st
##'
##' @description estimates parameters from hypothesis lambda = c_st * gamma
##'
##' @details There are 2*S*T free parameters under this hypothesis
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period    
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estCst <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    ## some numbers
    S <- ncol(Ydst); s <- seq_len(S)
    T <- nrow(Ydst); t <- seq_len(T)
    ST <- S*T; st <- seq_len(ST)

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
        ## cHat <- cHat_old <- matrix(0.5, ncol=S, nrow=T) # matrix(runif(ST), nrow=T)
        init <- est1(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma
        cHat <- cHat_old <- Xdst/(J*gammaHat) + pen
        
        ## iterate EM
        while ( TRUE ) {
            ## update lambda
            lambda <- cHat*gammaHat
            elambda <- exp(-lambda)
            EX <- lambda/(1-elambda)

            ## convenience
            ZEX <- Xdst*EX + pen      # pull ZEX away from zero

            ## iterate only once simultaneous eqs
            gammaHat <-  (ZEX + Ydst)/(J*cHat + I)                 
            cHat <- ZEX/(J*gammaHat)

            ## check convergence
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
                stop(sprintf("estCst: max EM iterations, %d, reached. Please adjust accordingly.", em_maxiter))
            }
                
        }
        
        ## since estimated 0s don't follow from calculations above
        ## lambdaHat[which(Xdst == 0, arr.ind=T)] <- 0

        ## calc standard error with est params
        ## SE <- seEM(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)
        Info <- diag(2*ST)             # initialize information matrix
        g2 <- gammaHat^2               # gamma^2
        l <- cHat*gammaHat
        expl <- exp(l)
        tmp <- J*expl/(expl-1)^2 # a common term
        
        ## fill in second derivatives
        diag(Info)[-st] <- unlist(Ydst/g2 + cHat^2*tmp) # gamma
        diag(Info)[st] <- unlist(g2*tmp)                # c
        Info[st,-st] <- unlist(-J/(expl-1) + l*tmp)     # fill off diags; upper tri
        lowmat <- lower.tri(Info)
        Info[lowmat] <- t(Info)[lowmat] # make symmetric from upper tri
        tryCatch(var <- solve(Info),
                 error=function(e) {
                     print("Variances not calculated; system is singular.")
                     var <- NULL
                 })

        ## calc loglik wit est params
        cv <- unlist(cHat)
        loglik <- llEM(Xdst, Ydst, NA, gammaHat, J, I, cv)
        
        list('gamma' = as.matrix(gammaHat), 'c' = cv,
             'em_iters' = em_iter, 'll' = loglik, 'var' = var)

    } else {
        ## some numbers
        XYdst <- Xdst + Ydst
        iter <- 1; maxiter <- 500
        
        ## not sure this is the right spot for these checks
        T <- nrow(Xdst)
        if ( length(J) != T ) {
            stop("J indexed oddly says estGen")
        }
        if ( length(I) != T ) {
            stop("I indexed oddly says estGen")
        }

        ## initialize some values
        ## cHat <- cHat_old <- matrix(0.5, ncol=S, nrow=T) # matrix(runif(ST), nrow=T)
        ## gammaHat <- gammaHat_old <- XYdst / (J*cHat + I)
        gammaHat <- gammaHat_old <- XYdst / (J + I)

        ## iteratively update; relies on concavity of log-lik
        while ( TRUE ) {
            
            ## update parameters
            cHat <- Xdst / (J*gammaHat)
            gammaHat <- XYdst / (J*cHat + I)

            ## check convergence
            if ( converged(gammaHat, gammaHat_old) &&
                converged(cHat, cHat_old) ) {
                break
            }

            ## if not converged, update estimates for next iteration
            gammaHat_old <- gammaHat
            cHat_old <- cHat
            iter <- iter+1

            ## limit iterations
            if ( iter > maxiter ) {
                stop(sprintf("estCst: %d not sufficient iterations for simultaneous equations.", maxiter))
            }
            
        }
        
        ## calc standard error with est params
        ## SE <- se(lambdaHat, gammaHat, NULL, Xdst, Ydst, J, I)
        Info <- diag(2*ST)             # initialize information matrix

        ## fill Info with second derivatives
        diag(Info)[-st] <- unlist(XYdst/gammaHat^2) # gamma
        diag(Info)[st] <- unlist(Xdst/cHat)         # c
        for ( i in t ) {
            Info[st,-st][,seq(0, (S-1)*T, by=T)+i] <- J[i]
        }
        lowmat <- lower.tri(Info)
        Info[lowmat] <- t(Info)[lowmat]       # make symmetric from upper tri

        ## calc loglik with est params
        cv <- unlist(cHat)
        loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cv)
        
        list(gamma=as.matrix(gammaHat), c=cv, iters=iter,
             ll=loglik, var=solve(Info))
    }

}
