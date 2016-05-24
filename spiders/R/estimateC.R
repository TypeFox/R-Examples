##' @title estimate hypothesis c_st = c
##'
##' @description estimate parameters from hypothesis lambda = c*gamma
##'
##' @details There are S*T + 1 free parameters under this hypothesis.
##' 
##' @param Xdst matrix of sums of number of eaten prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param Ydst matrix sum of number of caught prey species s during occurrence t; rows indexed by time, and cols indexed by prey species, TxS
##' @param J vector of predators caught in each time period
##' @param I vector of number of days all traps were left out in a given time period
##' @param EM boolean; whether or not EM algorithm is used
##' @param em_maxiter integer specifying max number of EM iterations
##' @param BALANCED boolean; whether or not data are BALANCED
estC <- function(Xdst, Ydst, J, I, EM, em_maxiter, BALANCED) {

    ## some numbers
    S <- ncol(Xdst); s <- seq_len(S)
    T <- nrow(Xdst); t <- seq_len(T)
    ST <- S*T

    if (EM) {
        ## initialize some values
        em_iter <- 1
        init <- est1(Xdst, Ydst, J, I, FALSE, em_maxiter, BALANCED)
        gammaHat <- gammaHat_old <- init$gamma 
        lambda <- elambda <- init$lambda
        cHat <- cHat_old <- as.numeric(sumST(Xdst) / sumT(J*sumSp(gammaHat)))

        ## iterate EM
        while ( TRUE ) {
            
            ## expected value of Xjst
            lambda <- cHat*gammaHat
            elambda <- exp(-lambda)
            EX <- lambda / (1-elambda)

            ## convenience
            ZEX <- Xdst*EX

            ## iterate only once simultaneous eqs
            gammaHat <- (ZEX + Ydst) / (cHat*J + I)
            cHat <- sumST(ZEX) / sumT(J*sumSp(gammaHat))

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
                stop(sprintf("estC: max EM iterations, %d, reached. Please adjust accordingly.", em_maxiter))
            }
            
        }
        
        ## calc standard error with est params
        ## SE <- seEM(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
        Info <- diag(ST+1)                 # initialize information matrix
        g2 <- gammaHat^2                   # gamma^2
        l <- cHat*gammaHat
        expl <- exp(l); tmp <- J*expl/(expl-1)^2 # a common term

        ## fill I with second derivatives
        diag(Info)[-1] <- unlist(Ydst/g2 + cHat^2*tmp) # gamma
        diag(Info)[1] <- sumST(tmp*g2)                # c
        Info[1,-1] <- unlist(-J/(expl-1) + l*tmp) # fill in off diags; upper tri only
        lowmat <- lower.tri(Info)
        Info[lowmat] <- t(Info)[lowmat]      # make symmetric from upper tri
        tryCatch(var <- solve(Info),
                 error=function(e) {
                     print("Variances not calculated; system is singular.")
                     var <- NULL
                 })

        ## calc log-lik with est params
        loglik <- llEM(Xdst, Ydst, NA, gammaHat, J, I, as.numeric(cHat))
        
        list(c=cHat, gamma=as.matrix(gammaHat), em_iters=em_iter,
             ll=loglik, var=var)
    } else {
        if (BALANCED) {

            XYdst <- Xdst + Ydst
            iXdY <- I[1]*sumST(Xdst)/sumST(Ydst)
            gammaHat <- XYdst / (iXdY + I[1])
            cHat <- iXdY/J[1]

            ## calc standard error with est params
            ## SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
            Info <- diag(ST+1)                         # initialize information matirx

            ## fill second derivatives
            diag(Info)[-1] <- unlist(XYdst/gammaHat^2) # gamma
            Info[1,] <- sumST(Xdst/cHat^2)             # c
            Info[1,-1] <- J             # off diags; upper tri only
            lowmat <- lower.tri(Info)
            Info[lowmat] <- t(Info)[lowmat] # make symmetric from upper tri
            

            ## calc log-lik with est params
            loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
            
            list(gamma=as.matrix(gammaHat), c=cHat,
                 ll=loglik, var=solve(Info))
        } else {
            
            ## some numbers
            XYdst <- Xdst + Ydst
            stXdst <- sumST(Xdst)
            iter <- 1; maxiter <- 500

            ## not sure this is the right spot for these checks
            ## ensure J & I have dimension T or 1
            if ( length(J) != T ) {
                stop("J indexed oddly says est0.")
            }
            if ( length(I) != T ) {
                stop("I indexed oddly says est0.")
            }

            ## initialize some values
            ## cHat <- cHat_old <- 0.5 # runif(1)
            ## gammaHat <- gammaHat_old <- XYdst / (J*cHat + I)
            gammaHat <- gammaHat_old <- XYdst / (J + I)

            ## iteratively update; relies on concavity of log-lik
            while ( TRUE ) {

                ## update parameters
                cHat <- stXdst / sumT(J*sumSp(gammaHat))            
                gammaHat <- XYdst / (J*cHat + I) # row-wise division

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
                    stop(sprintf("estC: %d not sufficient iterations for simultaneous equations.", maxiter))
                }
            }
            
            ## calc standard error with est params
            ## SE <- se(NULL, gammaHat, cHat, Xdst, Ydst, J, I)
            Info <- diag(ST+1)                         # initialize information matirx

            ## fill second derivatives
            diag(Info)[-1] <- unlist(XYdst/gammaHat^2) # gamma
            Info[1,] <- sumST(Xdst/cHat^2)            # c
            Info[1,-1] <- J             # off diags; upper tri only
            lowmat <- lower.tri(Info)
            Info[lowmat] <- t(Info)[lowmat] # make symmetric from upper tri

            ## calc log-lik with est params
            loglik <- ll(Xdst, Ydst, NA, gammaHat, J, I, cHat)
            
            list(gamma=as.matrix(gammaHat), c=cHat, iters=iter,
                 ll=loglik, var=solve(Info))
        }
    }
}
