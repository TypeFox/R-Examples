# initialize.msbvar.R -- Sets up initial values for the block
#                        optimization method used to find the mode for
#                        the MSVAR models
#
# 20120110 : Initial version -- PTB
# 20120504 : Updated to get initial regime parameters from a kmeans()
#            of the data
# 20120515 : Updated to use kmeans of VAR residuals to make more
#            robust
#

initialize.msbvar <- function(y, p, z=NULL, lambda0, lambda1, lambda3,
                              lambda4, lambda5, mu5, mu6, nu=NULL,
                              qm, prior, h, Q=NULL)

{
    # Set up the constants in the data matrices and the prior using a
    # call to szbvar()
    m <- ncol(y)
    if(is.null(nu)) nu <- m

    tmp <- szbvar(y, p, z=z, lambda0, lambda1, lambda3,
                  lambda4, lambda5, mu5, mu6, nu=m,
                  qm=qm, prior=prior,
                  posterior.fit=FALSE)

    # Define outputs for blkopt() and populate them
    thetahat.start <- array(NA, c(m, 1+m*p+m, h))

    # Use a kmeans of the data to separate across the regimes to
    # populate the regime-specific coefficients.

    #regimes <- kmeans(y, centers=h)$cluster
    regimes <- kmeans(tmp$residuals[(m+1):nrow(tmp$residuals),], centers=h)$cluster
    fptmp <- regimes[(p+1):length(regimes)]

    ################################################################
    # Now do regime-specific regressions for the data
    ################################################################

    # Parameters for regime-specific regressions
    df <- table(fptmp)
#    Sigma.draw <- array(NA, c(m,m,h))

    for(i in 1:h)
    {
        s <- c(rep(i,m+1), fptmp)
        Y1 <- tmp$Y[s==i,]
        X1 <- tmp$X[s==i,]
        XX <- crossprod(X1) + tmp$H0
        reg <- qr.coef(qr(crossprod(X1)),
                       crossprod(X1,Y1) + tmp$H0[,1:m])
        e <- Y1 - X1%*%reg

        # Sample the regime coefs
        # a) Inv-Wishart for Sigma for each regime
#        wisharts <- rwishart(h, df[i], diag(m))
        cSigma <- (tmp$S0 + tmp$H0[1:m,1:m] + crossprod(e))/(df[i]+nu)

        # Draw Sigma and tune it
#        Sigma.draw[,,i] <-
#            t(cSigma)%*%(df[i]*solve(wisharts[,,i]))%*%cSigma

        # Draw the regressors based on the tuned values
        ## bcoefs.se <- t(chol(chol2inv(chol(kronecker(cSigma, XX)))))


        ## bcoefs <- as.vector(reg) + crossprod(bcoefs.se,
        ##                                      rnorm((m^2)*p+m))

        ## # Reorder in MSBVAR format for VAR reg coefs
        ## Bhat <- matrix(bcoefs, ncol=m)

        # Populate the optimization object / array for blkopt

        # Intercepts
        thetahat.start[1:m,1,i] <- (reg[((m*p)+1),])

        # AR coefs
        thetahat.start[1:m, 2:((m*p)+1), i] <- t(reg[1:(m*p),])

        # Variances
        thetahat.start[1:m, (2+(m*p)):ncol(thetahat.start), i] <-
            t(cSigma)


    }

    # Now check the input Q or populate one.

    if(is.null(Q)==TRUE)
    {
        # Make Q from kmeans output
        lt <- length(fptmp)
        Qtmp <- table(fptmp[2:lt], fptmp[1:(lt-1)])
        Q <- Qtmp/rowSums(Qtmp)
    }

    # Now validate Q, user input

    if(is.null(Q)==FALSE)
    {

        if(sum(ifelse(rowSums(Q)==1, 1, 0))<h)
        {
            stop("initialize.msbvar(): Invalid user input: Improper initial Q, transition matrix, for the MS process.  Rows must sum to 1.\n")

        }
    }

    # Return the final object we need.
    return(list(init.model=tmp, thetahat.start=thetahat.start, Qhat.start=Q))
}
