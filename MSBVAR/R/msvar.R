### msvar.R -- Estimates the MLE for MSVAR and MS univariate
### models.

# 20120109 : Initial version by Ryan Davis


msvar <- function(Y, p, h, niterblkopt=10)
{

# Switching indicator: 'IAH' totally switching
# fixed for now
    indms <- 'IAH'

    n <- nrow(Y)
    m <- ncol(Y)

    # check value of h
    if(h<2) stop("h should be an integer >=2")

# Now do a baseline, non-regime model using szbvar() since this
# gives us all of the inputs we need for later.
# n.b. the below specification is equivalent to MLE
#      (but if mu5 or mu6 does not equal zero, then problems)
    init.model <- szbvar(ts(Y), p,
                         lambda0=1, lambda1=1, lambda3=1, lambda4=1,
                         lambda5=1, mu5=0, mu6=0, prior=2)

# set initial parameters for blockwise optimization
# initial Q
    Qhat.start <- (1-(h*0.1/(h-1)))*diag(h) + matrix(0.1/(h-1), h, h)

# array for storage
    thetahat.start <- array(NA, c(m, 1+m*p+m, h))

# set intercept and AR coef initial values all to zero
    thetahat.start[,1:(1+m*p),] <- 0

# set sigma initial values
# first, get residuals from initial model
# dummy obs are appended, so adjust for those
    res.im <- init.model$residuals[(m+2):n,]
    sig2.start <- (1/n)*crossprod(res.im, res.im)

# the sig2 starting values need to be different though,
# so adjust these by a small amount over regimes
    for (i in 1:h) { thetahat.start[,(1+m*p+1):(1+m*p+m),i] <-
                         sig2.start}

    blkopt.est <- blkopt(Y=Y, p=p, thetahat.start=thetahat.start,
                         Qhat.start=Qhat.start, niter=niterblkopt,
                         indms)


# now, setup hreg, adjusting for dummies
    hreg <- hregime.reg2.mle(h, m, p, TT=(n-p), fp=blkopt.est$fpH, init.model)


    output <- list(init.model=init.model,
                   hreg=hreg,
                   Q=blkopt.est$Qhat,
                   fp=blkopt.est$fpH,
                   m=m, p=p, h=h,
                   llfval=blkopt.est$llfval,
                   DirectBFGSLastSuccess=blkopt.est$DirectBFGSLastSuccess)
    class(output) <- "MSVAR"

return(output)

} # end mlemsvar() function




# Ryan adjusted several items from the original hregime.reg2 function
hregime.reg2.mle <- function(h, m, p, TT, fp, init.model)
{

    # Storage
    tmp <- vector(mode="list", length=h)
    Bk <- array(0, c(m*p+1, m, h))
    Sigmak <- array(0, c(m,m,h))
    df <- apply(fp, 2, sum)
    e <- array(0, c(TT, m, h))
    Y <- init.model$Y[(m+1+1):nrow(init.model$Y),]
    X <- init.model$X[(m+1+1):nrow(init.model$X),]

    # Loops to compute
    # 1) sums of squares for X and Y
    # 2) B(k) matrices
    # 3) Residuals
    # 4) Sigma(k) matrices

    for(i in 1:h)
    {
        # Note how the dummy obs. are appended to the moment matrices
        Sxy <- crossprod(X, diag(fp[,i]))%*%Y + crossprod(init.model$X[1:(m+1),], init.model$Y[1:(m+1),])
        Sxx <- crossprod(X, diag(fp[,i]))%*%X + crossprod(init.model$X[1:(m+1),])

        # Compute the regression coefficients
        hstar <- Sxx# + init.model$H0
        #Bk[,,i] <- solve(hstar,
       #                  (Sxy + init.model$H0[,1:m]))
        Bk[,,i] <- solve(hstar,Sxy,tol=1e-100)

        # Compute residuals and Sigma (based on Krolzig)

        # Get the full residuals -- need these for filtering
        e[,,i] <- Y - X%*%Bk[,,i]

        #Sigmak[,,i] <- (init.model$S0 + crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]
        Sigmak[,,i] <- (crossprod(e[,,i],diag(fp[,i]))%*%e[,,i])/df[i]

        # Save the moments
        tmp[[i]] <- list(Sxy=Sxy, Sxx=Sxx) #, ytmp=ytmp, xtmp=xtmp)
    }

    return(list(Bk=Bk, Sigmak=Sigmak, df=df, e=e, moment=tmp))
}

