## Function which fits an adaptive mixture of Student-t densities to a function KERNEL
## See Hoogerheide (2006, pp.46-49)
## __input__
## KERNEL     : [function] which computes the kernel !! must be vectorized !!
## mu0        : [kx1 vector] initial value for the location of the first component
##              (or location vector of the first component)
## Sigma0     : [kxk matrix] scaling matrix of the first component (default: NULL, i.e. estimated by 'AdMit')
## control    : [list] control parameters with the following components:
## $Ns        : [integer>100] number of draws used in the simulation (default: 1e5)
## $Np        : [integer>100] number of draws used when estimating the probabilities (default: 1e3)
## $Hmax      : [integer>0] maximum number of components (default: 10)
## $df        : [integer>0] degrees of freedom parameter (default: 1)
## $CVtol     : [double] relative decrease of the coefficient of variation (default: 0.1, i.e. 10%)
## $IS        : [logical] use importance sampling (default: FALSE)
## $ISpercent : [vector] of percentages of weights used to compute the scale matrix in IS sampling (default: c(.1, .15, .3))
## $ISscale   : [vector] of scaling coefficient for the covariance matrix in IS sampling (default: c(1,.25,4))
## $trace     : [logical] output printed during the fitting (default: FALSE)
## $trace.mu  : [logical] output printed of the optimizer (default: 0, i.e. no output)
## $maxit.mu  : [double] maximum number of iterations in the optimization (default: 1e4)
## $reltol.mu : [double] relative tolerance in the optimization (default: 1e-8)
## $trace.p   : [logical] output printed of the optimizer (default: 0, i.e. no output)
## $weightNC  : [double] weight of the new component (default: 0.1, i.e. 10%)
## $maxit.p   : [double] maximum number of iterations in the optimization (default: 1e4)
## $reltol.p  : [double] relative tolerance in the optimization (default: 1e-8)
## ...        : additional parameters used by 'KERNEL'
## __output__
## CV         : [Hx1 vector] of coefficient of variation
## mit        : [list] with the following components:
## $p         : [Hx1 vector] of probabilities
## $mu        : [Hxk matrix] of location vectors (in rows)
## $Sigma     : [Hxk^2 matrix] of scale matrices (in rows)
## $df        : [integer>0] degrees of freedom
## summary    : [data.frame] containing information on optimization algorithm, time and CV over fitting process
## __20081223__
'AdMit' <- function(KERNEL, mu0, Sigma0=NULL, control=list(), ...)
  {
    if (missing(KERNEL))
      stop ("'KERNEL' is missing in 'AdMit'")
    KERNEL <- match.fun(KERNEL)
    if (!any(names(formals(KERNEL))=="log"))
      stop ("'KERNEL' MUST have the logical argument 'log' which returns the (natural) logarithm of 'KERNEL' if 'log=TRUE'")  
    if (missing(mu0))
      stop ("'mu0' is missing in 'AdMit'")
    if (!is.vector(mu0))
      stop ("'mu0' must be a vector")
    if (!is.null(Sigma0))
      { ## if something is provided
        if (!is.matrix(Sigma0))
          { ## check if its a square matrix
            stop ("'Sigma0' must be a matrix")
          }
        else
          { ## if square matrix is provided
            if (!all(Sigma0==t(Sigma0))) ## check if the matrix is symmetric
              stop ("'Sigma0' is not symmetric")
            if (fn.isSingular(Sigma0)) ## check if the matrix is singular
              stop ("'Sigma0' is a singular matrix")
            if (!fn.isPD(Sigma0)) ## check is the matrix is positive definite
              stop ("'Sigma0' is not positive definite")
          }
      }
    Sigma0 <- as.vector(Sigma0) ## change square matrix into a vector
    
    con <- list(Ns=1e5, Np=1e3, Hmax=10, df=1, CVtol=.1,                           ## general control parameters
                IS=FALSE, ISpercent=c(.05,.15,.3), ISscale=c(1,.25,4),             ## importance sampling
                trace=FALSE,                                                       ## tracing information
                trace.mu=0, maxit.mu=5e2, reltol.mu=1e-8,                          ## mu optimization
                trace.p=0, weightNC=.1, maxit.p=5e2, reltol.p=1e-8)                ## p optimization
    con[names(control)] <- control
        
    if (con$Ns<100)
      stop ("'Ns' far too small.")
    if (con$Np<100)
      stop ("'Np' far too small")
    if (con$Np>con$Ns)
      stop ("'Np' must be lower or equal than 'Ns'")
    if (con$Hmax>15)
      warning ("'Hmax' larger than 15. May take some time and pose difficulties in the optimization")
    if (con$df<1)
      stop ("'df' must be greater or equal than 1")
    if (con$CVtol<=0 | con$CVtol>=1)
      stop ("'CVtol' must belong to ]0,1[")
    if (con$weightNC<=0 | con$weightNC>=1)
      stop ("'weightNC' must belong to ]0,1[")
    if (!is.logical(con$trace))
      stop ("'trace' must be logical")
    if (any(con$ISpercent<=0) | any(con$ISpercent>=1))
      stop ("components of 'ISpercent' must belong to ]0,1[")
    if (any(con$ISscale<0))
      stop ("components of 'ISscale' must be positive")

    controlIS <- list(percent=con$ISpercent, scale=con$ISscale)
    controloptmu <- list(trace=con$trace.mu, maxit=con$maxit.mu, reltol=con$reltol.mu)
    controloptp <- list(trace=con$trace.p, iter.max=con$maxit.p, rel.tol=con$reltol.p, weightNC=con$weightNC)
    
    fn.AdMit_sub(KERNEL, mu0, Sigma0, con$Ns, con$Np, con$Hmax, con$df, con$CVtol, con$trace,
                 con$IS, controlIS, controloptmu, controloptp, ...)
  }

'fn.AdMit_sub' <- function(KERNEL, mu0, Sigma0, Ns, Np, Hmax, df, CVtol, trace,
                           IS, controlIS, controloptmu, controloptp, ...)
  {
    ## initialisation 
    lnD <- lnK <- lnd <- h <- ph <- muh <- Sigmah <- NULL
    mit <- list(p=NULL, mu=NULL, Sigma=NULL, df=df)
    k <- length(mu0)
    Theta <- drawsh <- matrix(NA, Ns, k)
    nsummary <- c("H","METHOD.mu","TIME.mu","METHOD.p","TIME.p","CV")

    ## first component
    if (is.null(Sigma0)) 
      { ## optimize the first component
        optmu <- fn.optmu(FUN=KERNEL, mu0, control=controloptmu, ...)
        if (optmu$method=="IS")
          {
            stop ("Problem in the optimization. Try another starting values 'mu0'")
          }
      }
    else 
      { ## use the input by the user
        optmu <- list(mu=mu0, Sigma=Sigma0, method="USER", time=0)
      }

    ## first component
    h <- mit$p <- 1
    mit$mu <- matrix(optmu$mu, nrow=1)
    mit$Sigma <- matrix(optmu$Sigma, nrow=1)
    Theta <- drawsh <- fn.rmvt(Ns, mit$mu, mit$Sigma, df)
    lnK <- KERNEL(Theta, log=TRUE, ...)
    lnd <- fn.dmvt(Theta, mit$mu, mit$Sigma, df, log=TRUE)
    lnD <- as.matrix(lnd[1:Np])
    w <- fn.computeexpw(lnK, lnd)
    CV <- fn.CV(w)
    summary <- data.frame(h, optmu$method, optmu$time, "NONE", 0, CV)
    names(summary) <- nsummary
    if (trace)
      print(summary)
    
    hstop <- FALSE
    while (h<Hmax & hstop==FALSE)
      {
        h <- h+1
        if (!IS) ## usual optimization
          {
            mu00 <- drawsh[which.max(w),]
            tmp <- fn.muSigma(w, drawsh)
            theta <- fn.rmvt(Ns, tmp$mu, tmp$Sigma, 1)
            lnk <- KERNEL(theta, log=TRUE, ...)
            lnd <- fn.dmvt(theta, tmp$mu, tmp$Sigma, 1, log=TRUE)
            w1 <- fn.computeexpw(lnk, lnd)
            mu01 <- theta[which.max(w1),]
            mu0 <- rbind(mu00, mu01)
            tmpoptmu <- list()
            for (i in 1:2)
              { ## optimize using each starting value
                tmpoptmu[[i]] <- fn.optmu(FUN=fn.w, mu0[i,], control=controloptmu,
                                          KERNEL=KERNEL, mit=mit, ...)
              }
            tmptmpoptmu <- unlist(tmpoptmu)
            pos <- tmptmpoptmu[names(tmptmpoptmu)=="method"]!="IS"
            if (any(pos)) ## if at least one has converged
              {                
                v <- tmptmpoptmu[names(tmptmpoptmu)=="value"][pos]
                optmu <- tmpoptmu[[which.max(v)]]
              }
            else
              {
                ## or use importance sampling if not converged
                optmu <- fn.wIS(drawsh, w, controlIS)
              }
          }
        else
          { ## use importance sampling
            optmu <- fn.wIS(drawsh, w, controlIS)
          }

        ## loop over scaling factors
        tmpCV <- tmpw <- tmpph <- tmptheta <- tmpTheta <- tmplnK <- tmplnD <- pos <- NULL
        muh <- optmu$mu

        tmpCV <- Inf
        for (m in 1:nrow(optmu$Sigma))
          { ## loop over scaling matrix Sigma
            Sigmah <- optmu$Sigma[m,]

            ## draw from the new component
            tmptheta <- fn.rmvt(Ns, muh, Sigmah, df)
            tmpTheta <- cbind(Theta, tmptheta)
            tmplnK <- cbind(lnK, KERNEL(tmptheta, log=TRUE, ...))
            
            ## form matrix used in optimization of probabilities
            tmplnD <- matrix(NA, Np, h^2)
            for (i in 1:(h-1))
              {
                tmplnD[,(1+(i-1)*h):((h-1)+(i-1)*h)] <- lnD[,(1+(i-1)*(h-1)):((h-1)+(i-1)*(h-1))]
              }
            for (i in 1:h)
              {
                pos <- seq(from=1+(i-1)*k, length.out=k) 
                tmplnD[,i*h] <- fn.dmvt(tmpTheta[1:Np,pos], muh, Sigmah, df, log=TRUE)
              }
            for (i in 1:(h-1))
              {
                tmplnD[,(h-1)*h+i] <- fn.dmvt(tmptheta[1:Np,], mit$mu[i,], mit$Sigma[i,], df, log=TRUE)
              }
            ## optimization of the probabilities
            optp <- fn.optp(mit$p, tmplnK[1:Np,], tmplnD, control=controloptp)
            tmpph <- optp$p
            
            ## draw from the new mixture
            comp <- sample(1:h, Ns, prob=tmpph, replace=TRUE)
            tmpdrawsh <- tmptmplnK <- tmptmplnD <- tmpw <- pos <- NULL
            for (i in 1:h)
              {
                nh <- length(comp[comp==i])
                if (nh>0)
                  {
                    pos <- seq(from=1+(i-1)*k, length.out=k)
                    tmpdrawsh <- rbind(tmpdrawsh, matrix(tmpTheta[1:nh,pos], ncol=k))
                    tmptmplnK <- c(tmptmplnK, tmplnK[1:nh,i])      
                  }
              }
            tmptmplnD <- dMit(tmpdrawsh, mit=list(p=tmpph, mu=rbind(mit$mu, muh), ## log of mixture for these draws
                                           Sigma=rbind(mit$Sigma, Sigmah),
                                           df=mit$df), log=TRUE)
            tmpw <- fn.computeexpw(tmptmplnK, tmptmplnD)
            tmptmpCV <- fn.CV(tmpw)
            
            if (tmptmpCV<tmpCV) 
              { ## if coefficient of variation better than before
                newCV <- tmpCV <- tmptmpCV
                neww <- tmpw
                newph <- tmpph
                newSigmah <- Sigmah
                newdrawsh <- tmpdrawsh
                newTheta <- tmpTheta
                newlnK <- tmplnK
                newlnD <- tmplnD
              }          
            
            summaryh <- data.frame(h, optmu$method[m], optmu$time, optp$method, optp$time, tmptmpCV)
            names(summaryh) <- nsummary
            if (trace)
              print(summaryh)
            summary <- rbind(summary, summaryh)
          }

        ## add the new component to the mixture
        CV[h] <- newCV
        w <- neww
        ph <- newph
        Sigmah <- newSigmah
        drawsh <- newdrawsh
        Theta <- newTheta
        lnK <- newlnK
        lnD <- newlnD

        mit$p <- ph
        mit$mu <- rbind(mit$mu, muh)
        mit$Sigma <- rbind(mit$Sigma, Sigmah)

        ## stopping criterion
        hstop <- (abs((CV[h]-CV[h-1])/CV[h-1]) <= CVtol)
      }

    ## form the output (add labels)
    names(mit$p) <- paste("cmp", 1:h, sep="")
    dimnames(mit$mu) <- list(paste("cmp", 1:h, sep=""), paste("k", 1:k, sep=""))
    lab <- NULL
    for (i in 1:k)
      for (j in 1:k)
        lab <- c(lab, paste("k", i, "k", j, sep=""))
    dimnames(mit$Sigma) <- list(paste("cmp", 1:h, sep=""), lab)
    
    ## output
    list(CV=as.vector(CV),
         mit=as.list(mit),
         summary=as.data.frame(summary))
  }
