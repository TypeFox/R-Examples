# Generate the impulse responses -- requires a posterior sample of the
# A0 to be drawn already.  Could probably put an option / overload in here that
# if there are no A0 supplied that it draws them.

"mc.irf" <- function(varobj, nsteps, draws=1000, A0.posterior=NULL, sign.list=rep(1,ncol(varobj$Y)))
{
    if(inherits(varobj, "VAR")){
        return(mc.irf.VAR(varobj=varobj, nsteps=nsteps, draws=draws))
    }
    if(inherits(varobj, "BVAR")){
        return(mc.irf.BVAR(varobj=varobj, nsteps=nsteps, draws=draws))
    }
    if(inherits(varobj, "BSVAR")){
        return(mc.irf.BSVAR(varobj=varobj, nsteps=nsteps,
                            A0.posterior=A0.posterior, sign.list=sign.list))
    }
    if(inherits(varobj, "MSBVAR")){
        return(mc.irf.MSBVAR(varobj=varobj, nsteps=nsteps,
                             draws=length(varobj$ss.sample)))
    }
}

"mc.irf.VAR" <- function(varobj, nsteps, draws)
{ output <- .Call("mc.irf.var.cpp", varobj, as.integer(nsteps),
                  as.integer(draws))
  attr(output, "class") <- c("mc.irf", "mc.irf.VAR")
  attr(output, "eqnames") <- attr(varobj, "eqnames")
  return(output)
}


"mc.irf.BVAR" <- function(varobj, nsteps, draws)
{
    output <- mc.irf.VAR(varobj, nsteps, draws)
    attr(output, "class") <- c("mc.irf", "mc.irf.BVAR")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}


"mc.irf.BSVAR" <- function(varobj, nsteps, A0.posterior, sign.list)
{ m<-dim(varobj$ar.coefs)[1]  # Capture the number of variablesbcoefs <- varobj$Bhat

  # Check length of sign.list
  if(length(sign.list)!=m) {
      stop("sign.list has wrong number of elements.  Must be m")
  }

  p<-dim(varobj$ar.coefs)[3]    # Capture the number of lags
  ncoef <- dim(varobj$B.posterior)[1]
  n0 <- varobj$n0
  n0cum <- c(0,cumsum(n0))
  N2 <- A0.posterior$N2

  # Get the covar for the coefficients
  XXinv <- chol(solve(varobj$Hpinv.posterior[[1]]))

  # storage for the impulses and the sampled coefficients.
  impulse <- matrix(0,nrow=N2, ncol=(m^2*nsteps))

  output <- .Call("mc.irf.bsvar.cpp", A0.posterior$A0.posterior,
                  as.integer(nsteps), as.integer(N2), as.integer(m),
                  as.integer(p), as.integer(ncoef), as.integer(n0),
                  as.integer(n0cum), XXinv, varobj$Ui, varobj$P.posterior,
                  sign.list)
  attr(output, "class") <- c("mc.irf", "mc.irf.BSVAR")
  attr(output, "eqnames") <- attr(varobj, "eqnames")
  return(output)
}

######################################
# Helper functions for mc.irf.MSBVAR
######################################

# Get an A0 -- function to get the reduced form VAR A(0) to capture
# the contemp effects.

getA0 <- function(z, m, h)
{
    tmp <- apply(matrix(z, ncol=h), 2, xpnd)
    tmp <- array(tmp, c(m,m,h))
    tmp <- array(apply(tmp, 3, chol), c(m,m,h))
    aperm(tmp, c(2,1,3))
}




"irf.var.from.beta" <-
function(A0,bvec,nsteps)
{   m <- ncol(A0)
    p <- (length(bvec))/m^2;
    bmtx <- matrix(bvec,ncol=m)
    ar.coefs<-t(bmtx)     # extract the ar coefficients
    dim(ar.coefs)<-c(m,m,p)              # push ar coefs into M x M x P array
    ar.coefs<-aperm(ar.coefs,c(2,1,3))   # reorder array so columns
                                         # are for eqn
    impulses <- irf.VAR(list(ar.coefs=ar.coefs),nsteps,A0=A0)$mhat
    impulses <- aperm(impulses, c(3,1,2)) # flips around the responses
                                          # to stack them as series
                                          # for each response, so the
                                          # first nstep elements are
                                          # the responses to the first
                                          # shock, the second are the
                                          # responses of the first
                                          # variable to the next
                                          # shock, etc.

    dim(impulses)<-c((m^2)*nsteps,1)
    return(impulses)
  }


"mc.irf.MSBVAR" <- function(varobj, nsteps, draws)
{
    # Get constants
    m <- varobj$m
    p <- varobj$p
    h <- varobj$h
    N2 <- draws

    # Get the A0 or structural shocks
    Sigmavec <- varobj$Sigma.sample

    cat("Running setup tasks to sort regime parameters.\n")

    # Set them up in regime specific arrays
    A0s <- array(sapply(1:N2, function(i) {getA0(Sigmavec[i,], m, h)}),
                 c(m, m, h, N2))

    # Now set up the AR coefficients vector.  Need to pluck out the
    # intercepts from the posterior and assemble the AR coefs in a vector
    # for each regime

    # Intercept indices
    mpplus1 <- m*p + 1
    nc <- m*mpplus1
    ii <- seq(mpplus1, by=mpplus1, length=m)  # Intercept indices

    # Split Beta by regimes

    Beta <- array(varobj$Beta.sample, c(N2, nc, h))

    # Remove the intercepts -- just the AR part
    AR <- Beta[,-ii,]
    rm(Beta)

    cat("Computing long run regime probabilities.\n")
    # Get summaries of the long run ergodic probabiities
    lrQ <- sapply(1:N2,
                  function(i) {
                      steady.Q(matrix(varobj$Q.sample[i,],h,h))})

    # Storage for the impulse responses: draws, steps, response-shock
    # combinations, regimes array

    tmp <- array(0, c(N2, nsteps, m^2, h))
    outputavg <- array(0, c(N2, nsteps, m^2))

    # Loop over the draws and regimes to compute the IRF for each
    # regime and in equilibrium

    # Computing the IRF for each regime
    cat("Beginning to compute the IRFs for each regime, and averaged over the regime probabilities.\n")
    for(i in 1:N2)
    {
        for(j in 1:h)
        {
            tmp[i,,,j] <- irf.var.from.beta(A0s[,,j,i], AR[i,,j], nsteps)

            # Now take long run averages and store them
            outputavg[i,,] <- tmp[i,,,j] + tmp[i,,,j]*lrQ[j,i]
        }
        if(i%%1000==0) cat("Finished ", 100*i/N2, " percent\n", sep="")
    }

    output <- list(shortrun=tmp, longrun=outputavg)
    attr(output, "class") <- c("mc.irf", "mc.irf.MSBVAR")
    attr(output, "eqnames") <- attr(varobj, "eqnames")
    return(output)
}

######################################################################
# Main classes plotting function for impulse responses computed with
# mc.irf and it classes cousins.
######################################################################

"plot.mc.irf" <- function(x, method=c("Sims-Zha2"), component=1,
                          probs=c(0.16,0.84),
                          varnames=attr(x, "eqnames"),
                          regimelabels=NULL, ask=TRUE, ...)
{
    # VAR
    if(inherits(x, "mc.irf.VAR"))
    { tmp <- plot.mc.irf.VAR(x, method=method, component=component,
                             probs=probs, varnames=varnames,...) }
    # BVAR
    if(inherits(x, "mc.irf.BVAR"))
    { tmp <- plot.mc.irf.BVAR(x, method=method, component=component,
                              probs=probs, varnames=varnames,...) }
    # BSVAR
    if(inherits(x, "mc.irf.BSVAR"))
    { tmp <- plot.mc.irf.BSVAR(x, method=method, component=component,
                               probs=probs, varnames=varnames,...) }
    # MSBVAR
    if(inherits(x, "mc.irf.MSBVAR"))
    { tmp <- plot.mc.irf.MSBVAR(x, method=method,
                                component=component, probs=probs,
                                varnames=varnames,
                                regimelabels=regimelabels, ask,...) }
    # Get the eigen-fractions, if relevant and available.
    return(invisible(tmp))
 }

# Use a single function to compute the IRF across each kind of VAR.
# Then the computation of the IRF can be separated from the plotting
# functions for more user control.  This is hidden from the user, but
# called inside of the various functions needed to draw the IRF plots
# with error bands.

"compute.plot.mc.irf" <- function(x, method, component, probs)
{

    mc.impulse <- x
    m <- sqrt(dim(mc.impulse)[3])
    nsteps <- dim(mc.impulse)[2]
    draws <- dim(mc.impulse)[1]

  # Storage
    irf.ci <- array(0,c(nsteps,length(probs)+1,m^2))

  # Compute the IRF confidence intervals for each element of MxM
  # array

  # There are multiple methods for doing this.

  # Monte Carlo / Bootstrap percentile method
    if (method=="Percentile")
    {
        eigen.sum <- 0
        for(i in 1:m^2)
        {
            irf.bands <- t(apply(mc.impulse[,,i], 2, quantile, probs))
            irf.mean <- apply(mc.impulse[,,i], 2, mean)
            irf.ci[,,i] <- cbind(irf.bands, irf.mean)
        }
    }
    if (method=="Normal Approximation")
    {
        eigen.sum <- 0
        for (i in 1:m^2)
        {
            irf.mean <- apply(mc.impulse[,,i], 2, mean)
            irf.var <- apply(mc.impulse[,,i], 2, var)
            irf.bands <- irf.mean + matrix(rep(qnorm(probs), each=nsteps), nrow=nsteps)*irf.var
            irf.ci[,,i] <- cbind(irf.bands, irf.mean)
        }
    }

  # Sims and Zha symmetric eigen decomposition (assumes normality
  # approximation) with no accounting for the correlation across the
  # responses.  Method does account for the correlation over time.

    if (method=="Sims-Zha1")
    {
        eigen.sum <- matrix(0, m^2, nsteps)
        for(i in 1:m^2)
        {
            decomp <- eigen(var(mc.impulse[,,i]), symmetric=T)
            W <- decomp$vectors
            lambda <- decomp$values
            irf.mean <- apply(mc.impulse[,,i], 2, mean)
            irf.bands <- irf.mean + W[,component]*matrix(rep(qnorm(probs), each=nsteps), nrow=nsteps)*sqrt(lambda[component])
            irf.ci[,,i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i,] <- 100*lambda/sum(lambda)
        }
    }

  # Sims and Zha asymmetric eigen decomposition (no normality
  # assumption) with no accounting for the correlation across the
  # responses.  Method does account for correlation over time.
    if (method=="Sims-Zha2")
    {
        eigen.sum <- matrix(0, m^2, nsteps)
        for(i in 1:m^2)
        {
            decomp <- eigen(var(mc.impulse[,,i]), symmetric=T)
            W <- decomp$vectors
            lambda <- decomp$values
            gammak <- mc.impulse[,,i]*(W[component,])
            gammak.quantiles <- t(apply(gammak, 2, quantile, probs=probs))
            irf.mean <- apply(mc.impulse[,,i], 2, mean)
            irf.bands <- irf.mean + gammak.quantiles
            irf.ci[,,i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i,] <- 100*lambda/sum(lambda)
        }
    }

  # Sims and Zha asymmetric eigen decomposition with no normality
  # assumption and an accounting of the temporal and cross response
  # correlations.

    if (method=="Sims-Zha3")
    {
        eigen.sum <- matrix(0, m^2, m^2*nsteps)

      # Stack all responses and compute one eigen decomposition.
        stacked.irf <- array(mc.impulse, c(draws, m^2*nsteps))
        decomp <- eigen(var(stacked.irf), symmetric=T)
        W <- decomp$vectors
        lambda <- decomp$values
        gammak <- stacked.irf*W[component,]
        gammak.quantiles <- apply(gammak, 2, quantile, probs)
        irf.mean <- matrix(apply(stacked.irf, 2, mean),
                           nrow=length(probs),
                           ncol=dim(stacked.irf)[2],
                           byrow=T)
        irf.bands <- irf.mean + gammak.quantiles

      # Reshape these....
        irf.ci <- array((rbind(irf.bands,irf.mean[1,])),
                        c(length(probs)+1, nsteps, m^2))
        irf.ci <- aperm(irf.ci, c(2, 1, 3))
        eigen.sum <- 100*lambda/sum(lambda)
    }
    return(list(responses=irf.ci, eigenvector.fractions=eigen.sum,
                m=m, nsteps=nsteps))
}

"plot.mc.irf.VAR" <- function(x, method=method, component=component,
                              probs=probs, varnames=varnames,...)
{
    tmp <- compute.plot.mc.irf(x, method, component, probs)

    # Get out the main components we need
    m <- tmp$m
    nsteps <- tmp$nsteps
    irf.ci <- tmp$responses
    eigen.sum <- tmp$eigenvector.fractions

    # Compute the bounds for the plots
    minmax <- matrix(0, nrow=m, ncol=2)
    within.plots <- apply(irf.ci, 3, range)

    tmp <- (c(1:m^2)%%m)
    tmp[tmp==0] <- m
    indices <- sort(tmp, index.return=T)$ix
    dim(indices) <- c(m, m)

    for(i in 1:m){minmax[i,] <- range(within.plots[,indices[,i]])}
    # Now loop over each m columns to find the minmax for each column
    # responses in the MAR plot.
    j <- 1
    # Plot the results
    par(mfcol=c(m,m), mai=c(0.25,0.25,0.15,0.25),
        omi=c(0.15,0.75,1,0.15))

    for(i in 1:m^2)
    {
        lims <- ifelse((i-m)%%m==0, m, (i-m)%%m)
        ts.plot(irf.ci[,,i],
                gpars=list(xlab="",ylab="",ylim=minmax[lims,]), ...)
        abline(h=0)

        if(i<=m){ mtext(varnames[i], side=2, line=3)}
        if((i-1)%%m==0){
            mtext(varnames[j], side=3, line=2)
            j <- j+1
        }
    }

    mtext("Response in", side=2, line=3, outer=T)
    mtext("Shock to", side=3, line=3, outer=T)

    # Put response names on the eigenvector fractions
    if(method == "Sims-Zha1" | method == "Sims-Zha2")
    {
        if(is.null(varnames)==T) varnames <- paste("V", seq(1:m), sep = "")
        shock.name <- rep(varnames, m)
        response.name <- rep(varnames, each=m)
        eigen.sum <- cbind(shock.name, response.name, as.data.frame(eigen.sum))
        colnames(eigen.sum) <- c("Shock","Response", paste("Component", seq(1:nsteps)))
    }
    if(method == "Sims-Zha3")
    { names(eigen.sum) <- c(paste("Component", seq(1:m^2*nsteps))) }

    # Return
    return(list(responses=irf.ci, eigenvector.fractions=eigen.sum))
}

"plot.mc.irf.BVAR" <- function(x, method=method, component=component,
                               probs=probs, varnames=varnames, ...)
{
    plot.mc.irf.VAR(x, method, component,
                    probs, varnames, ...)
}


# BSVAR model IRFs

"plot.mc.irf.BSVAR" <- function(x, method=method, component=component,
                                probs=probs, varnames=varnames, ...)
{
    m <- sqrt(dim(x)[3])
    nsteps <- dim(x)[2]
    draws <- dim(x)[1]
    irf.ci <- array(0, c(nsteps, length(probs) + 1, m^2))
    if (method == "Percentile") {
        eigen.sum <- 0
        for (i in 1:m^2) {
            irf.bands <- t(apply(x[, , i], 2, quantile,
                probs))
            irf.mean <- apply(x[, , i], 2, mean)
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
        }
    }
    if (method == "Normal Approximation") {
        eigen.sum <- 0
        for (i in 1:m^2) {
            irf.mean <- apply(x[, , i], 2, mean)
            irf.var <- apply(x[, , i], 2, var)
            irf.bands <- irf.mean + matrix(rep(qnorm(probs),
                each = nsteps), nrow = nsteps) * irf.var
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
        }
    }
    if (method == "Sims-Zha1") {
        eigen.sum <- matrix(0, m^2, nsteps)
        for (i in 1:m^2) {
            decomp <- eigen(var(x[, , i]), symmetric = T)
            W <- decomp$vectors
            lambda <- decomp$values
            irf.mean <- apply(x[, , i], 2, mean)
            irf.bands <- irf.mean + W[, component] * matrix(rep(qnorm(probs),
                each = nsteps), nrow = nsteps) * sqrt(lambda[component])
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i, ] <- 100 * lambda/sum(lambda)
        }
    }
    if (method == "Sims-Zha2") {
        eigen.sum <- matrix(0, m^2, nsteps)
        for (i in 1:m^2) {
            decomp <- eigen(var(x[, , i]), symmetric = T)
            W <- decomp$vectors
            lambda <- decomp$values
            gammak <- x[, , i] * (W[component, ])
            gammak.quantiles <- t(apply(gammak, 2, quantile,
                probs = probs))
            irf.mean <- apply(x[, , i], 2, mean)
            irf.bands <- irf.mean + gammak.quantiles
            irf.ci[, , i] <- cbind(irf.bands, irf.mean)
            eigen.sum[i, ] <- 100 * lambda/sum(lambda)
        }
    }
    if (method == "Sims-Zha3") {
        eigen.sum <- matrix(0, m^2, m^2 * nsteps)
        stacked.irf <- array(x, c(draws, m^2 * nsteps))
        decomp <- eigen(var(stacked.irf), symmetric = T)
        W <- decomp$vectors
        lambda <- decomp$values
        gammak <- stacked.irf * W[component, ]
        gammak.quantiles <- apply(gammak, 2, quantile, probs)
        irf.mean <- matrix(apply(stacked.irf, 2, mean), nrow = length(probs),
            ncol = dim(stacked.irf)[2], byrow = T)
        irf.bands <- irf.mean + gammak.quantiles
        irf.ci <- array((rbind(irf.bands, irf.mean[1, ])), c(length(probs) +
            1, nsteps, m^2))
        irf.ci <- aperm(irf.ci, c(2, 1, 3))
        eigen.sum <- 100 * lambda/sum(lambda)
    }

    minmax <- matrix(0, nrow = m, ncol = 2)
    within.plots <- apply(irf.ci, 3, range)

    tmp <- (c(1:m^2)%%m)
    tmp[tmp==0] <- m
    indices <- sort(tmp, index.return=T)$ix
    dim(indices) <- c(m, m)

    for(i in 1:m)
      { minmax[i,] <- range(within.plots[,indices[,i]])
      }

    # Now loop over each m columns to find the minmax for each column
    # responses in the MAR plot.
    j <- 1
    # Plot the results
    par(mfcol=c(m,m),mai=c(0.15,0.2,0.15,0.15), omi=c(0.15,0.75,1,0.15))
    for(i in 1:m^2)
      {
          lims <- ifelse((i-m)%%m==0, m, (i-m)%%m)
          ts.plot(irf.ci[,,i],
                  gpars=list(xlab="",ylab="",ylim=minmax[lims,]))

          abline(h=0)

          if(i<=m)
            { mtext(varnames[i], side=2, line=3) }
          if((i-1)%%m==0)
            { mtext(varnames[j], side=3, line=2)
              j <- j+1
            }
        }

    # Add row and column labels for graph

    mtext("Response in", side = 2, line = 3, outer = T)
    mtext("Shock to", side = 3, line = 3, outer = T)

    # Return eigendecomposition pieces if necessary

    if (method == "Sims-Zha1" | method == "Sims-Zha2") {
        if (is.null(varnames) == T)
            varnames <- paste("V", seq(1:m), sep = "")
        shock.name <- rep(varnames, m)
        response.name <- rep(varnames, each = m)
        eigen.sum <- cbind(shock.name, response.name, as.data.frame(eigen.sum))
        colnames(eigen.sum) <- c("Shock", "Response", paste("Component",
            seq(1:nsteps)))
    }
    if (method == "Sims-Zha3") {
        names(eigen.sum) <- c(paste("Component", seq(1:m^2 *
            nsteps)))
    }

    return(list(responses = irf.ci, eigenvector.fractions = eigen.sum))
}


# MSBVAR IRF plot
"plot.mc.irf.MSBVAR" <- function(x, method=method,
                                 component=component, probs=probs,
                                 varnames=varnames,
                                 regimelabels=regimelabels, ask=ask, ...)
{
    # Subset out the irfs by the number of regimes and
    # variables. Recall that the dimensions of the irf array in "x"
    # are N2 x nsteps x (shock*response) x no. regimes.

    # Start by plotting each of the shortrun plots
    d <- dim(x$shortrun)
    h <- d[4]

    # Storage / list for the IRF output summaries
    out <- vector(mode="list", length=(h+1))

    # Make regime labels if they have not been provided
    if(is.null(regimelabels)) regimelabels <- paste("Regime", 1:h)

    # Main loop for the plot.
    for (i in 1:h)
    {
        # plot each short run irf, per regime, and a label as
        # such.

       # Store the summary / CI / eigendecomposition if it is done.
        out[[i]] <- plot.mc.irf.VAR(x=x$shortrun[,,,i], method=method,
                                    component=component,
                                    probs=probs, varnames=varnames, ...)
        # Add regime label to plot
        mtext(regimelabels[i], side=3, line=4.5, outer=T)

        devAskNewPage(ask=ask)
    }

    # Now plot the long-run of regime probability averaged IRF, with
    # error bands

    tmp <- plot.mc.irf.VAR(x$longrun, method=method, component=component,
                           probs=probs, varnames=varnames)
    mtext("Regime averaged IRF", side=3, line=4.5, outer=TRUE)

    out[[h+1]] <- tmp

    # Add regime names to output
    names(out) <- c(regimelabels, "Regime averaged")

    return(out)
}
