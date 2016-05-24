# plot.ms.irf()
# Plots a combined, colored MSBVAR set of IRFs for comparison.  This
# is in contrast to the one in plot.mc.irf which displays the
# regime-specific IRFs as separate plots.
#
# Patrick T. Brandt
#
# 20140611 : Initial version
# 20140619 : Revised name to omit S3 clashes with plot.mc.irf

"plot.ms.irf" <- function(x, method="Sims-Zha2",
                              component=1, probs=c(0.16, 0.84),
                              varnames=attr(x,"eqnames"),...)
{
    # Get constants
    dd <- dim(x$shortrun)
    N2 <- dd[1]; nsteps <- dd[2]; h <- dd[4]
    m <- sqrt(dd[3])

# output for the irfs and error bands
    out <- vector(mode="list", length=(h+1))
# storage for the ranges
    rgs <- array(0, c(2, m^2, h))

    # Compute the IRF plot quantities and the bounds for the plots.
    for(i in 1:h)
    {
        out[[i]] <- compute.plot.mc.irf(x$shortrun[,,,i], method,
                                    component, probs)
    # Compute the plot ranges we will need
        rgs[,,i] <- apply(out[[i]]$responses, 3, range)

    }

# Find the range for each IRF equation
    minmax <- matrix(0, nrow=m, ncol=2)

# Setup which IRFs go with which eqn for the ranges
    tmp <- (c(1:m^2)%%m)
    tmp[tmp==0] <- m
    indices <- sort(tmp, index.return=T)$ix
    dim(indices) <- c(m, m)

    for(i in 1:m){minmax[i,] <- range(rgs[,indices[,i], 1:h])}

# Set up the color vector for each regime
    col <- rep(1:h, each=3)

# Now loop over each m columns to find the minmax for each column
# responses in the MAR plot.
    j <- 1

# Now do the plotting
    par(mfcol=c(m,m), mai=c(0.25,0.25,0.15,0.25),
        omi=c(0.15,0.75,1,0.15))

    for(i in 1:m^2)
    {
        lims <- ifelse((i-m)%%m==0, m, (i-m)%%m)

    # Combine together the items to be plotted
        for(k in 1:h)
        {
            if(k==1)
            {
                tmp <- out[[k]]$responses[,,i]
            } else {
                tmp <- cbind(tmp, out[[k]]$responses[,,i])
            }
        }
        # plot
        ts.plot(tmp,
                gpars=list(xlab="",ylab="",ylim=minmax[lims,], col=col, ...))

        abline(h=0)

        # shock / response labeling
        if(i<=m){ mtext(varnames[i], side=2, line=3)}
        if((i-1)%%m==0){
            mtext(varnames[j], side=3, line=2)
            j <- j+1
        }
    }

    mtext("Response in", side=2, line=3, outer=T)
    mtext("Shock to", side=3, line=3, outer=T)

    return(invisible(out))
}
