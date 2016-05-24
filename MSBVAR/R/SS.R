## SS.R:  R-side functions for the SS class objects that describe the
##        state-space for an MSBVAR model
##
## Date:  20071115 : Initial version
##        20081014 : Revised to use the bit package for storage.
##        20100617 : Updated to move sampler to gibbs.msbvar.  This
##                   file only now handles the posterior
##                   post-processing of the SS classed objects.
## Description:
##
## State spaces computed and stored in compressed
## form to reduce the computation and memory overhead involved in
## computing and storing state space draws.  R side storage of these
## state spaces is done using the bit package.  Functions here are
## used to summarize and plot the posterior state spaces.


######################################################################
# Summary functions for MS state-spaces from the BHLK filter or other
# estimators for the state-space
######################################################################

######################################################################
# Sums of the number of draws in each state as a function of the number
# of periods TT
######################################################################

sum.SS <- function(x, ...)
{
    h <- x$h
    N <- length(x$ss.sample)
    ss <- x$ss.sample
    TTh <- virtual(x$ss.sample[[1]])$Length
    TT <- TTh/(h-1)

    # Now find the sums
    sums <- apply(matrix(unlist(lapply(x$ss.sample, as.integer)), nrow=TTh, ncol=N), 1, sum)
    sums <- matrix(sums, TT, h-1)
    sums <- cbind(sums, rep(N, TT) - rowSums(sums))  # Recover the h'th regime
    return(sums)
}


######################################################################
# Mean regime probability for each state
######################################################################
mean.SS <- function(x, ...){
    sums <- sum.SS(x)
    N <- length(x$ss.sample)
    return(sums/N)
}

######################################################################
# Plot the mean posterior regime probabilities
######################################################################

plot.SS <- function(x, ylab="State Probabilities", ...)
{
    tmp <- mean.SS(x)
    shift <- x$p/attr(x, "freq")
    plot(ts(tmp,
            start=attr(x, "start")+shift,
            end=attr(x, "end"), frequency=attr(x, "freq")),
###            deltat=attr(x, "freq")),
         plot.type="single", col=1:ncol(tmp),
         ylim=c(0,1), ylab=ylab, ...)
    abline(h=0.5, lty=2, ...)
}

# Legacy R code version of some state-space generation functions --
# for legacy and exposition purposes.

# Function to generate the vector of binary indicators of the state
# space.  This is a state-space generation function for one
# observation.
#
# p = filtered probabilities of each regime for an observation
# m = number of regimes.
#
# Returns a vector of 0-1 that indicate the regime.  A one is returned
# for the regime that the observation falls into.

## bingen <- function(p, Q, st1)
## {
##     h <- dim(Q)[1]
##     i <- 1

##     while(i<h)
##     { pr0 <- p[i]*Q[st1,i]/sum(p[i:h]*Q[st1,i:h])

##       if(runif(1)<=pr0)
##       { return(diag(h)[i,]) } else { st1 <- i <- i+1 }
##   }
##     return(diag(h)[h,])
##   }

## # Multi-move Gibbs sampler for the state space
## # filtered.prob  = BHLK.filter probabilities from BHLK.filter
## # P = Markov transition matrix

## generate.states <- function(filtered.prob, Q)
##   { TT <- nrow(filtered.prob)
##     h <- ncol(filtered.prob)

##     # storage
##     ss <- matrix(0, nrow=h, ncol=TT)

##     # generate the TT state
##     ss[,TT] <- bingen(filtered.prob[TT,], Q, 1)

##     for (t in (TT-1):1)
##       {
##           ss[,t] <- bingen(filtered.prob[t,], Q, which(ss[,t+1]==1))
##       }

##     return(t(ss))
##   }


## SS.draw.R <- function(xi, b, Q, m, h, n0, init.model, fp)
## {

##     SS <- hregime.SS(h, m, fp, b, xi, n0, init.model, method="update.Q")

##     # Form Xi for the regimes
##     xi.tmp <- matrix(xi^2, m, h)
##     Xik <- array(0, c(m,m,h))

##     for(i in 1:h)
##     {

##         Xik[,,i] <- solve(diag(xi.tmp[,i]))
##     }

##     fp <- BHLK.filter(SS$e, Xik, Q)

##     ss <- generate.states(fp, Q)
##     return(ss)
## }
