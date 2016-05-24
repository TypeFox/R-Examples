#####################################################################
## PROGRAM: sim.basic.R
## PURPOSE: Simulate future data for PIPs and combine with current data
## 
## INPUT:
##  data.type: "binary" or "t.test"
##  group.names: Names of groups to simulate
##  iters: Number of iterations (simulations) to run
##  sim.n: N to simulate in each group
##  true.y:  True parameter to simulate for each group. For binary,
##   this is the true proportion.  For t.test, this is the true
##   mean
##  obs.count: (for data.type="binary") Counts successes in observed data 
##  N: (for data.type="t.test") Overall N for each group 
##   (redundant, just sim.n + obs.n)
##  obs.n, obs.mean, obs.sd: (for data.type="t.test") the N, mean
##   and std for each observed group
##
## OUTPUT:        List containing the simulated resulst
## MACROS USED:   None
## CALLED BY:     pred.int
## AUTHOR:        Daniel Muenz
## CREATION DATE: 2010
## NOTES:
## MODIFICATIONS:
##  [RG20120107] Ray Griner standardized stop messages
##  [RG20120830] Ray Griner stdize program header
#####################################################################
sim.basic <- function(data.type, group.names,
                      obs.n, sim.n, N,
                      true.y, obs.count,
                      obs.mean, obs.sd, iters, ...)
{
    n.group<-length(group.names)
    if (data.type == "binary") {
        ## Initialize a matrix to hold the sum of the success counts
        ## from the observed data and our simulations.  The matrix
        ## will have as many columns as there are groups, and the
        ## number of rows will be the number of simulations.
        pred.count <- matrix(nrow=iters, ncol=n.group,
                             dimnames=list(1:iters,group.names))

        ## Simulate future data.
        for(name in group.names) {
            if (sim.n[[name]] > 0)
                pred.count[,name] <- obs.count[[name]] +
                                     rbinom(n=iters,
                                            size=sim.n[[name]],
                                            prob=true.y[[name]])
            else
                pred.count[,name] <- obs.count[[name]]
        }

        ## Return results.
        return(pred.count)
    }
    else if (data.type == "t.test") {
        if (any(N!=(sim.n+obs.n))) stop("'N' should equal 'sim.n'+'obs.n'")

        ## Initialize 4 matrices: the 1st two to hold the means and
        ## sums of squared deviations (ssd) for our simulations, and
        ## the next two to hold the means and variances obtained from
        ## combining the observed data with each set of simulated
        ## data.  Each matrix will have as many columns as there are
        ## groups, and the number of rows will be the number of
        ## simulations.
        sim.mean <- matrix(0, nrow=iters, ncol=n.group,
                           dimnames=list(1:iters,group.names))
        sim.ssd  <- matrix(0, nrow=iters, ncol=n.group,
                           dimnames=list(1:iters,group.names))
        pred.mean <- matrix(nrow=iters, ncol=n.group,
                            dimnames=list(1:iters,group.names))
        pred.var  <- matrix(nrow=iters, ncol=n.group,
                            dimnames=list(1:iters,group.names))

        ## Get the subset of groups for which data are to be simulated.
        sim.group.names <- group.names[sim.n > 0]

        ## Simulate new data "iters" times.
        for(i in 1:iters) {
            ## Generate random data for each group, and calculate the mean
            ## and sum of squared deviations (we don't actually need to
            ## calculate the variance).
            for(name in sim.group.names) {
                sim.data <- rnorm(n=sim.n[[name]], mean=true.y[[name]],
                                  sd=obs.sd[[name]])
                sim.mean[i,name] <- mean(sim.data)
                sim.ssd[i,name] <- sum((sim.data - sim.mean[i,name])^2)
            }
        }

        ## Calculate the predicted mean by taking a weighted average
        ## of the observed and simulated means.
        pred.mean <- t((obs.mean*obs.n + t(sim.mean)*sim.n) / N)

        ## Calculate the predicted variance.
        pred.var <- t( ( obs.sd^2*(obs.n-1) +
                         (obs.mean-t(sim.mean))^2*obs.n*sim.n/N +
                         t(sim.ssd) ) / (N-1) )

        ## Return list containing predicted means as first item and
        ## predicted variance as second.
        return( list(pred.mean, pred.var) )
    }
    else if (data.type == "time.to.event") {
    }
}
