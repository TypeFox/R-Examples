####################################################################################
## PROGRAM: pred.int.R
## PURPOSE: R function to generate predicted intervals based on observed and 
##   simulated data.
## 
## INPUT:
## y (REQUIRED): a numeric vector of outcomes (with at least 2 elements and no 
##   missing values) 
## group: an optional vector of groups. If it exists, it must be the same length as 
##   y with no missing values. If this doesnt exist, the program puts all the y's in 
##   a single group
## N (REQUIRED): a required vector with length equal to the number of distinct groups.
##   The ith element is what we want the number of observations to be for the ith 
##   group after simulation.  (So if y is length 150 and only one group, then to 
##   simulate 50 outcomes we pass N=200.)
## true.y: Either {"observed","no.diff", or a vector of constants}.  Define 
##   mean/proportion used when simulating the data. 
## var.equal: TRUE/FALSE whether assume variance equal in t-test. Default is FALSE
## ref: An optional group name that will serve as the reference group.  Default is 
##   the first alphabetically.   
## data.type (REQUIRED): Either {"t.test","binary"}
## conf.level: Confidence level (between 0 and 1).  Default is 0.95.  This is the
##   confidence level used for the predicted intervals and will also define the 
##   confidence level used for the observed interval unless obs.conf.level is 
##   also used.
## obs.conf.level: Confidence level for the observed intervals.  Default is the
##   same level specified for the predicted intervals in conf.level parameter.
## iters: Number of predicted intervals to generate.  Default is 100
## 
## OUTPUT:
##
## The function returns a list of class "pred.int" having these elements:
##
##  obs.mean, obs.n: Observed mean,N for each group (vector with length = n(groups))
##  sim.n: N simulated for each group (vector with length = n(groups))
##  data.type: Passed through from input parameters
##  conf.level: Confidence level
##  ci: A list of vectors of length 3 that contain the point estimate, lower, and 
##   upper confidence intervals for the observed effect.  There are n(groups)-1 
##   elements in the list (one for each comparison/graph.
##  pi: A list of matrices with 3 columns and iters rows.  The columns are the point 
##   estimate and lower/upper confidence limit for each predicted interval.  There 
##   are (n(groups)-1) matrices in the list (one for each comparison/graph).
##
## MACROS USED:   prop.test.mult.R, sim.basic.R, t.test.mult.R
## CALLED BY:     None
## AUTHOR:        Daniel Muenz (2010) / Ray Griner (2011 - some clean-up)
## CREATION DATE: 2010
## NOTES:
## MODIFICATIONS: 
# [RG20111202] Ray Griner added line to test version and exit if < 2.12.1
# [RG20120107] Ray Griner standardized stop messages per guidelines in Writing R
#  Extensions (R-exts.pdf) since we are now going to submit to CRAN. Mostly this was
#  putting parameter names in single quotes instead of double, but there were also
#  a few cases where we were pasting the message together, and this is discouraged.
# [RG20120201] Ray Griner disabled printing of run-time optional (so that repeated 
#  runs could be identical).  Useful for submitting to CRAN because now test cases 
#  will always have same output when repeated.
# [RG20120830] Ray Griner added support for obs.conf.level and stdized header
####################################################################################

pred.int <- function(y, group=NA, N, true.y="observed", ref=NA,
                     data.type=c("t.test","binary"), var.equal=FALSE,
                     conf.level=0.95, obs.conf.level=NA, iters=100)
{
    ## Time of invokation
    start.time <- as.POSIXlt(Sys.time(), "EST")

    ## If group is NA, we assume that all values in y belong to the
    ## same group.
    if (no.group <- (is.null(group) || (length(group)==1 && is.na(group))))
        group <- rep('A', times=length(y))
    group.names <- as.character( sort(unique(group)) )
    
    ## Find which endpoint type the argument data.type specifies
    ## (allowing for partial matching).
    valid.types <- c("t.test","binary","time.to.event")

    # Actually, I should figure out exactly which version is required 1.8.1 isnt
    # good enough, but 2.12.1 may not be necessary.  These are the two versions we
    # have installed [RG20111202]
    if ((version$major==2 && version$minor<12.1) || version$major<2) stop(paste('Need R version 2.12.1 at least. This is ',R.version.string, '. Try running on daisy.',sep=''))

    data.type <- match.arg(data.type, several.ok=FALSE)
 
    ## Confidence level for observed data should default to level used for predicted data
    if (is.na(obs.conf.level)) 
      obs.conf.level<-conf.level

    ## Argument validation
    ## [RG20120107]S
    if (!is.numeric(y))
        stop("'y' must contain numeric data") 
    else if (any(is.na(y)))
        stop("'y' must not have missing values") 
    else if (length(y) < 2)
        stop("'y' must contain at least 2 elements")
    else if (length(y) != length(group))
        stop("'group', if specified, must have the same length as 'y'")
    else if (any(is.na(group)))
        stop("'group', if specified, must not have missing values")
    else if (!is.numeric(N))
        stop("'N' must contain numeric data")
    else if (length(group.names) != length(N))
        stop("length of 'N' must equal number of distinct values in 'group'")
    else if (is.numeric(true.y) && length(true.y) != length(N))
        stop("if numeric, length of 'true.y' must equal number of distinct values in 'group'")
    else if (data.type=="binary" && length(true.y)>1 && (any(true.y>1) || any(true.y<0)))
        stop("for binary outcomes, if 'true.y' is a number it must be between 0 and 1 (inclusive) for each group")
    else if (class(true.y) == "character" &&
             (length(true.y) != 1 || !(true.y %in% c("observed","no.diff"))))
        stop("invalid value passed to 'true.y'")
    else if (!any(is.na(ref)) && (length(ref) > 1 || !(ref %in% group.names)))
        stop("'ref', if specified, must be a single value from 'group'")
    else if ((length(conf.level) != 1) || is.na(conf.level) || (conf.level <= 0)
             || (conf.level >= 1))
        stop("'conf.level' must be a single number between 0 and 1")
    else if ((length(obs.conf.level) != 1) || is.na(obs.conf.level) || (obs.conf.level <= 0)
             || (obs.conf.level >= 1))
        stop("'obs.conf.level' must be a single number between 0 and 1")
    else if ( !(abs(iters-round(iters))<.Machine$double.eps^0.5) || (iters <= 0))   # Got the idea for this from (is.integer) 
        stop("'iters' (number of simulations) must be a positive integer")
    ##[RG20120107]E

    ## Convert data into a list where each element is a vector of
    ## values for a group.
    data <- list()
    for(name in group.names)
        data[[name]] <- y[group==name]

    ## Get the observed mean, standard dev, and sample size for each
    ## group.
    obs.mean <- sapply(data, mean)
    obs.n <- sapply(data, length)

    ## Number of groups
    n.group <- length(group.names)

    ## Number of data points to generate for each group
    sim.n <- N - obs.n
    names(sim.n) <- group.names
    names(N) <- group.names

    ## Make sure that the projected total sample size is >= the
    ## current observed sample size for each group.
    if ( !is.na(pos <- match(TRUE, sim.n < 0)) ) {
        ## [RG20120107]S
        if (n.group == 1) {
            msg<-sprintf("'N' should be greater than the number of observations in 'y' (currently %d)",obs.n[[pos]]) 
            stop(msg)
        }
        else {
            msg<-sprintf("invalid value in position %d in 'N': should be at least the number of observations in group %s (currently %d)", 
                         pos, group.names[pos], obs.n[[pos]])
            stop(msg)
        }
        ## [RG20120107]E
    }
    ## Also make sure that the projected total sample size is > the
    ## current observed sample size for at least 1 group.
    else if (all(sim.n == 0))
        stop("invalid value in 'N': there is nothing to simulate!") ## [RG20120107] 

    if ((n.group==1) && (true.y=="no.diff")) stop("'true.y' cannot be no.diff when there is only one group") ## [RG20120107]

    ## If true.y is a character string, replace it with numbers.
    if ( true.y == "observed" )
        true.y <- obs.mean
    else if ( true.y == "no.diff" ) {
        true.y <- rep(mean(y), times=n.group)   ## mean(obs.mean) ??? ##
    }
    names(true.y) <- group.names

    ## If no reference group is specified, use the group that comes
    ## first alphanumerically.
    if (is.na(ref)) ref <- group.names[1]

    ## Create a variable, comp.group.names, which is simply the group name
    ## if there's one group, but holds all the group names save the
    ## reference if there are multiple groups.  For example, if there
    ## are 3 groups, A, B, and C, and A is chosen as the reference,
    ## then comp.group.names = c("B", "C").
    if (n.group == 1)
        comp.group.names <- group.names
    else {
        ref.index <- match(ref, group.names)
        comp.group.names <- group.names[-ref.index]
    }

    ## Initialize a list of matrices to hold the predicted intervals.
    ## Each matrix will have 3 columns: the 1st to hold the point
    ## estimates from the simulations, and the 2nd and 3rd to hold the
    ## lower and upper bounds of the predicted intervals,
    ## respectively.  The number of rows in each matrix will be the
    ## number of simulations.
    ci <- list()
    pi <- list()
    for(name in comp.group.names) {
        ci[[name]] <- vector(mode="numeric", length=3)
        names(ci[[name]]) <- c("point","lower","upper")
        pi[[name]] <- matrix(nrow=iters, ncol=3,
                             dimnames=list(1:iters,
                             c("point","lower","upper")))
    }


    ## Now we go separate ways depending on the type of input data.
    if (data.type == "binary") {

        ## Count the number of "successes" in each group.
        obs.count <- sapply(data, sum)

        ## Time of simulation start
        sim.start.time <- as.POSIXlt(Sys.time(), "EST")

        ## Simulate future data, and put it in a matrix which holds
        ## the sum of the success counts from the observed data and
        ## our simulations.  The matrix will have as many columns as
        ## there are groups, and the number of rows will be the number
        ## of simulations.
        pred.count <- do.call(sim.basic,
                              list(data.type=data.type,
                                   group.names=group.names,
                                   obs.n=obs.n,
                                   sim.n=sim.n,
                                   true.y=true.y,
                                   obs.count=obs.count,
                                   iters=iters))

        ## Time of simulation stop
        sim.stop.time <- as.POSIXlt(Sys.time(), "EST")

        ## Calculate confidence intervals based on the observed data
        ## and then predicted intervals based on the observed and
        ## simulated data.
        for(name in comp.group.names) {
            ci[[name]] <- prop.test.mult(x1=obs.count[[name]],
                                         n1=obs.n[[name]],
                                         x2=obs.count[[ref]],
                                         n2=obs.n[[ref]],
                                         samples=n.group,
                                         conf.level=obs.conf.level,
                                         correct=FALSE,
                                         as.vector=TRUE)

            pi[[name]] <- prop.test.mult(x1=pred.count[,name],
                                         n1=N[name],
                                         x2=pred.count[,ref],
                                         n2=N[ref],
                                         samples=n.group,
                                         conf.level=conf.level,
                                         correct=FALSE)
        }

    }
    else if (data.type == "t.test") {

        ## Get the standard deviation of the outcome for each group.
        obs.sd <- sapply(data, sd)

        ## Time of simulation start
        sim.start.time <- as.POSIXlt(Sys.time(), "EST")

        ## Simulate future data, and put it into a list of 2 matrices 
        ## which hold the means and variances obtained from combining 
        ## the observed data with each set of simulated data.  Each matrix 
        ## will have as many columns as there are groups, and the number of
        ## rows will be the number of simulations.
        pred.res <- do.call(sim.basic,
                            list(data.type=data.type,
                                 group.names=group.names,
                                 obs.n=obs.n,
                                 sim.n=sim.n,
                                 N=N,
                                 true.y=true.y,
                                 obs.mean=obs.mean,
                                 obs.sd=obs.sd,
                                 iters=iters))
        pred.mean <- pred.res[[1]]
        pred.var <- pred.res[[2]]

        ## Time of simulation stop
        sim.stop.time <- as.POSIXlt(Sys.time(), "EST")

        ## Calculate confidnce intervals based on the observed data
        ## and then predicted intervals based on the observed and
        ## simulated data.
        for(name in comp.group.names) {
  ####################################################################################
  ## obs.mean, obs.sd, obs.n are all vectors.  So here we are only calculating one
  ##  t-test per comparison (which makes sense, because the output graph will have
  ##  only one confidence interval for the ACTUAL data
  ####################################################################################
            ci[[name]] <- t.test.mult(mean1=obs.mean[name],
                                      var1=obs.sd[name]^2,
                                      n1=obs.n[name],
                                      mean2=obs.mean[ref],
                                      var2=obs.sd[ref]^2,
                                      n2=obs.n[ref],
                                      samples=n.group,
                                      conf.level=obs.conf.level,
                                      var.equal=var.equal,
                                      as.vector=TRUE)

  ####################################################################################
  ## pred.mean1, pred.mean2, and pred.var are matrices.  So this means for the mean1,
  ##  mean2, var1, var2 parameter, we are passing vectors (one element per pred int)
  ## Note that N is a vector with elements equal to the number of groups.  So in n1, 
  ## n2 we are passing a single value that is the N for that group.
  ####################################################################################
            pi[[name]] <- t.test.mult(mean1=pred.mean[,name],
                                      var1=pred.var[,name],
                                      n1=N[name],
                                      mean2=pred.mean[,ref],
                                      var2=pred.var[,ref],
                                      n2=N[ref],
                                      samples=n.group,
                                      var.equal=var.equal,
                                      conf.level=conf.level)
        }

    }
    else if (data.type == "time.to.event") {
    }


    ## Sort all the predicted interval matrices by the values of the
    ## estimates.
    for(name in comp.group.names) {
        order <- order( pi[[name]][,"point"] )
        pi[[name]] <- pi[[name]][order,]
        dimnames(pi[[name]])[[1]] <- 1:iters
    }

    if (no.group) {
        ## Remove the group name if we made it up.
        names(obs.mean) <- NULL
        names(obs.n) <- NULL
        names(ci) <- NULL
        names(sim.n) <- NULL
        names(pi) <- NULL
    }
    else if (n.group != 1) {
        ## Construct names for all the comparisons.
        comp.names <- vector()
        for(name in comp.group.names)
            comp.names[name] <- paste(name, "vs", ref)

        names(pi) <- comp.names
        names(ci) <- comp.names
    }

    ## Put all the important information we've collected into a list.
    result <- list(obs.mean=obs.mean,
                   obs.n=obs.n, ci=ci,
                   sim.n=sim.n, pi=pi,
                   data.type=data.type,
                   conf.level=conf.level,
                   obs.conf.level=obs.conf.level)

    ## Make result of class pred.int.  This lets us use the plot and
    ## print methods for this class.
    class(result) <- "pred.int"

    ## Time of completion
    stop.time <- as.POSIXlt(Sys.time(), "EST")

    whole.time <- stop.time - start.time
    sim.time <- sim.stop.time - sim.start.time

    ## [RG20120107]S 
    msg<-sprintf("Whole procedure took %f secs. Simulation took %f secs.", round(whole.time, digits=3), round(sim.time, digits=3))
    ## message(msg)  ## [RG20120201]
    ## [RG20120107]E 

    ## The end.
    return(result)
}
