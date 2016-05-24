## Code for solving sequential decision problems for bdpopt

## The code in this file implements a version of the gridding method for solving
## sequential decision problems suggested by Anthony E. Brockwell and Joseph B. Kadane in
## "A Gridding Method for Bayesian Sequential Decision Problems" (2003).

## Construct an object representing a decision problem.

## n.stages = number of stages, an integer greater than or equal to 1.

## post.sample(stage, s, n.sims) should return a list of length n.sims consisting of
## random sample of the parameters from the distribution defined by the
## state value s when being at the given stage.

## pred.sample(stage, thetas, d) should return a list consisting of random samples from the predictive
## distribution of the observations given the values in a list of parameter values thetas and
## given that the decision d is taken. The first element of the list returned should correspond
## to the first element in thetas, etc., and the length of the list returned must equal length(thetas).

## update.state(stage, s, d, Xs) takes a state value s, a decision d and
## a list of observations Xs, and should return a list of updated state
## values obtained when combining the elements of Xs in order with s, given d.

## term.decisions must be a list of length equal to dp$n.stages, the i:th element of which is a list
## that specifies the available terminal decisions at stage i.

## term.obs.decisions must be a list of length equal to dp$n.stages, the i:th element of which is a list
## that specifies the available terminal decisions (with terminal observations) at stage i.

## cont.decisions must be a list of length equal to dp$n.stages, the i:th element of which is a list
## that specifies the available continuation decisions at stage i.

## For any stage i, at least one of the elements of the decision lists must be nonempty, i.e.,
## length(term.decisions[[i]]) + length(term.obs.decisions[[i]]) + length(cont.decisions[[i]]) >= 1, i = 1, ..., dp$n.stages.

## Note that all decisions of the last stage must be terminal decisions, i.e.,
## length(cont.decisions[[dp$n.stages]]) = 0 and
## length(term.decisions[[dp$n.stages]]) + length(term.obs.decisions[[dp$n.stages]]) >= 0.

## term.utility.fun, term.obs.utility.fun and cont.utility.fun are lists of length dp$n.stages which
## givens the stage-wise utility functions.
## The elements of term.utility.fun are functions mapping (d, theta) to a numeric value.
## The elements of term.obs.utility.fun are functions mapping (d, X, theta) to a numeric value.
## The elements of cont.utility.fun are functions mapping (d, X) to a numeric value.

sequential.dp <-
    function(n.stages,
             post.sample,
             pred.sample,
             update.state,
             term.decisions,
             term.obs.decisions,
             cont.decisions,             
             term.utility.fun,
             term.obs.utility.fun,
             cont.utility.fun) {
    list(n.stages = n.stages,
         post.sample = post.sample,
         pred.sample = pred.sample,
         update.state = update.state,
         term.decisions = term.decisions,
         term.obs.decisions = term.obs.decisions,
         cont.decisions = cont.decisions,
         term.utility.fun = term.utility.fun,
         term.obs.utility.fun = term.obs.utility.fun,
         cont.utility.fun = cont.utility.fun)
}

## Check the consistency of a given decision problem.

check.dp.structure <- function(dp) {
    
    if ( !(is.numeric(dp$n.stages) &&
           length(dp$n.stages) == 1 &&
           dp$n.stages >= 1) )
        stop("'n.stages' must be an integer greater than or equal to 1")
    
    if (!(is.list(dp$term.decisions) && length(dp$term.decisions) == dp$n.stages &&
          is.list(dp$term.obs.decisions) && length(dp$term.obs.decisions) == dp$n.stages &&
          is.list(dp$cont.decisions) && length(dp$cont.decisions) == dp$n.stages))
        stop("each of 'term.decisions', 'term.obs.decisions' and 'cont.decisions' must be a list of length 'n.stages'")

    for (i in 1:dp$n.stages)
        if (length(dp$term.decisions[[i]]) + length(dp$term.obs.decisions[[i]]) + length(dp$cont.decisions[[i]]) < 1)
            stop("there must be at least one decision available in each stage")
        
    if (length(dp$cont.decisions[[dp$n.stage]]) > 0)
        stop("the last stage may not contain any continuation decisions")

    if (!(is.list(dp$term.utility.fun) && length(dp$term.utility.fun) == dp$n.stages &&
          is.list(dp$term.obs.utility.fun) && length(dp$term.obs.utility.fun) == dp$n.stages &&
          is.list(dp$cont.utility.fun) && length(dp$cont.utility.fun) == dp$n.stages))
        stop("each of 'term.utility.fun', 'term.obs.utility.fun' and 'cont.utility.fun' must be a list of length 'n.stages'")
}

## Solve a sequential decision problem dp with a grid for the state
## determined by the minimum values, the maximum values and the steps between these.

## n.sims is a vector of sample sizes for the simulation in each stage, its length must equal dp$n.stages.

## state.start, if not equal to NA, will fix the value of the state for the first stage.
## Hence, specifying state.start is like specifying a prior for the parameters.
## Otherwise, stage 1 is handled as all other stages.

optimise.sequential.eu <- function(dp, mins, maxs, steps, n.sims, state.start = NA) {
    check.dp.structure(dp)

    if ( !(is.numeric(mins) && is.atomic(mins) &&
           is.numeric(maxs) && is.atomic(maxs) &&
           is.numeric(steps) && is.atomic(steps) &&
           length(mins) == length(maxs) && 
           length(mins) == length(steps)) )
        stop("each of 'mins', 'maxs' and 'steps' must be numeric, atomic vectors of equal length")

    if (!(is.numeric(n.sims) && is.atomic(n.sims) && length(n.sims) == dp$n.stages))
         stop("'n.sims' must be a numeric, atomic vector of length equal to the number of stages")
    
    if (!identical(state.start, NA) &&
        !(is.numeric(state.start) &&
          is.atomic(state.start) &&
          length(state.start) == length(mins)))
        stop("'state.start' must be a numeric, atomic vector of length equal to 'mins'")
            
    ## Construct a matrix of grid points (each row a point)    
    S.component.grids <- Map(function(x, y, z) seq(x, y, z), mins, maxs, steps)
    S.lengths <- vapply(S.component.grids, length, 0)
    S.grid <- as.matrix(expand.grid(S.component.grids))
    colnames(S.grid) <- NULL
    
    ## Convert a value of the state into the closest grid point index
    to.indices <- function(ss)
        pmin.int(pmax.int(round((ss - mins) / steps) + 1, 1), S.lengths)      

    opt.d <- vector(mode = "list", dp$n.stages)
    opt.eu <- vector(mode = "list", dp$n.stages)
    
    ## Compute the expected utility at a given stage for a terminal decision (without terminal observation)
    term.eu <- function(d, stage, thetas) {       
        mean(vapply(thetas, function(theta) dp$term.utility.fun[[stage]](d, theta), 0))
    }                            

    ## Compute the expected utility at a given stage for a terminal decision (with terminal observation)
    term.obs.eu <- function(d, stage, thetas) {       
        mean(as.numeric(Map(function(X, theta) dp$term.obs.utility.fun[[stage]](d, X, theta),
                            dp$pred.sample(stage, thetas, d), thetas)))
    } 
    
    ## Compute the expected utility at a given stage for a continuation decision
    cont.eu <- function(s, d, stage, thetas) {
        Xs <- dp$pred.sample(stage, thetas, d)
        
        ## Convert the updated state values to grid indices.
        ## These are stored in a matrix, with each column corresponding to one set of grid indices
        s.indices <- t(sapply(dp$update.state(stage, s, d, Xs), to.indices))
        
        ## Compute eu for the current stage plus the eu of continuing optimally
        current.stage.eu <- mean(vapply(Xs, function(X) dp$cont.utility.fun[[stage]](d, X), 0))
        continue.eu <- mean(opt.eu[[stage + 1]][s.indices])  
        current.stage.eu + continue.eu
    }
  
    ## Perform the backward induction
    for (i in (dp$n.stages):1) {
        ## Special handling of first stage if state.start has been provided
        if (i == 1 && !identical(state.start, NA)) {
            opt.d[[i]] <- vector(mode = "numeric", length = 1)
            opt.eu[[i]] <- vector(mode = "numeric", length = 1)

            s <- state.start
            
            ## Sample from the posterior for the parameters theta given the state
            thetas <- dp$post.sample(i, s, n.sims[i])        
        
            eu <- c(vapply(dp$term.decisions[[i]], function(d) term.eu(d, i, thetas), 0),
                    vapply(dp$term.obs.decisions[[i]], function(d) term.obs.eu(d, i, thetas), 0),
                    vapply(dp$cont.decisions[[i]], function(d) cont.eu(s, d, i, thetas), 0))
        
            opt.d[[i]] <- which.max(eu)
            opt.eu[[i]] <- max(eu)     

        }  else {        
            opt.d[[i]] <- array(vector(mode = "numeric"), dim = S.lengths)
            opt.eu[[i]] <- array(vector(mode = "numeric"), dim = S.lengths)
            
            for (j in 1:nrow(S.grid)) {
                s <- S.grid[j,]            
                
                ## Sample from the posterior for the parameters theta given the state
                thetas <- dp$post.sample(i, s, n.sims[i])        
                
                eu <- c(vapply(dp$term.decisions[[i]], function(d) term.eu(d, i, thetas), 0),
                        vapply(dp$term.obs.decisions[[i]], function(d) term.obs.eu(d, i, thetas), 0),
                        vapply(dp$cont.decisions[[i]], function(d) cont.eu(s, d, i, thetas), 0))
        
                opt.d[[i]][j] <- which.max(eu)
                opt.eu[[i]][j] <- max(eu)  
            }
        }
    }    
    
    ## Return (an approximation of) the optimal policy.
    ## This is defined to be a pair of functions, opt.decision and opt.utility,
    ## which take a stage and a state value as input and returns the optimal decision
    ## and optimal utility corresponding to the closest state grid point.
    
    term.max.index <- sapply(dp$term.decisions, length)
    term.obs.max.index <- sapply(dp$term.obs.decisions, length)

    convert.decision.index <- function(stage, di) {    
        if (di <= term.max.index[stage])
            return(dp$term.decisions[[stage]][[di]])

        di <- di - term.max.index[stage]        
        if (di <= term.obs.max.index[stage])
            return(dp$term.obs.decisions[[stage]][[di]])

        di <- di - term.obs.max.index[stage]
        return(dp$cont.decisions[[stage]][[di]])
    }
    
    opt.decision <- function(stage, s) {
        if (stage > dp$n.stages)
            stop("'stage' must be less than or equal to ", dp$n.stages)

        if (identical(state.start, NA) && stage < 1)
            stop("'stage' must be greater than or equal to ", 1)
        
        if (!identical(state.start, NA) && stage <= 1)
            stop("'stage' must be greater than or equal to ", 2)            
            
        convert.decision.index(stage, opt.d[[stage]][t(to.indices(s))])     
    }
    
    opt.utility <- function(stage, s) opt.eu[[stage]][t(to.indices(s))]        

    if (!identical(state.start, NA)) {
         list(opt.stage1.decision = convert.decision.index(1, opt.d[[1]]),
              opt.stage1.utility = opt.eu[[1]],
              opt.decision = opt.decision,
              opt.utility = opt.utility)    
    } else {
        list(opt.decision = opt.decision,
             opt.utility = opt.utility)   
    }
}
