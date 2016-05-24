## Simulation (using JAGS) model object for bdpopt

## Constructor for a probabilistic model using simulation.

## sim.points, sim.means and sim.SEs are all set to "NA" when a new model is created.
## When an evaluation has been performed, these become
## sim.points = A matrix of grid point locations (the columns) for the data.
##              The row names gives the corresponding variables.
## sim.means = A vector of sample means.
## sim.SEs = A vector of standard errors corresponding to the sample means.

## regression.fun and gpr.hyper.params are both set to "NA" when a new model is created.
## regression.fun is set equal to a regression approximation function after fit.gpr or fit.loess has been called,
## and gpr.hyper.params is set equal to the hyper.parameters of the gpr model.

## The typical working sequence is thus:

## 1. Create a new simulation model object.
## 2. Perform simulation over a grid.
## 3. Optionally, fit a regression model to the simulated values.
## 4. Plot the results and optimise over the decision variables.   
    
sim.model <- function(model.file, data) {
    ## Check that model.file is a valid file name
    if( !(is.character(model.file) && length(model.file) == 1) )
        stop("'model.file' must be a character vector consisting of a single string")

    if (!file.exists(model.file))
        stop("could not find the file specified by 'model.file'")
    
    ## If data is a list, use it directly. If a file name, read it.   
    if (is.character(data) && length(data) == 1) {
        if (!file.exists(data))
            stop("could not find the file specified by the 'data'")

        data <- rjags::read.jagsdata(data)
        
    } else if (!is.list(data)) {
        stop("'data' must be a character vector consisting of a single string or a named list")
    }         
    
    structure(        
        list(model.file = model.file,
             data = data,
             grid.spec.list = NA,
             sim.points = NA,
             sim.means = NA,             
             sim.SEs = NA,
             regression.fun = NA,
             gpr.hyper.params = NA),
        class = c("sim.model")
        )
}

## Print a summary description of a sim.model object
print.sim.model <- function(x, ...) {
    cat(paste("Model file:", x$model.file), "\n\n")
    
    if (identical(x$sim.points, NA))        
        cat("No simulations performed.\n\n")
    else
        cat(paste("Simulation performed. Number of grid points:", length(x$sim.means), "\n\n"))

    if (identical(x$regression.fun, NA))        
        cat("No regression performed.\n")
    else {
        cat("Regression performed.\n")
        
        if (!identical(x$gpr.hyper.params, NA)) {
            cat("Hyperparameters (standard deviation and lengths):\n")

            cat(paste("sd:", x$gpr.hyper.params[1], "\n"))

            hps <- x$gpr.hyper.params[-1]
            rns <- rownames(x$sim.points)
            for (i in 1:length(hps))
                cat(paste(rns[i], ": ", hps[i], "\n", sep = ""))
        }
    }
}

## Generic function for diagnosing JAGS output for estimating expected utility.
diag.mcmc.list <- function(model, utility.fun, data, n.iter, n.burn.in, n.adapt = 1000, n.chains = 1, inits = NULL) {
    UseMethod("diag.mcmc.list", model)
}

## Method for diagnosing JAGS output for estimating expected utility.
diag.mcmc.list.sim.model <- function(model, utility.fun, data, n.iter, n.burn.in, n.adapt = 1000, n.chains = 1, inits = NULL) {    
    ## Get the total set of data and the names
    all.data <- c(model$data, data)
    all.data.names <- names(all.data)

    ## Extract the argument names of the utility function,
    ## and construct separate lists of the param and data names used
    arg.names <- names(formals(utility.fun))    
    arg.names.params <- setdiff(arg.names, all.data.names)
    arg.names.data <- setdiff(arg.names, arg.names.params)   
    
    ## Load jags model from file, update for "burn in", draw coda samples and return the result
    jm <-
        if (identical(inits, NULL))
            rjags::jags.model(model$model.file, data = all.data,
                              n.chains = n.chains, n.adapt = n.adapt, quiet = TRUE)
        else
            rjags::jags.model(model$model.file, data = all.data,
                              n.chains = n.chains, n.adapt = n.adapt, quiet = TRUE,
                              inits = inits)
    
    update(jm, n.iter = n.burn.in, progress.bar = "none")   
    rjags::coda.samples(jm, arg.names.params, n.iter = n.iter, progress.bar = "none")   
}

## Generic function for evaluation of expected utility.
eval.eu <- function(model, utility.fun, data,
                    n.iter, n.burn.in, n.adapt = 1000, inits = NULL,
                    independent.SE = FALSE) {
    UseMethod("eval.eu", model)
}

## Method for evaluation of expected utility using simulation.
## The arguments of the utility function must be named as the "parameters" to draw
## or the "data" in the data file of JAGS.
eval.eu.sim.model <- function(model, utility.fun, data,
                              n.iter, n.burn.in, n.adapt = 1000, inits = NULL,
                              independent.SE = FALSE) {       
    ## Get the total set of data and the names
    all.data <- c(model$data, data)
    all.data.names <- names(all.data)

    ## Extract the argument names of the utility function,
    ## and construct separate lists of the param and data names used
    arg.names <- names(formals(utility.fun))    
    arg.names.params <- setdiff(arg.names, all.data.names)
    arg.names.data <- setdiff(arg.names, arg.names.params)   
    
    ## Load jags model from file, update for "burn in", draw coda samples and return the result as a matrix
    jm <-
        if (identical(inits, NULL))
            rjags::jags.model(model$model.file, data = all.data, n.adapt = n.adapt, quiet = TRUE)
        else
            rjags::jags.model(model$model.file, data = all.data, n.adapt = n.adapt, quiet = TRUE, inits = inits)
    
    update(jm, n.iter = n.burn.in, progress.bar = "none")   
    M <- as.matrix(rjags::coda.samples(jm, arg.names.params, n.iter = n.iter, progress.bar = "none"))
    
    ## Get the column names of the coda samples, then remove them
    M.colnames <- colnames(M)
    colnames(M) <- NULL

    ## Create a list of named dimensions to be converted into a named array list of params,
    ## and construct the list of data containing the names used by the utility function.
    named.dims <- named.dims.from.elements(M.colnames)
    data.used <- all.data[arg.names.data]        

    ## Sum over the rows of M (the individual samples)
    us <- vector(mode = "numeric", length = nrow(M))    
    for (i in 1:nrow(M)) {
        named.arrays <- construct.named.arrays(named.dims, M[i,])
        us[i] <- do.call(utility.fun, c(named.arrays, data.used))
    }
    
    ## Return the observed sample mean and the standard error of the sample mean
    if (independent.SE)        
        list(mean = mean(us), SE = sqrt(var(us) / length(us)))
    else
        list(mean = mean(us), SE = sqrt(coda::spectrum0.ar(us)$spec / length(us)))   
}

## Generic function for expected utility evaluation on a grid
eval.on.grid <- function(model, utility.fun, grid.spec.list,
                         n.iter, n.burn.in, n.adapt = 1000,
                         independent.SE = FALSE,
                         parallel = FALSE) {
    UseMethod("eval.on.grid", model)
}

## Evaluate (simulate) on a grid and return a new object containing the results saved in a matrix
eval.on.grid.sim.model <- function(model, utility.fun, grid.spec.list,
                                   n.iter, n.burn.in, n.adapt = 1000,
                                   independent.SE = FALSE,
                                   parallel = FALSE) {        
    ## Construct a grid point list of simulation points from a grid.spec.list
    sim.point.list <- make.grid.point.list(grid.spec.list)       
    grid.size <- length(sim.point.list)

    ## Do parallel or serial simulation using JAGS    
    means.and.SEs <- matrix(vector(mode = "numeric"), nrow = 2, ncol = grid.size)
    
    if (identical(parallel, TRUE)) {
        
        ## Initialise worker processes for parallel simulation
        detected.cores <- parallel::detectCores()
        n.cores <- if (is.na(detected.cores)) 1 else detected.cores
        cl <- parallel::makeCluster(n.cores, "PSOCK")
        
        cl.split <- parallel::clusterSplit(cl, 1:grid.size)
        par.seeds <- rjags::parallel.seeds("base::BaseRNG", grid.size)
        work.list <- lapply(cl.split, function(is) list(par.seeds[is], sim.point.list[is]))
        
        ## Export necessary functions
        parallel::clusterExport(cl, c("eval.eu",
                                      "eval.eu.sim.model",
                                      "construct.named.arrays",
                                      "named.dims.from.elements"),
                                environment())
        
        parallel.eval.eu <- function(chunk) {   
            mapply(
                function(ps, sp) {
                    out <- eval.eu(model, utility.fun, sp,
                                   n.iter, n.burn.in, n.adapt = n.adapt, inits = ps,
                                   independent.SE = independent.SE)
                    c(out$mean, out$SE)                    
                }, chunk[[1]], chunk[[2]])                       
        }   
 
        means.and.SEs <- do.call(cbind, parallel::clusterApply(cl, work.list, parallel.eval.eu))      
        parallel::stopCluster(cl)
        
    } else {
        means.and.SEs <- mapply(
            function(sp) {
                out <- eval.eu(model, utility.fun, sp, n.iter, n.burn.in, n.adapt = n.adapt,
                               independent.SE = independent.SE)
                c(out$mean, out$SE)                
            }, sim.point.list)
    }
    
    ## Return updated model object
    model$sim.points <- gpl.to.matrix(sim.point.list)    
    model$sim.means <- means.and.SEs[1,]
    model$sim.SEs <- means.and.SEs[2,]
    model$grid.spec.list <- grid.spec.list
    model$regression.fun <- NA
    model$gpr.hyper.params <- NA

    model
}

## Plot method for a sim.model -------------------------------------------------------

## Make a plot of the simulation results for a model.
## This function should be called with a model for which some simulation has actually been done.
## main.var.name provides the name (as a string) of a choice variable of main interest.
## Plotting will be done for points in the choice set satisfying mian.var.min <= main.var <= main.var.max.
## fixed is a list of vectors with named entries. Each such vector defines a set of fixed values
## for the remaining choice variables. Hence, the number of curves in the plot will be equal to the length of fixed.
## In the special case when there is only one choice variable, only the model needs to be specified.
## In the special case when there are only two choice variables, fixed may also be given as a vector.
## It then specifies the values of the secondary variable and one curve will be drawn for each value.
## By default, legends are included, with numbers corresponding to the entries of fixed. Set no.legends = TRUE to remove it.
## The default behaviour is to also plot the fitted gpr function if it is available.
plot.sim.model <- function(x, main.var.name = NULL,
                           main.var.min = -Inf,
                           main.var.max = Inf,
                           fixed = list(),
                           no.legends = FALSE,
                           no.reg = FALSE,
                           reg.steps = 100, ...) {
    ## Check that we have something to plot
    if (identical(x$sim.points, NA))
        stop("no simulation points found (evaluate on a grid before plotting)")

    ## Get the names of all decision variables
    var.names <- rownames(x$sim.points)
    
    ## Handle the case when there is only one variable, which is then taken to be the main variable
    if (identical(main.var.name, NULL)) {
        if (length(var.names) != 1)
            stop("a name for the main variable must specified be unless there is exactly one decision variable")
        else
            main.var.name <- var.names
    }
       
    ## Check that main.var.name is valid
    if (!(main.var.name %in% var.names)) {
        stop(main.var.name, " is not a valid main variable name")
    }        
    
    ## Handle the case when there are only two variables and fixed is given as a numeric vector
    if (length(var.names) == 2 && is.numeric(fixed)) {
        other.var <- setdiff(var.names, main.var.name)       
        fixed <- lapply(as.list(fixed), function(y) {names(y) <- other.var; y})
    }
    
    ## Adjust the range
    main.var.min <- max(main.var.min, min(x$sim.points[main.var.name,]))
    main.var.max <- min(main.var.max, max(x$sim.points[main.var.name,]))
    if (main.var.min > main.var.max)
        stop("min value of main var must be less than or equal to the max value")

    ## Check that the names of each element of fixed equals the names of the choice variables
    ## with main.var.name removed
    if (!(all(sapply(fixed, length) == length(var.names) - 1) &&
          all(sapply(fixed, function(y) identical(setdiff(var.names, names(y)), main.var.name))))) {
        stop("the names of each element of 'fixed' must equal the total set of decision variables with the main variable removed")
    }

    ## Construct matrix column selection masks.
    ## The or.mask defines all columns (or points) to be included in the plot. 
    range.mask <- x$sim.points[main.var.name,] >= main.var.min & x$sim.points[main.var.name,] <= main.var.max   
    
    fixed.masks <-
        if (length(var.names) == 1)
            list(range.mask)
        else            
            Map(function(y)
                Reduce('&', c(list(range.mask),
                              Map(function(n) elementwise.all.equal(x$sim.points[n,], as.numeric(y[n])),
                                              names(y)))),
                fixed)

    ## Remove all empty masks (having all false values)
    length.before.filter <- length(fixed.masks)
    fixed.masks <- Filter(any, fixed.masks)

    if (length(fixed.masks) == 0)
        stop("no combination of values in 'fixed' matched a point set")
    
    if (length.before.filter > length(fixed.masks))
        warning("some combination of values in 'fixed' did not match a point set")
   
    or.mask <- Reduce('|', fixed.masks)
    
    ## Set up plot
    y.range.min <- min(x$sim.means[or.mask])
    y.range.max <- max(x$sim.means[or.mask])

    plot(x = c(main.var.min, main.var.max),
         y = c(y.range.min, y.range.max),
         type = "n",
         xlab = main.var.name,
         ylab = "Expected utility")

    if (!no.legends)
        legend(main.var.min,
               y.range.max,
               1:length(fixed.masks),
               pch = 1:length(fixed.masks))
    
    ## Plot selected values as points
    for (i in 1:length(fixed.masks)) {
        points(x$sim.points[main.var.name, fixed.masks[[i]]],
               x$sim.means[fixed.masks[[i]]],
               pch = i)   
    }

    ## Plot fitted gpr function as lines
    if (!identical(x$regression.fun, NA) && !no.reg) {

        if (!no.legends)
            legend(main.var.max - (main.var.max - main.var.min) / 8,
                   y.range.max,
                   1:length(fixed.masks),
                   lty = 1:length(fixed.masks))

        reg.x <- seq(main.var.min, main.var.max, length.out = reg.steps)
        for (i in 1:length(fixed.masks)) {
            M <-
                if (length(var.names) == 1)
                    matrix(reg.x, nrow = 1, dimnames = list(main.var.name))
                else                    
                    do.call(cbind,
                            lapply(reg.x, function(y) {
                                z <- c(y, fixed[[i]]);
                                names(z)[1] <- main.var.name;
                                z}))

            lines(reg.x,
                  apply(M[var.names,, drop = FALSE], 2, x$regression.fun),
                  lty = i)   
        }
    }

    NULL
}    
## Generic function for gpr fitting
fit.gpr <- function(model, start, gr = TRUE, method = "L-BFGS-B", lower = 0, upper = Inf, control = list()) {
    UseMethod("fit.gpr", model)
}

fit.gpr.sim.model <- function(model, start, gr = TRUE, method = "L-BFGS-B", lower = 0, upper = Inf, control = list()) {    
    
    ## Check that a simulation has been performed
    if (identical(model$sim.points, NA))        
        stop("perform a simulation before doing Gaussian process regression")

    ## Check that the method chosen is supported by optim
    if (!(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")))
        stop("the chosen otimisation method is not supported by optim")

    ## Check that lower satisfies lower >= 0
    if (any(lower < 0))
        stop("each element of 'lower' must be >= 0")
    
    ## Check that lower is below upper
    if ( any(lower >= upper) )
        stop("all entries of 'lower' must be strictly below the corresponding entry of 'upper'")    

    if (any(start < lower) || any(start > upper))
        stop("'start' must satisfy 'lower' <= 'start' <= 'upper'")
    
    ## The variances of the simulation errors, obtained by squaring the estimated SEs.    
    vars <- model$sim.SEs^2
    
    ## Proceed according to method to use for optim
    if (identical(method, "L-BFGS-B")) {      
        f <- function(hyper.params)
            squared.exp.lh.obj(model$sim.points, model$sim.means, vars,
                               hyper.params[1], hyper.params[-1])
    } else {        
        f <- function(hyper.params)
            if (any(hyper.params < lower) || any(hyper.params > upper))
                Inf
            else
                squared.exp.lh.obj(model$sim.points, model$sim.means, vars,
                                   hyper.params[1], hyper.params[-1])

        lower <- -Inf
        upper <- Inf
    }    
    
    ## Check if a gradient should be used
    gr <-
        if (identical(gr, TRUE) && !identical(method, "SANN"))
            function(hyper.params)
                squared.exp.lh.obj.grad(model$sim.points, model$sim.means, vars, hyper.params[1], hyper.params[-1])
        else
            NULL

    ## Fit the hyperparameters to the observations by calling the standard optim routine
    optim.out <- optim(start, f, gr = gr, method = method, lower = lower, upper = upper, control = control)
    
    ## Check and report on output of optim
    if (!identical(optim.out$convergence, 0)) {
        if(identical(optim.out$convergence, 1))            
            warning("iteration limit 'maxit' reached in optim")
        
        if(identical(optim.out$convergence, 10))
            warning("degeneracy of the Nelder-Mead simplex in optim")

        if(identical(optim.out$convergence, 51) || identical(optim.out$convergence, 52))
            warning(optim.out$message)
    }

    ## Return updated model object       
    model$gpr.hyper.params <- optim.out$par
    model$regression.fun <- gpr.given.hyperparams(model$sim.points, model$sim.means, vars, model$gpr.hyper.params[1], model$gpr.hyper.params[-1])    
    
    model
}

## Generic function for loess fitting
fit.loess <- function(model, span = 0.75, degree = 2) {
    UseMethod("fit.loess", model)
}

## Fitting a loess model to the observed simulation values
fit.loess.sim.model <- function(model, span = 0.75, degree = 2) {        
    
    ## Check that a simulation has been performed
    if (identical(model$sim.points, NA))        
        stop("perform a simulation before attempting regression")
    
    ## Get the number of variables
    n.vars <- nrow(model$sim.points)

    ## Check that the number of variables is accepted by loess
    if (n.vars < 1 || n.vars > 4)
        stop("the number of variables for loess regression must be >= 1 and <= 4")
       
    ## Create formula object for loess regression
    xnam <- paste0("x", 1:n.vars)
    fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))

    ## Create data frame form the grid points and simulated means
    M <- t(model$sim.points)
    colnames(M) <- xnam   
    df <- data.frame(cbind(y = model$sim.means, M))

    ## Perform loess regression
    loess.out <- loess(fmla, df, span = span, degree = degree, surface = "direct")   
    
    ## Return updated model     
    model$regression.fun <- function(x) {        
        names(x) <- NULL
        predict(loess.out, matrix(x, 1, length(x)))
    }

    model$gpr.hyper.params <- NA
    
    model
}


## Generic function for expected utility optimisation
optimise.eu <- function(model, start, method = "L-BFGS-B", lower = -Inf, upper = Inf, control = list()) {
    UseMethod("optimise.eu", model)
}

## Optimisation of expected utility.
## method is either "Grid" or one of the methods supported by optim.
## If "Grid", no argument is used except model.
optimise.eu.sim.model <- function(model, start, method = "L-BFGS-B", lower = -Inf, upper = Inf, control = list()) {
    ## Check that the method chosen is supported by optim or is "Grid"
    if (!(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "Grid")))
        stop("the chosen otimisation method is not supported")    

    ## Check that model has been evaluated over a grid
    if (identical(model$sim.points, NA))
        stop("simulate on a grid before optimising")
    
    if (identical(method, "Grid")) {
        ## Find column index of grid point corresponding to max utility
        max.index <- which.max(model$sim.means)
        
        opt.arg <- model$sim.points[,max.index]
        opt.eu <- model$sim.means[max.index]

        return(list(opt.arg = opt.arg, opt.eu = opt.eu))
    }

    if (identical(model$regression.fun, NA))
        stop("regression function not available (simulate and fit before optimising)")

    ## Check that lower is below upper
    if ( any(lower >= upper) )
        stop("all entries of 'lower' must be strictly below the corresponding entry of 'upper'")    

    if (any(start < lower) || any(start > upper))
        stop("'start' must satisfy 'lower' <= 'start' <= 'upper'")        

    ## Proceed according to method to use for optim
    if (identical(method, "L-BFGS-B")) {      
        f <- model$regression.fun          
    } else {        
        f <- function(x)
            if (any(x < lower) || any(x > upper))
                -Inf
            else
                model$regression.fun(x)

        lower <- -Inf
        upper <- Inf
    }    
    
    ## Optimise the fitted regression function.
    ## Any control values supplied by the user is appended to (and may override!) fnscale = -1.
    optim.out <- optim(start, model$regression.fun, gr = NULL, method = method,
                       lower = lower, upper = upper, control = c(list(fnscale = -1), control))
        
    opt.arg <- optim.out$par
    opt.eu <- optim.out$value
    names(opt.arg) <- rownames(model$sim.points)
        
    return(list(opt.arg = opt.arg, opt.eu = opt.eu))
}
