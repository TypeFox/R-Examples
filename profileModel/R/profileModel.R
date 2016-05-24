`prelim.profiling` <-
function (fitted, quantile = qchisq(0.95, 1), objective = stop("'objective' is missing."),
    verbose = TRUE, which = 1:length(coef(fitted)), stepsize = 0.5,
    stdn = 5, agreement = TRUE, trace.prelim = FALSE, stdErrors = NULL,
    ...)
    ## which should be a vector of integers
{
    if (is.null(stdErrors))
        stdErrors <- summary(fitted)$coefficients[, 2]
    Betas <- coef(fitted)
    BetasNames <- names(Betas)
    noNA <- !is.na(Betas) # aliased out
    fm.call <- fitted$call
    mf <- model.frame(fitted$terms,data=eval(fm.call$data))
    if (inherits(fitted,"BTm")) {
        Y <- fitted$y0
        X <- fitted$x0
    }
    else {
        Y <- model.response(mf)
        X <- model.matrix(fitted$terms,mf,contrasts = fitted$contrasts)
    }
    if (is.null(fitted$X.max.scaleFit) & inherits(fitted,"polr"))
        X <- X[, -1, drop = FALSE]
    O <- model.offset(mf)
    if (is.null(O))
        O <- rep(0, nrow(X))
#    fitted.formula <- formula(fitted)
    fm.call$offset <- NULL
    LP.or <- fitted$linear.predictor
    if (length(Betas) == 1)
        fm.call$formula <- Y ~ -1 + offset(o)
    else {
      if (inherits(fitted,"BTm")) {
        fm.call$formula <- Y ~ -1 + offset(o) + Xnoi
      }
      else {
        XnoiNams <- paste(paste("Xnoi[,", 1:(ncol(X) - sum(!noNA) -1),
                                "]", sep=""))
        fm.call$formula <- as.formula(paste("Y ~ -1 + offset(o) + ",
                                      paste(XnoiNams, collapse = "+")))
      }
    }
    fm.call$offset <- NULL
    fitted.formals <- names(formals(as.character(fm.call)))
    test1 <- "etastart" %in% fitted.formals
    if (test1) {
        LP.or <- fitted$linear.predictor
        fm.call$etastart <- as.symbol("LP")
    }
    ObjValue.or <- objective(fitted, ...)
    p <- length(which)
    grid.bounds <- intersects <- matrix(NA, p, 2, dimnames = list(BetasNames[which],
        c("Left", "Right")))
    ## The information in prelim.profile could be used... TBD
    ##  if (return.profiles) {
    ##    res <- as.list(rep(NA,p))
    ##    names(res) <- BetasNames[which]
    ##  }
    should.intersect <- intersects.temp <- rep(NA, 2)
    numberofsteps <- stdn/stepsize
    if (verbose & !trace.prelim)
        cat("Preliminary iteration ")
    for (i in which) {
        if (verbose & !trace.prelim)
            cat(".")
        if (!noNA[i])  # aliased out
            next
        tb.included <- noNA
        tb.included[i] <- FALSE # without the aliased
        profiledName <- BetasNames[i]
        if (trace.prelim)
            cat(profiledName, "\n")
        stepsize.temp <- c(stepsize, stepsize)
        stdErrors.i <- stdErrors[i]
        Xnoi <- X[, tb.included, drop = FALSE] # without the aliased
        Xonlyi <- X[, i]
        ################
        ### left
        ################
        cc <- 1
        tempDiff <- quantile + 10
        test.intersections <- TRUE
        b <- Betas[i]
        b <- b.old <- sign(b) * min(abs(b), 30)
        LP <- (abs(b) <= 30) * LP.or
        while (test.intersections) {
            b <- b.old
            ## if (trace.prelim) {
            ##   title1 <- paste(profiledName,"left","with stepsize",stepsize.temp[1]) #TR
            ##   plot(1,1,xlim=c(Betas[i]-5*min(30,stdErrors.i),Betas[i]+5*min(30,stdErrors.i)), #TR
            ##        ylim=c(0,100),type="n",main=title1) #TR
            ##   points(x=c(Betas[i]-5*min(30,stdErrors.i),Betas[i]+5*min(30,stdErrors.i)), #TR
            ##          y=c(quantile,quantile),type="l") #TR
            ## }
            test <- TRUE
            curPoint <- 0
            slope.pp <- 1
            slopes.pp <- rep(1, numberofsteps + 5)
            while (test & curPoint < numberofsteps) {
                tempDiff.old <- tempDiff
                b.old <- b
                ## if (trace.prelim) {
                ##  if (!curPoint) {points(x=b,y=tempDiff,pch='s');Sys.sleep(1)} #TR
                ##  else {points(x=b,y=tempDiff);Sys.sleep(1)} #TR
                ## }
                curPoint <- curPoint + 1
                b <- b - min(30, stdErrors.i)/slope.pp * curPoint *
                  stepsize.temp[1]
                o <- O + Xonlyi * b
                suppressWarnings(fm <- eval(fm.call))
                LP <- fm$linear.predictor
                tempDiff <- objective(fm, ...) - ObjValue.or
                if (is.na(tempDiff))
                  stop("Profiling failed. NA's introduced by the objective.")
                if (is.infinite(tempDiff)) {
                  warning("Infinite values were introduced by the objective.")
                  slope.pp <- 1
                }
                else {
                  slope.pp <- abs(ss <- (tempDiff.old - tempDiff)/(b.old -
                    b))
                  slopes.pp[curPoint] <- slope.pp
                  ## set to give the objective a chance to increase and at the
                   # same time  to avoid huge steps while being conservative
                  if (slope.pp < 1)
                    slope.pp <- 1
                  if (slope.pp > 500)
                    slope.pp <- 500
                }
                ## if you have done the first three iterations and nothing
                 # is found stop
                if (curPoint < 4)
                  nonzero.slopes <- TRUE
                else nonzero.slopes <- !all(slopes.pp[(curPoint -
                  3):curPoint] < 1e-08)
                test <- (tempDiff < quantile | ss > 1e-08) &
                  nonzero.slopes
                if (trace.prelim)
                  cat("<-- iteration:", curPoint, "\t", paste(c("CPV:",
                    "SL:", "OV:"), format(round(c(b, ss, tempDiff),
                    digits = 3), zero.print = TRUE)), "SS:",
                    stepsize.temp[1], "\n")
            }
            grid.bounds[profiledName, 1] <- b
            should.intersect[1] <- nonzero.slopes
            intersects.temp[1] <- tempDiff > quantile
            ## test.intersections is interpreted as should intersect but it does not
             # ss>0 stands for the case where the search on the left side starts outside
             # the right end of the profiled objective (if agreement is false then this is
             # necessary. if agreement is true then ss>0 does not violate the facts)
            test.intersections <- (should.intersect[1] | ss >
                1e-08) & (!intersects.temp[1] | ss > 1e-08)
        }
        ################
        ### right
        ################
        cc <- 1
        tempDiff <- quantile + 10
        test.intersections <- TRUE
        b <- Betas[i]
        b <- b.old <- sign(b) * min(abs(b), 30)
        LP <- (abs(b) <= 30) * LP.or
        while (test.intersections) {
            b <- b.old
            ## if (trace.prelim) {
            ## title1 <- paste(profiledName,"right","with stepsize",stepsize.temp[2]) #TR
            ##  plot(1,1,xlim=c(Betas[i]-5*min(30,stdErrors.i),Betas[i]+5*min(30,stdErrors.i)), #TR
            ##       ylim=c(0,200),type="n",main=title1) #TR
            ##  points(x=c(Betas[i]-5*min(30,stdErrors.i),Betas[i]+5*min(30,stdErrors.i)), #TR
            ##        y=c(quantile,quantile),type="l") #TR
            ## }
            test <- TRUE
            curPoint <- 0
            slope.pp <- 1
            slopes.pp <- rep(1, numberofsteps + 5)
            while (test & curPoint < numberofsteps) {
                tempDiff.old <- tempDiff
                b.old <- b
                ## if (trace.prelim) {
                ##  if (!curPoint) {points(x=b,y=tempDiff,pch='s');Sys.sleep(1)} #TR
                ##  else {points(x=b,y=tempDiff);Sys.sleep(1)} #TR
                ## }
                curPoint <- curPoint + 1
                b <- b + min(30, stdErrors.i)/slope.pp * curPoint *
                  stepsize.temp[2]
                o <- O + Xonlyi * b
                suppressWarnings(fm <- eval(fm.call))
                LP <- fm$linear.predictor
                tempDiff <- objective(fm, ...) - ObjValue.or
                if (is.na(tempDiff))
                  stop("Profiling failed. NA's introduced by the objective.")
                if (is.infinite(tempDiff)) {
                  warning("Infinite values were introduced by the objective.")
                  slope.pp <- 1
                }
                else {
                  slope.pp <- abs(ss <- (tempDiff.old - tempDiff)/(b.old -
                    b))
                  slopes.pp[curPoint] <- slope.pp
                  ## set to give the objective a chance to increase and at the
                   # same time  to avoid huge steps while being conservative
                  if (slope.pp < 1)
                    slope.pp <- 1
                  if (slope.pp > 500)
                    slope.pp <- 500
                }
                ## if you have done the first three iterations and nothing
                # is found stop
                if (curPoint < 4)
                  nonzero.slopes <- TRUE
                else nonzero.slopes <- !all(slopes.pp[(curPoint -
                  3):curPoint] < 1e-08)
                test <- (tempDiff < quantile | ss < -1e-08) &
                  nonzero.slopes
                if (trace.prelim)
                  cat("--> iteration:", curPoint, "\t", paste(c("CPV:",
                    "SL:", "OV:"), format(round(c(b, ss, tempDiff),
                    digits = 3), zero.print = TRUE)), "SS:",
                    stepsize.temp[1], "\n")
            }
            grid.bounds[profiledName, 2] <- b
            should.intersect[2] <- nonzero.slopes
            intersects.temp[2] <- tempDiff > quantile
            ## test.intersections is interpreted as should intersect but it does not
             # ss>0 stands for the case where the search on the left side starts outside
             # the right end of the profiled objective (if agreement is false then this is
             # necessary. if agreement is true then ss>0 does not violate the facts)
            test.intersections <- (should.intersect[2] | ss <
                -1e-08) & (!intersects.temp[2] | ss < -1e-08)
            stepsize.temp[2] <- stepsize.temp[2] + stepsize
            cc <- cc + 1
        }
        intersects[profiledName, ] <- intersects.temp
        if (i == which[p] & verbose & !trace.prelim)
            cat(" Done\n\n")
    }
    if (trace.prelim) {
        cat("<--: Left | -->: Right | CPV: Current Parameter value\n")
        cat("SL: slope | OV: Objective Value | SS: StepSize\n")
    }
    list(grid.bounds = grid.bounds, intersects = intersects)
}

`profileModel` <-
function (fitted, gridsize = 20, stdn = 5, stepsize = 0.5, grid.bounds = NULL,
    quantile = NULL, objective = stop("'objective' is missing."),
    agreement = TRUE, verbose = TRUE, trace.prelim = FALSE, which = 1:length(coef(fitted)),
    profTraces = TRUE, zero.bound = 1e-08, scale = FALSE, stdErrors = NULL,
    ...)
{
    Betas <- coef(fitted)
    BetasNames <- names(Betas)
    noNA <- !is.na(Betas)
    if (is.null(stdErrors)) {
      stdErrors <- rep(NA, length(Betas))
      stdErrors[noNA] <- summary(fitted)$coefficients[BetasNames[noNA], 2]
    }
    if (scale) {
        fitted <- scaleFit(fitted)
        Xmax <- fitted$X.max.scaleFit
    }
    if ((zero.bound < 0) | (zero.bound > 1e-06)) {
        stop("zero.bound takes values between 0 and 1e-6.")
    }
    if (is.character(which)) {
        which <- match(which, BetasNames)
        ttt <- is.na(which)
        if (any(ttt))
            stop("A least a parameter name specified in 'which' does not exist in the fitted model.")
    }
    if (any(duplicated(which))) {
        warning("At least a parameter was specified more than once in 'which'. Profiling for the duplicated parameter(s) was done only once.")
        which <- unique(which)
    }
    if (min(which) < 1 | max(which) > length(Betas)) {
        stop("At least a parameter position specified in 'which' is not valid.")
    }
    p <- length(which)
    na.in.which <- !noNA[which]
    if (all(na.in.which)) {
        stop("'which' refers to parameters which have value 'NA' in the original fit.")
    }
    if (any(na.in.which)) {
        warning("At least a parameter with value 'NA' exists in the original fit.  Profiling did not take place for these parameters.")
    }
    if (!is.null(grid.bounds))
        if (length(grid.bounds) != 2 * p)
            stop("The dimension of 'grid.bounds' is not compatible with the length of 'which'.")
    objective <- match.fun(objective)
    if (is.null(grid.bounds)) {
        if (is.null(quantile)) {
            grid.bounds <- cbind(Betas[which] - stdn * stdErrors[which],
                Betas[which] + stdn * stdErrors[which])
            if (scale)
                grid.bounds <- grid.bounds * Xmax[which]
            result <- profiling(fitted, grid.bounds = grid.bounds,
                gridsize = gridsize, verbose = verbose, objective = objective,
                which = which, agreement = agreement, profTraces = profTraces,
                zero.bound = zero.bound, ...)
            intersects <- NULL
            attr(grid.bounds, "from.prelim") <- FALSE
        }
        else {
           if (scale)
                stdErrors <- stdErrors * Xmax
            prelim.res <- prelim.profiling(fitted, quantile = quantile,
                objective = objective, verbose = verbose, which = which,
                stepsize = stepsize, stdn = stdn, agreement = agreement,
                trace.prelim = trace.prelim, stdErrors = stdErrors,
                ...)
            grid.bounds <- prelim.res$grid.bounds
            result <- profiling(fitted, grid.bounds = grid.bounds,
                gridsize = gridsize, verbose = verbose, objective = objective,
                which = which, agreement = agreement, profTraces = profTraces,
                zero.bound = zero.bound, ...)
            intersects <- prelim.res$intersects
            rownames(intersects) <- BetasNames[which]
            attr(grid.bounds, "from.prelim") <- TRUE
        }
    }
    else {
        if (is.null(dim(grid.bounds)))
            grid.bounds <- matrix(grid.bounds, ncol = 2, byrow = TRUE)
        if (scale)
            grid.bounds <- grid.bounds * Xmax[which]
        result <- profiling(fitted, grid.bounds = grid.bounds,
            gridsize = gridsize, verbose = verbose, objective = objective,
            which = which, agreement = agreement, profTraces = profTraces,
            zero.bound = zero.bound, ...)
        intersects <- NULL
        attr(grid.bounds, "from.prelim") <- FALSE
    }
    names(result) <- rownames(grid.bounds) <- BetasNames[which]
    if (scale) {
        grid.bounds <- grid.bounds/Xmax[which]
        for (i in 1:p) {
            if (!noNA[which[i]])
                next
            result[[i]][, 1] <- result[[i]][, 1]/Xmax[which[i]]
            colnames(result[[i]])[1] <- BetasNames[which[i]]
        }
        if (profTraces) {
            for (i in 1:p) {
                if (!noNA[which[i]])
                  next
                tb.included <- noNA
                tb.included[which[i]] <- FALSE
                result[[i]][, -c(1, 2)] <- sweep(result[[i]][,
                  -c(1, 2), drop = FALSE], 2, Xmax[tb.included], "/")
                colnames(result[[i]])[-c(1, 2)] <- BetasNames[tb.included]
            }
        }
    }
    dotss <- match.call(expand.dots = FALSE)[["..."]]
    dotssNames <- names(dotss)
    for (i in dotssNames) formals(objective)[[i]] <- eval(dotss[[i]])
    result <- list(profiles = result, fit = fitted, quantile = quantile,
        gridsize = gridsize, intersects = intersects, profiled.parameters = which,
        profiled.objective = objective, isNA = !noNA[which],
        agreement = agreement, zero.bound = zero.bound, call = match.call(),
        grid.bounds = grid.bounds)
    attr(result, "includes.traces") <- profTraces
    class(result) <- "profileModel"
    result
}


`profiling` <-
function (fitted, grid.bounds, gridsize = 20, verbose = TRUE,
    objective = stop("'objective' is missing."), agreement = TRUE,
    which = 1:length(coef(fitted)), profTraces = TRUE, zero.bound = 1e-08,
    ...)
    ## which should be a vector of integers
    ## grid.bounds should be a 2*length(which) vector of reals or
     # a 2 by length(which) matrix of reals
{
    if (is.null(dim(grid.bounds)))
        grid.bounds <- matrix(grid.bounds, ncol = 2, byrow = TRUE)
    Betas <- coef(fitted)
    p.or <- length(Betas)
    BetasNames <- names(Betas)
    noNA <- !is.na(Betas)
    p <- length(which)
    fm.call <- fitted$call
    mf <- model.frame(fitted$terms,data=eval(fm.call$data))
    Y <- model.response(mf)
    if (inherits(fitted,"BTm")) {
        Y <- fitted$y0
        X <- fitted$x0
    }
    else {
        Y <- model.response(mf)
        X <- model.matrix(fitted$terms,mf,contrasts = fitted$contrasts)
    }
    if (is.null(fitted$X.max.scaleFit) & inherits(fitted, "polr"))
        X <- X[, -1, drop = FALSE]
    O <- model.offset(mf)
    if (is.null(O))
        O <- rep(0, nrow(X))
#    fitted.formula <- formula(fitted)
    if (p.or == 1)
        fm.call$formula <- Y ~ -1 + offset(o)
    else {
      if (inherits(fitted,"BTm"))
        fm.call$formula <- Y ~ -1 + offset(o) + Xnoi
      else {
        XnoiNams <- paste(paste("Xnoi[,", 1:(ncol(X) - sum(!noNA) -1),
                                "]", sep=""))
        fm.call$formula <- as.formula(paste("Y ~ -1 + offset(o) + ",
                                      paste(XnoiNams, collapse = "+")))
      }
    }
    fm.call$offset <- NULL
    fitted.formals <- names(formals(as.character(fm.call)))
    test1 <- "etastart" %in% fitted.formals
    if (test1) {
        LP.or <- fitted$linear.predictor
        fm.call$etastart <- as.symbol("LP")
    }
    ObjValue.or <- objective(fitted, ...)
    res <- as.list(rep(NA, p))
    names(res) <- BetasNames[which]
    for (i in 1:p) {
        iprof <- which[i]
        if (!noNA[iprof]) # aliased out
            next
        tb.included <- noNA
        tb.included[iprof] <- FALSE # without the aliased
        profiledName <- BetasNames[iprof]
        gridd <- seq(grid.bounds[i, 1], grid.bounds[i, 2], length = gridsize)
        curPoint <- 0
        Xnoi <- X[, tb.included, drop = FALSE] # without the aliased
        Xonlyi <- X[, iprof]
        inds.right <- which(gridd >= Betas[iprof])
        inds.left <- which(gridd < Betas[iprof])
        # Make sure you start as close as possible to the estimate
        if (grid.bounds[i, 1] <= grid.bounds[i, 2])
            inds.left <- inds.left[order(inds.left, decreasing = TRUE)]
        else inds.right <- inds.right[order(inds.right, decreasing = TRUE)]
        inds <- list(inds.left, inds.right)
        ObjValues <- cbind(gridd, 0)
        if (profTraces) {
            tracesNames <- BetasNames[tb.included]
            profile.traces <- matrix(0, nrow = gridsize, ncol = sum(noNA) -
                1)
            colnames(profile.traces) <- tracesNames
        }
        colnames(ObjValues) <- c(profiledName, "Differences")
        if (verbose)
            cat("Profiling for parameter", profiledName, "...")
        for (k in 1:2) {
            if (test1)
                LP <- LP.or
            ## else supply no starting valiues...
            ## maybe an argument to control starting values??? TBD
            for (curPoint in inds[[k]]) {
                bp <- gridd[curPoint]
                o <- O + Xonlyi * bp
                suppressWarnings(fm <- eval(fm.call))
                ## LP will be NULL if fm$linear.predictor does not exist...OK
                LP <- fm$linear.predictor
                ObjValue.current <- (objective(fm, ...) - ObjValue.or)
                if (is.na(ObjValue.current))
                  stop("Profiling failed. NA's introduced by the objective.")
                if (is.infinite(ObjValue.current))
                  warning("Infinite values were introduced by the objective.")
                if (agreement) {
                  if (ObjValue.current < -(zero.bound * 1000)) {
                    stop("Profiling has found a better solution. Original fit had not converged.")
                  }
                  ObjValues[curPoint, 2] <- (ObjValue.current >=
                    zero.bound) * ObjValue.current
                }
                else {
                  ObjValues[curPoint, 2] <- ObjValue.current
                }
                if (profTraces)
                  profile.traces[curPoint, ] <- coef(fm)
            }
        }
        if (verbose)
            cat(" Done\n")
        res[[profiledName]] <- if (profTraces)
            cbind(ObjValues, profile.traces)
        else ObjValues
    }
    res
}
