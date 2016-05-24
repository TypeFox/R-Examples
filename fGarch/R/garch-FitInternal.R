
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                DESCRIPTION:
#  .garchFit                Internal GARCH Fit
#  .garchArgsParser         Parses formula and data for garchFit
#  .garchModelSeries        Composes model series like in lm fits
#  .garchOptimizerControl   Sets default values for Garch Optimizer
#  .garchNames              Slot names, @fit slot, parameters and controls
################################################################################


.garchFit <-
function(
    formula.mean = ~arma(0, 0),
    formula.var = ~garch(1, 1),
    series,
    init.rec = c("mci", "uev"),
    delta = 2,
    skew = 1,
    shape = 4,
    cond.dist = c("norm", "snorm", "ged", "sged", "std", "sstd", "QMLE"),
    include.mean = TRUE,
    include.delta = NULL,
    include.skew = NULL,
    include.shape = NULL,
    leverage = NULL,
    trace = TRUE,
    algorithm = c("sqp", "nlminb", "lbfgsb", "nlminb+nm", "lbfgsb+nm"),
    hessian = c("ropt", "rcd"),
    robust.cvar,
    control = list(),
    title = NULL,
    description = NULL,
    ...)
{
    # A function implemented by Diethelm Wuertz

    # Description
    #   Fit parameters to a ARMA-GARCH model

    # Arguments:
    #   formula.mean - ARMA(m,n) mean specification
    #   formula.var - GARCH/APARCH(p,q) variance specification
    #   series - time series
    #   init.rec - names type of initialization of recurrence
    #       mci = mu-current-iteration, or
    #       uev = unconditional-expected-variances
    #   delta - numeric value of the exponent delta
    #   skew - optional skewness parameter
    #   shape - optional shape parameter
    #   cond.dist - name of the conditional distribution
    #   include.mean - should the mean value be estimated ?
    #   include.delta - should the exponent be estimated ?
    #   leverage - should the leverage factors be estimated ?
    #   trace - should the optimization be traced ?
    #   algorithm -
    #   control - list of additional control parameters for solver
    #   title - an optional title string
    #   description - an optional project description string

    # Note:
    #   This is the old version of garchFit, we keep it for backward
    #   compatibility.

    # FUNCTION:

    # Debug Mode:
    DEBUG <- FALSE

    # Allow only full formula specification:
    if(DEBUG) print("Formula Specification ...")
    fcheck = rev(all.names(formula.mean))[1]
    if (fcheck == "ma") {
        stop("Use full formula: arma(0,q) for ma(q)")
    } else if (fcheck == "ar") {
        stop("Use full formula expression: arma(p,0) for ar(p)")
    }

    # Check for Recursion Initialization:
    if(DEBUG) print("Recursion Initialization ...")
    if(init.rec[1] != "mci" & algorithm[1] != "sqp") {
        stop("Algorithm only supported for mci Recursion")
    }

    # Get Start Time:
    .StartFit <- Sys.time()

    # Generate Control List - Define Default Settings:
    if(DEBUG) print("Generate Control List ...")
    con <- .garchOptimizerControl(algorithm, cond.dist)
    con[(namc <- names(control))] <- control

    # Initialize Time Series Information - Save Globally:
    # keep copy of input data
    if(DEBUG) print("Initialize Time Series ...")
    data <- series
    # scale time series
    scale <- if (con$xscale) sd(series) else 1
    series <- series/scale
    .series <- .garchInitSeries(
        formula.mean = formula.mean,
        formula.var = formula.var,
        cond.dist = cond.dist[1],
        series = series,
        scale = scale,
        init.rec = init.rec[1],
        h.start = NULL,
        llh.start = NULL,
        trace = trace)
    .setfGarchEnv(.series = .series)

    # Initialize Model Parameters - Save Globally:
    if(DEBUG) print("Initialize Model Parameters ...")
    .params <- .garchInitParameters(
        formula.mean = formula.mean,
        formula.var = formula.var,
        delta = delta,
        skew = skew,
        shape = shape,
        cond.dist = cond.dist[1],
        include.mean = include.mean,
        include.delta = include.delta,
        include.skew = include.skew,
        include.shape = include.shape,
        leverage = leverage,
        algorithm = algorithm[1],
        control = con,
        trace = trace)
    .setfGarchEnv(.params = .params)

    # Select Conditional Distribution Function:
    if(DEBUG) print("Select Conditional Distribution ...")
    .setfGarchEnv(.garchDist = .garchSetCondDist(cond.dist[1]))

    # Estimate Model Parameters - Minimize llh, start from big value:
    if(DEBUG) print("Estimate Model Parameters ...")
    .setfGarchEnv(.llh = 1.0e99)
    .llh <- .getfGarchEnv(".llh")
    fit = .garchOptimizeLLH(hessian, robust.cvar, trace)

    # Add to Fit:
    if (DEBUG) print("Add to fit ...")
    .series <- .getfGarchEnv(".series")
    .params <- .getfGarchEnv(".params")
    names(.series$h) <- NULL
    fit$series = .series
    fit$params = .params

    # Retrieve Residuals and Fitted Values:
    if (DEBUG) print("Retrieve Residuals and Fitted Values ...")
    residuals = .series$z
    fitted.values = .series$x - residuals
    h.t = .series$h
    if (.params$includes["delta"])
        deltainv = 1/fit$par["delta"]
    else
        deltainv = 1/fit$params$delta
    sigma.t = (.series$h)^deltainv

    # Standard Errors and t-Values:
    if (DEBUG) print("Standard Errors and t-Values ...")
    fit$cvar <-
        if (robust.cvar)
            (solve(fit$hessian) %*% (t(fit$gradient) %*% fit$gradient) %*%
             solve(fit$hessian))
        else
            - solve(fit$hessian)
    fit$se.coef = sqrt(diag(fit$cvar))
    fit$tval = fit$coef/fit$se.coef
    fit$matcoef = cbind(fit$coef, fit$se.coef,
    fit$tval, 2*(1-pnorm(abs(fit$tval))))
    dimnames(fit$matcoef) = list(names(fit$tval), c(" Estimate",
            " Std. Error", " t value", "Pr(>|t|)"))

    # Add Title and Description:
    if (DEBUG) print("Add Title and Description ...")
    if(is.null(title)) title = "GARCH Modelling"
    if(is.null(description)) description = description()

    # Total Execution Time:
    Time =  Sys.time() - .StartFit
    if(trace) {
        cat("\nTime to Estimate Parameters:\n ")
        print(Time)
    }

    # Return Value:
    new("fGARCH",
        call = as.call(match.call()),
        formula = as.formula(paste("~", formula.mean, "+", formula.var)),
        method = "Max Log-Likelihood Estimation",
        data = data,
        fit = fit,
        residuals = residuals,
        fitted = fitted.values,
        h.t = h.t,
        sigma.t = as.vector(sigma.t),
        title = as.character(title),
        description = as.character(description)
    )
}


# ------------------------------------------------------------------------------


.garchArgsParser <-
function(formula, data, trace = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Parses formula and data for garchFit

    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance
    #       specification
    #   data - time series input as a timeSeries

    # Note:
    #   This function returns the input formula and input data in
    #   proper formats. Two cases are deistinguished

    # FUNCTION:

    # Get Data:
    allVars = unique(sort(all.vars(formula)))
    allVarsTest =  mean(allVars %in% colnames(data))
    if (allVarsTest != 1) {
        print(allVars)
        print(colnames(data))
        stop ("Formula and data units do not match.")
    }
    formula.lhs = as.character(formula)[2]

    # Model frame:
    mf = match.call(expand.dots = FALSE)
    if(trace) {
        cat("\nMatched Function Call:\n ")
        print(mf)
    }
    m = match(c("formula", "data"), names(mf), 0)
    mf = mf[c(1, m)]

    # Model the timeSeries - Have a look on the function .garchModelSeries() ...
    #   here we cant use "model/frame" !
    mf[[1]] = as.name(".garchModelSeries")
    mf$fake = FALSE
    mf$lhs = TRUE
    if(trace) {
        cat("\nModelSeries Call:\n ")
        print(mf)
    }
    x = eval(mf, parent.frame())
    if(trace) print(x)

    # Now extract the modelled series ...
    x = as.vector(x[, 1])
    names(x) = rownames(data)
    if(trace) print(x)

    # Compose Mean and Variance Formula:
    allLabels = attr(terms(formula), "term.labels")
    if(trace) {
        cat("\nAll Term Labels:\n ")
        print(allLabels)
    }
    if(length(allLabels) == 2) {
        formula.mean = as.formula(paste("~", allLabels[1]))
        formula.var = as.formula(paste("~", allLabels[2]))
    } else if(length(allLabels) == 1) {
        formula.mean = as.formula("~ arma(0, 0)")
        formula.var = as.formula(paste("~", allLabels[1]))
    }
    if(trace) {
        cat("\nMean Formula:\n ")
        print(formula.mean)
        cat("\nVariance Formula:\n ")
        print(formula.var)
    }

    # Result:
    ans <- list(
        formula.mean = formula.mean,
        formula.var = formula.var,
        formula.lhs = formula.lhs,
        series = x)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.garchModelSeries <-
function (formula, data, fake = FALSE, lhs = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Composes model series like in lm fits

    # Arguments:

    # Note:
    #   ... have a look on model.frame()

    # FUNCTION:

    # Formula:
    if (length(formula) == 2) {
        formula = as.formula(paste("x", formula[1], formula[2],
            collapse = ""))
        stopifnot(!missing(data))
    }

    # Missing Data ?
    if (missing(data)) {
        data = eval(parse(text = search()[2]), parent.frame())
    }

    # Numeric Data ?
    if (is.numeric(data)) {
        data = data.frame(data)
        colnames(data) = all.vars(formula)[1]
        lhs = TRUE
    }

    # Faked Formula ?
    if (fake) {
        response = as.character(formula)[2]
        Call = as.character(match.call()[[2]][[3]])
        method = Call[1]
        predictors = Call[2]
        formula = as.formula(paste(response, "~", predictors))
    }

    # Left-Hand-Side Formula ?
    if (lhs) {
        response = as.character(formula)[2]
        formula = as.formula(paste(response, "~", 1))
    }

    # Compose Model Frame:
    x = model.frame(formula, data)

    # timeSeries ?
    if (class(data) == "timeSeries") {
        x = timeSeries(x)
    }

    # Add control atrribute:
    if (fake) {
        attr(x, "control") <- method
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


.garchOptimizerControl <-
function(algorithm, cond.dist)
# function(algorithm, cond.dist)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets default values for Garch Optimizer

    # Arguments:
    #   none

    # FUNCTION:

    # Check llh for the standardized NIG Distribution:
    llh = "internal"
    if (cond.dist == "snig") llh = "filter"

    # Generate Control List with Default Settings:
    con <- list(

        # In General:
        fscale = TRUE,
        xscale = TRUE,
        algorithm = algorithm,
        llh = llh,

        # BFGS - NLMINB Algorithm:
        tol1 = 1,
        tol2 = 1,

        # SQP Algorithm:
        MIT = 2000,    # maximum number of iterations (200)
        MFV = 5000,    # maximum number of function evaluations (500)
        MET = 5,       # specifies scaling strategy:
                       #  MET=1 - no scaling
                       #  MET=2 - preliminary scaling in 1st iteration (default)
                       #  MET=3 - controlled scaling
                       #  MET=4 - interval scaling
                       #  MET=5 - permanent scaling in all iterations
        MEC = 2,       # correction for negative curvature:
                       #  MEC=1 - no correction
                       #  MEC=2 - Powell correction (default)
        MER = 1,       # restarts after unsuccessful variable metric updates:
                       #  MER=0 - no restarts
                       #  MER=1 - standard restart
        MES = 4,       # interpolation method selection in a line search:
                       #  MES=1 - bisection
                       #  MES=2 - two point quadratic interpolation
                       #  MES=3 - three point quadratic interpolation
                       #  MES=4 - three point cubic interpolation (default)
        XMAX = 1.0e3,
        TOLX = 1.0e-10,
        TOLC = 1.0e-6,
        TOLG = 1.0e-6,
        TOLD = 1.0e-6,
        TOLS = 1.0e-4,
        RPF  = 1.0e-2) # 1.0e-4)

    # Return Value:
    con
}


# ------------------------------------------------------------------------------


.garchNames <-
function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print slot names, @fit slot, parameters and controls

    # Arguments:
    #   object - an object of class 'fGARCH'

    # FUNCTION:

    # Slot Names:
    cat("\nNames - @ \n")
    print(slotNames(object))

    # @fit Slot:
    cat("\nNames - @fit \n")
    print(names(object@fit))

    # Parameters:
    cat("\nNames - @fit$params \n")
    print(names(object@fit$params))

    # Control:
    cat("\nNames - @fit$params$control \n")
    print(names(object@fit$params$control))

    # Return Value:
    invisible()
}


################################################################################

