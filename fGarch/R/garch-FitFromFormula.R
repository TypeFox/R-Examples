
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
# FUNCTION:               PARAMETER ESTIMATION:
#  garchFit                Fits the parameters of GARCH process
#  .garchArgsParser        Parses formula and data for garchFit
#  .garchOptimizerControl  Sets default values for Garch Optimizer
#  .garchFit               ... old Version, still in use by garchFit()
#  .garchNames             Slot names, @fit slot, parameters and controls
################################################################################


garchFit <-
function(
    formula = ~ garch(1,1), 
    data = dem2gbp,
    init.rec = c("mci", "uev"),
    delta = 2, 
    skew = 1, 
    shape = 4,
    cond.dist = c("norm", "snorm", "ged", "sged", "std", "sstd", 
        "snig", "QMLE"),
    include.mean = TRUE, 
    include.delta = NULL, 
    include.skew = NULL,
    include.shape = NULL, 
    leverage = NULL,
    trace = TRUE,
    algorithm = c("nlminb", "lbfgsb", "nlminb+nm", "lbfgsb+nm"),
    hessian = c("ropt", "rcd"), 
    control = list(),
    title = NULL,
    description = NULL, 
    ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fit parameters to a ARMA-GARCH model by Formula Specification

    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance
    #       specification
    #   data - any univariate time series which can be converted
    #       into a timeSeries using the generic function as.timeSeries
    #   init.rec - names type of initialization of recurrence
    #       mci = mu-current-iteration, or
    #       uev = unconditional-expected-variances
    #   delta - numeric value of the exponent delta
    #   skew - optional skewness or skewness related parameter
    #   shape - optional shape parameter
    #   cond.dist - name of the conditional distribution, one of
    #       norm, snorm, ged, sged, std, sstd, snig, QMLE   
    #   include.mean - a logical, should the mean value be estimated ?
    #   include.delta - should the exponent be estimated ?
    #   leverage - should the leverage factors be estimated ?
    #   trace - should the optimization be traced ?
    #   control - list of additional control parameters for the solver
    #   title - an optional title string
    #   description - an optional project description string
    
    # Example:
    #   garchFit() 

    # FUNCTION:
    
    # DEBUG:
    DEBUG = FALSE

    # Match arguments:
    init.rec = match.arg(init.rec)
    cond.dist = match.arg(cond.dist)
    hessian = match.arg(hessian)
    algorithm = match.arg(algorithm)

    # Call:
    CALL = match.call()

    # Parse formula and data for garchFit ...
    #   Note in the new version we are working with timeSeries ...
    Name = capture.output(substitute(data))
    if(is.character(data)) {
        eval(parse(text = paste("data(", data, ")")))
        data = eval(parse(text = data))
    }
    # data <- if (inherits(data, "timeSeries") data else as.timeSeries(data)
    data <- as.data.frame(data)

    # Column Names:
    if (isUnivariate(data)) {
        colnames(data) <- "data"
    } else {
        # Check unique column Names:
        uniqueNames = unique(sort(colnames(data)))
        if (is.null(colnames(data))) {
            stop("Column names of data are missing.")
        }
        if (length(colnames(data)) != length(uniqueNames)) {
            stop("Column names of data are not unique.")
        }
    }

    # Handle if we have no left-hand-side for the formula ...
    #   Note in this case the length of the formula is 2 (else 3):
    if (length(formula) == 3 && isUnivariate(data) ) formula[2] <- NULL
    if (length(formula) == 2) {
        if (isUnivariate(data)) {
            # Missing lhs -- we substitute the data file name as lhs ...
            formula = as.formula(paste("data", paste(formula, collapse = " ")))
        } else {
            stop("Multivariate data inputs require lhs for the formula.")
        }
    }

    # Robust Covariance ?
    robust.cvar <- (cond.dist == "QMLE")

    # Parse Arguments:
    args = .garchArgsParser(formula = formula, data = data, trace = FALSE)

    # DEBUG - Print Arguments:
    if (DEBUG) print(list(
        formula.mean = args$formula.mean,
        formula.var = args$formula.var,
        series = args$series,
        init.rec = init.rec, 
        delta = delta, 
        skew = skew, 
        shape = shape, 
        cond.dist = cond.dist, 
        include.mean = include.mean,
        include.delta = include.delta, 
        include.skew = include.skew, 
        include.shape = include.shape, 
        leverage = leverage,
        trace = trace,
        algorithm = algorithm,
        hessian = hessian,
        robust.cvar = robust.cvar,
        control = control,
        title = title, 
        description = description))
        
    # Fit:
    ans = .garchFit(
        formula.mean = args$formula.mean,
        formula.var = args$formula.var,
        series = args$series,
        init.rec, 
        delta, 
        skew, 
        shape, 
        cond.dist, 
        include.mean,
        include.delta, 
        include.skew, 
        include.shape, 
        leverage,
        trace,
        algorithm,
        hessian,
        robust.cvar,
        control,
        title, 
        description, 
        ...)
    ans@call = CALL
    attr(formula, "data") <- paste("data = ", Name, sep = "")
    ans@formula = formula

    # Return Value:
    ans
}


################################################################################

