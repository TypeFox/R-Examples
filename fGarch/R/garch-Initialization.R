
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
# FUNCTION:               DESCRIPTION:
#  .garchInitSeries        Initializes Series
#  .garchInitParameters    Initializes Parameters
################################################################################


.garchInitSeries <-
function(formula.mean, formula.var, cond.dist, series, scale, init.rec,
    h.start, llh.start, trace)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Initialize time series

    # Arguments:
    #   see function garchFit()

    # FUNCTION:

    # Check Mean Formula ARMA - Is it Valid ?
    mm = length(formula.mean)
    if(mm != 2) stop("Mean Formula misspecified")
    end = regexpr("\\(", as.character(formula.mean[mm])) - 1
    model.mean = substr(as.character(formula.mean[mm]), 1, end)
    if(!any( c("ar", "ma", "arma") == model.mean))
        stop("formula.mean must be one of: ar, ma, arma")

    # Check Variance Formula GARCH - Is it Valid ?
    mv = length(formula.var)
    if(mv != 2) stop("Variance Formula misspecified")
    end = regexpr("\\(", as.character(formula.var[mv])) - 1
    model.var = substr(as.character(formula.var[mv]), 1, end)
    if(!any( c("garch", "aparch") == model.var))
        stop("formula.var must be one of: garch, aparch")

    # Determine Mean Order from ARMA Formula:
    model.order = as.numeric(strsplit(strsplit(strsplit(as.character(
        formula.mean), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
    u = model.order[1]
    v = 0
    if(length(model.order) == 2) v = model.order[2]
    maxuv = max(u, v)
    if(u < 0 | v < 0) stop("*** ARMA orders must be positive.")

    # Determine Variance Order from GARCH Formula:
    model.order = as.numeric(strsplit(strsplit(strsplit(as.character(
        formula.var), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
    p = model.order[1]
    q = 0
    if(length(model.order) == 2) q = model.order[2]
    if(p+q == 0)
        stop("Misspecified GARCH Model: Both Orders are zero!")
    maxpq = max(p, q)
    if(p < 0 | q < 0) stop("*** GARCH orders must be positive.")

    # Fix Start Position of Series "h" and for Likelihood Calculation:
    max.order = max(maxuv, maxpq)
    if(is.null(h.start)) h.start = max.order + 1
    if(is.null(llh.start)) llh.start = 1

    # Check for Recursion Initialization:
    if(init.rec != "mci" & model.var != "garch") {
        stop("Algorithm only supported for mci Recursion")
    }

    # Trace the Result:
    if(trace) {
        cat("\nSeries Initialization:")
        cat("\n ARMA Model:               ", model.mean)
        cat("\n Formula Mean:             ", as.character(formula.mean))
        cat("\n GARCH Model:              ", model.var)
        cat("\n Formula Variance:         ", as.character(formula.var))
        cat("\n ARMA Order:               ", u, v)
        cat("\n Max ARMA Order:           ", maxuv)
        cat("\n GARCH Order:              ", p, q)
        cat("\n Max GARCH Order:          ", maxpq)
        cat("\n Maximum Order:            ", max.order)
        cat("\n Conditional Dist:         ", cond.dist)
        cat("\n h.start:                  ", h.start)
        cat("\n llh.start:                ", llh.start)
        cat("\n Length of Series:         ", length(series))
        cat("\n Recursion Init:           ", init.rec)
        cat("\n Series Scale:             ", scale)
        cat("\n\n")
    }

    # Result:
    ans  = list(
        model = c(model.mean, model.var),
        order = c(u = u, v = v, p = p, q = q),
        max.order = max.order,
        z = rep(0, times = length(series)),
        h = rep(var(series), times = length(series)),
        x = series,
        scale = scale,
        init.rec = init.rec,
        h.start = h.start,
        llh.start = llh.start)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.garchInitParameters <-
function(formula.mean, formula.var, delta, skew, shape, cond.dist,
    include.mean, include.delta, include.skew, include.shape, leverage,
    algorithm, control, trace)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Initialize model parameters

    # Arguments:
    #   see function garchFit()

    # FUNCTION:

    # DEBUG:
    .DEBUG = FALSE

    # global variables
    .series <- .getfGarchEnv(".series")

    # Determine Mean Order from ARMA Formula:
    model.order = as.numeric(strsplit(strsplit(strsplit(as.character(
        formula.mean), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
    u = model.order[1]
    v = 0
    if(length(model.order) == 2) v = model.order[2]

    # Determine Variance Order from GARCH Formula:
    model.order = as.numeric(strsplit(strsplit(strsplit(as.character(
        formula.var), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
    p = model.order[1]
    if (p == 0) stop("The order p must be > 0 in GARCH/APARCH(p,q)")
    q = 0
    if(length(model.order) == 2) q = model.order[2]

    # Includes:
    model.var = .series$model[2]
    if(is.null(include.delta)) {
        if(model.var == "garch") {
            include.delta = FALSE
        } else {
            include.delta = TRUE
        }
    }
    if(is.null(leverage)) {
        if(model.var == "garch") {
            leverage = FALSE
        } else {
            leverage = TRUE
        }
    }

    # Distributional Includes:
    if(cond.dist == "t") cond.dist = "std"
    skewed.dists = c("snorm", "sged", "sstd", "snig")
    if(is.null(include.skew)) {
        if(any(skewed.dists == cond.dist)) {
            include.skew = TRUE
        } else {
            include.skew = FALSE
        }
    }
    shaped.dists = c("ged", "sged", "std", "sstd", "snig")
    if(is.null(include.shape)) {
        if(any(shaped.dists == cond.dist)) {
            include.shape = TRUE
        } else {
            include.shape = FALSE
        }
    }

    # Set Names for Parameters:
    Names = c(
        "mu",
        if(u > 0) paste("ar", 1:u, sep = ""),
        if(v > 0) paste("ma", 1:v, sep = ""),
        "omega",
        if(p > 0) paste("alpha", 1:p, sep = ""),
        if(p > 0) paste("gamma", 1:p, sep = ""),
        if(q > 0) paste("beta",  1:q, sep = ""),
        "delta",
        "skew",
        "shape")
    if(.DEBUG) { cat("\nDEBUG - Names: \n"); print(Names) }

    # Initialize Model Parameters to be Estimated:
    fit.mean = arima(.series$x, order = c(u, 0, v),
        include.mean = include.mean)$coef
    alpha.start = 0.1
    beta.start = 0.8
    ## if(include.delta) delta = 1.5
    params = c(
        if(include.mean) fit.mean[length(fit.mean)] else 0,
        if(u > 0) fit.mean[1:u],
        if(v > 0) fit.mean[(u+1):(length(fit.mean)-as.integer(include.mean))],
        var(.series$x, na.rm = TRUE)*(1-alpha.start-beta.start),
        if(p > 0) rep(alpha.start/p, times = p),
        if(p > 0) rep(0.1, times = p),
        if(q > 0) rep(beta.start/q, times = q),
        delta,
        skew,
        shape)
    names(params) = Names
    if(.DEBUG) { cat("\nDEBUG - params: \n"); print(params) }

    # Set Lower Limits of Parameters to be Estimated:
    TINY = 1.0e-8
    USKEW = 1/10; USHAPE = 1
    if (cond.dist == "snig") USKEW = -0.99
    U = c(
        -10*abs(mean(.series$x)),
        if(u > 0) rep(-1+TINY, times = u),
        if(v > 0) rep(-1+TINY, times = v),
        1.0e-6*var(.series$x),
        if(p > 0) rep( 0+TINY, times = p),
        if(p > 0) rep(-1+TINY, times = p),
        if(q > 0) rep( 0+TINY, times = q),
        0,          # delta
        USKEW,      # skew
        USHAPE)     # shape
    names(U) = Names
    if(.DEBUG) { cat("\nDEBUG - U: \n"); print(U) }

    # Set Upper Limits of Parameters to be Estimated:
    VSKEW = 10; VSHAPE = 10
    if (cond.dist == "snig") VSKEW = 0.99
    V = c(
        10*abs(mean(.series$x)),
        if(u > 0) rep(1-TINY, times = u),
        if(v > 0) rep(1-TINY, times = v),
        100*var(.series$x),
        if(p > 0) rep(1-TINY, times = p),
        if(p > 0) rep(1-TINY, times = p),
        if(q > 0) rep(1-TINY, times = q),
        2,          # delta
        VSKEW,      # skew
        VSHAPE)     # shape
    names(V) = Names
    if(.DEBUG) { cat("\nDEBUG - V: \n"); print(V) }

    # Includes:
    includes = c(
        include.mean,
        if(u > 0) rep(TRUE, times = u),
        if(v > 0) rep(TRUE, times = v),
        TRUE,
        if(p > 0) rep(TRUE, times = p),
        if(p > 0) rep(leverage, times = p),
        if(q > 0) rep(TRUE, times = q),
        include.delta,
        include.skew,
        include.shape)
    names(includes) = Names
    if(.DEBUG) { cat("\nDEBUG - V: \n"); print(includes) }

    # Index List of Parameters to be Optimized:
    index = (1:length(params))[includes == TRUE]
    names(index) = names(params)[includes == TRUE]
    if(.DEBUG) { cat("\nDEBUG - fixed: \n"); print(index) }

    # Persistence:
    alpha <- beta <- NULL
    if(p > 0) alpha = params[substr(Names, 1, 5) == "alpha"]
    if(p > 0 & leverage) gamma = params[substr(Names, 1, 5) == "gamma"]
    if(p > 0 & !leverage) gamma = rep(0, times = p)
    if(q > 0) beta  = params[substr(Names, 1, 4) == "beta"]
    if(.series$model[2] == "garch") {
        persistence = sum(alpha) + sum(beta)
    } else if(.series$model[2] == "aparch") {
        persistence = sum(beta)
        for (i in 1:p)
            persistence = persistence + alpha[i]*garchKappa(cond.dist,
                gamma[i], params["delta"], params["skew"], params["shape"])
    }
    names(persistence) = "persistence"

    # Trace the Result:
    if(trace) {
        cat("Parameter Initialization:")
        cat("\n Initial Parameters:          $params")
        cat("\n Limits of Transformations:   $U, $V")
        cat("\n Which Parameters are Fixed?  $includes")
        cat("\n Parameter Matrix:\n")
        ans = data.frame(U, V, params, includes)
        rownames(ans) = paste("   ", names(params))
        print(ans)
        cat(" Index List of Parameters to be Optimized:\n")
        print(index)
        cat(" Persistence:                 ", persistence, "\n")
    }

    # Return Value:
    list(params = params, 
        U = U, 
        V = V, 
        includes = includes,
        index = index, 
        mu = params[1], 
        delta = delta, 
        skew = skew,
        shape = shape, 
        cond.dist = cond.dist, 
        leverage = leverage,
        persistence = persistence, 
        control = control)
}


################################################################################

