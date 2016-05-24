################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Space-time K-function analysis of "epidataCS" objects
### along the lines of Diggle et al (1995):
### "Second-order analysis of space-time clustering" (Stat Methods Med Res)
###
### Copyright (C) 2015 Sebastian Meyer
### $Revision: 1347 $
### $Date: 2015-05-29 11:45:51 +0200 (Fre, 29. Mai 2015) $
################################################################################

## call K-function methods in package "splancs"
stKcall <- function (which = c("stkhat", "stsecal", "stmctest"),
                     object, eps.s, eps.t, ...)
{
    stopifnot(inherits(object, "epidataCS"))
    
    ## get the function
    which <- match.arg(which)
    FUN <- get(which, mode = "function", envir = getNamespace("splancs"))
    
    ## default arguments
    commonArgs <- list(
        pts = coordinates(object$events), times = object$events$time,
        poly = NULL,                      tlimits = summary(object)$timeRange,
        s = eps.s,                        tm = eps.t
    )
    args <- modifyList(commonArgs, list(...))
    if (is.null(args$poly)) { # use coordinates of first polygon
        if (length(object$W) > 1L || length(object$W@polygons[[1]]@Polygons) > 1L)
            stop("package \"splancs\" does not support multi-'poly'gons")
        args$poly <- coordinates(object$W@polygons[[1L]]@Polygons[[1L]])
    }
    if (which == "stmctest" && is.null(args[["nsim"]])) {
        args$nsim <- 199L
    }
    
    ## unfortunately, argument names are not consistent across functions
    if (which == "stsecal")
        names(args)[names(args) == "tlimits"] <- "tlim"
    if (which == "stmctest")
        names(args)[names(args) == "tm"] <- "tt"

    ## call the selected splancs function
    do.call(FUN, args)
}

## Monte-Carlo test for space-time interaction 
stKtest <- function (object, eps.s = NULL, eps.t = NULL, B = 199,
                     cores = 1, seed = NULL, poly = object$W)
{
    stopifnot(inherits(object, "epidataCS"),
              isScalar(cores), cores > 0, isScalar(B), B > 0)
    cores <- as.integer(cores)
    B <- as.integer(B)
    
    ## naive default grids
    if (is.null(eps.s))
        eps.s <- seq(0, min(object$events$eps.s, apply(bbox(object$W), 1, diff)/2),
                     length.out = 10)
    if (is.null(eps.t))
        eps.t <- seq(0, min(object$events$eps.t, tail(object$stgrid$stop,1L)/2),
                     length.out = 10)
    
    ## extract coordinates of the polygon
    polycoordslist <- xylist(poly)
    if (length(polycoordslist) > 1L) {
        stop("package \"splancs\" does not support multi-'poly'gons")
    }
    Wcoords <- as.matrix(as.data.frame(polycoordslist[[1L]]))

    ## calculate K-function
    stK <- stKcall("stkhat", object = object, eps.s = eps.s, eps.t = eps.t,
                   poly = Wcoords)

    ## calculate standard error
    seD <- stKcall("stsecal", object = object, eps.s = eps.s, eps.t = eps.t,
                   poly = Wcoords)

    ## perform Monte Carlo permutation test (parallelized)
    permt <- plapply(
        X = diff(round(seq(from = 0, to = B, length.out = cores + 1L))),
        FUN = function (nsim) {
            stKcall("stmctest", object = object, eps.s = eps.s, eps.t = eps.t,
                    poly = Wcoords, nsim = nsim, quiet = TRUE)[["t"]]
        },
        .parallel = cores, .seed = seed, .verbose = FALSE
    )
    mctest <- list(
        "t0" = sum(stK$kst - outer(stK$ks, stK$kt)),
        "t" = unlist(permt, recursive = FALSE, use.names = FALSE)
    )
    PVAL <- mean(c(mctest[["t0"]], mctest[["t"]]) >= mctest[["t0"]])

    ## return test results
    structure(
        list(method = "Diggle et al (1995) K-function test for space-time clustering",
             data.name = deparse(substitute(object)),
             statistic = setNames(mctest$t0, "U"), # sum of residuals
             parameter = setNames(B, "B"), p.value = PVAL,
             pts = coordinates(object$events),
             stK = stK, seD = seD, mctest = mctest),
        class = c("stKtest", "htest")
    )
}

## diagnostic plots related to space-time K-function analysis
## inspired by splancs::stdiagn authored by Barry Rowlingson and Peter Diggle
plot.stKtest <- function (x, which = c("D", "R", "MC"),
                          args.D = list(), args.D0 = args.D,
                          args.R = list(), args.MC = list(),
                          mfrow = sort(n2mfrow(length(which))), ...)
{
    stkh <- x$stK
    stse <- x$seD
    stmc <- x$mctest
    
    if (identical(which, "stdiagn")) {
        splancs::stdiagn(pts = x$pts, stkh = stkh, stse = stse, stmc = stmc)
        return(invisible())
    }
    
    which <- match.arg(which, several.ok = TRUE)
    stopifnot(is.list(args.D), is.list(args.D0), is.list(args.R), is.list(args.MC))

    ## K_0(s,t) = K(s) * K(t)
    K0 <- outer(stkh$ks, stkh$kt)
    ## D(s,t) = K(s,t) - K_0(s,t)
    st.D <- stkh$kst - K0

    if (!is.null(mfrow)) {
        omfrow <- par(mfrow = mfrow)
        on.exit(par(omfrow))
    }

    ## D plots
    Dzero <- which[which %in% c("D", "D0")] == "D0"
    whichDzero <- match(Dzero, c(FALSE, TRUE))
    omar <- par(mar = if (is.null(args.D[["mar"]])) c(2,2,par("mar")[3L],1) else args.D[["mar"]])
    mapply(
        FUN = function (z, Dzero, args) {
            defaultArgs <- list(
                x = stkh$s, y = stkh$t, z = z,
                main = if (Dzero) "Excess risk" else "D plot",
                xlab = "Distance", ylab = "Time lag", zlab = "",
                ticktype = "detailed", shade = 0.5, col = "lavender",
                theta = -30, phi = 15, expand = 0.5
            )
            do.call("persp", modifyList(defaultArgs, args))
        },
        z = list(st.D, st.D/K0)[whichDzero],
        Dzero = Dzero, args = list(args.D, args.D0)[whichDzero],
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    par(omar)

    ## Residual plot
    if ("R" %in% which) {
        st.R <- st.D/stse
        defaultArgs.R <- list(
            x = K0, y = st.R,
            panel.first = quote(abline(h = c(-2,0,2), lty = c(2,1,2))),
            xlab = "K(s)K(t)", ylab = "R", main = "Standardized residuals",
            ylim = range(0, st.R, finite = TRUE)
        )
        do.call("plot.default", modifyList(defaultArgs.R, args.R))
    }

    ## MC permutation test plot
    if ("MC" %in% which) {
        defaultArgs.MC <- list(
            permstats = stmc$t,
            xmarks = setNames(stmc$t0, "observed"),
            main = "MC permutation test"
        )
        do.call("permtestplot", modifyList(defaultArgs.MC, args.MC))
    }

    invisible()
}
