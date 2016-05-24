# grid.R: functions for creating the grid of values to be plotted in plotmo graphs

# Get the x matrix (actually a data.frame) with median values (or
# first level for factors), ngrid1 rows, all rows identical.

get.degree1.xgrid <- function(x, grid.func, grid.levels, pred.names, ngrid1)
{
    stopifnot(!is.null(pred.names))
    check.grid.levels.arg(x, grid.levels, pred.names)
    xgrid <- data.frame(matrix(0, ngrid1, ncol(x), byrow=TRUE))
    for(ipred in seq_len(ncol(x)))
        xgrid[[ipred]] <- get.fixed.grid.lev(x[[ipred]], pred.names[ipred],
                                             grid.func, grid.levels)
    warn.if.not.all.finite(xgrid, "'xgrid' for degree1 plots")
    colnames(xgrid) <- pred.names
    xgrid
}
# Update xgrid for the predictor currently being plotted.
# That is, replace this predictor's column with a range of values.
# For factors or discrete variables, we shorten the frame to match the nbr of levels.

get.degree1.xframe <- function(xgrid, x, ipred, ngrid1,
                               ndiscrete, ux.list, extend)
{
    x1 <- x[[ipred]]
    u1 <- ux.list[[ipred]]
    if(is.factor(x1) && length(u1) > ngrid1)
        stop0("ngrid1=", ngrid1, " is less than the number ",
              "of levels in", colnames(x)[ipred])
    if(is.factor(x1) || length(u1) <= ndiscrete) {
        levels <- get.all.levs(x1, u1)
        xframe <- xgrid[1:length(levels), , drop=FALSE] # shorten xframe
        xframe[[ipred]] <- levels
    } else {
        xframe <- xgrid
        xrange <- range1(x1)
        if(extend != 0) {
            # extend the range of x (TODO consider allowing extend with discrete vars)
            stopifnot(xrange[2] >= xrange[1])
            ext <- extend * (xrange[2] - xrange[1])
            xrange[1] <- xrange[1] - ext
            xrange[2] <- xrange[2] + ext
        }
        xframe[[ipred]] <- seq(from=xrange[1], to=xrange[2], length=ngrid1)
    }
    xframe
}
# We want to display discrete variables in degree1 plots as quantized.
# (Factors get handled elsewhere).  So if a variable is discrete, then
# modify the xframe and yhat to do so.  For example, an xframe that was
#
#      pclass yhat
#           1  1.1
#           2  2.2
#           3  3.3
#
# becomes
#
#      pclass yhat
#         1.0  1.1
#         1.5  1.1
#         1.5  2.2
#         2.5  2.2
#         2.5  3.3
#         3.0  3.3

blockify.degree1.frame <- function(xframe, yhat, intervals,
                                   ipred, ux.list, ndiscrete)
{
    u1 <- ux.list[[ipred]]
    # TODO the integral check is necessary for compatibility with blockify.degree2.frame
    # (the code here can handle non integers but the code in blockify.degree2.frame can't)
    if(!is.factor(xframe[[ipred]]) && length(u1) <= ndiscrete && is.integral(u1)) {
        # discrete, so duplicate each elem in yhat
        yhat <- rep(yhat, each=2)
        if(!is.null(intervals)) {
            new.intervals <- data.frame(
                                lwr = rep(intervals$lwr, each=2),
                                upr = rep(intervals$upr, each=2))
            if(!is.null(intervals$fit))
                new.intervals$fit <- rep(intervals$fit, each=2)
            if(!is.null(intervals$cint.lwr)) {
                new.intervals$cint.lwr <- rep(intervals$cint.lwr, each=2)
                new.intervals$cint.upr <- rep(intervals$cint.upr, each=2)
            }
            intervals <- new.intervals
        }
        # duplicate each row of xframe, except the first and last row
        xframe <- xframe[rep(seq_len(nrow(xframe)), each=2), , drop=FALSE]
        if(nrow(xframe) >= 4) {
            x1 <- xframe[[ipred]]
            for(i in seq(2, length(x1)-1, by=2))
                x1[i] <- x1[i+1] <- (x1[i] + x1[i+1]) / 2
            xframe[[ipred]] <- x1
        }
    }
    list(xframe=xframe, yhat=yhat, intervals=intervals)
}
# Get the x matrix (actually a data.frame) to plot in degree2 plots.
# Each row of xgrid is identical (the medians).

get.degree2.xgrid <- function(x, grid.func, grid.levels, pred.names, ngrid2)
{
    check.grid.levels.arg(x, grid.levels, pred.names)
    xgrid <- list(ncol(x))
    for(ipred in seq_len(ncol(x)))
        xgrid[[ipred]] <- get.fixed.grid.lev(x[[ipred]], pred.names[ipred],
                                             grid.func, grid.levels)
    warn.if.not.all.finite(xgrid, "'xgrid' for degree2 plots")
    names(xgrid) <- pred.names
    xgrid <- as.data.frame(xgrid)
    xgrid[seq_len(ngrid2^2), ] <- xgrid
    xgrid
}
# Update xgrid for the predictor pair currently being plotted (ipred1
# and ipred2 are column numbers in x).  That is, replace two columns
# with a range of values.

get.degree2.xframe <- function(xgrid, x, ipred1, ipred2,
                               ngrid2, xranges, ux.list, ndiscrete)
{
    check.faclen <- function(ipred, nlevs)
    {
        if(is.factor(x[[ipred]]) && nlevs > ngrid2)
            stop0("ngrid2=", ngrid2, " is less than the number",
                  " of levels in ", colnames(x)[ipred])
    }
    n1 <- ngrid2 # will change only if ipred1 is treated as discrete
    u1 <- ux.list[[ipred1]]
    nlevs1 <- length(u1)
    check.faclen(ipred1, nlevs1)
    # TODO I'm not convinced the following code is correct (Apr 2015)
    if(is.factor(x[[ipred1]]) || nlevs1 <= ndiscrete) { # discrete?
        grid1 <- get.all.levs(x[[ipred1]], u1)
        n1 <- nlevs1
        xgrid <- xgrid[1:(n1 * ngrid2), , drop=FALSE]   # shorten xgrid
    } else
        grid1 <- seq(from=xranges[1,ipred1], to=xranges[2,ipred1], length=ngrid2)
    xgrid[[ipred1]] <- grid1 # will recycle

    n2 <- ngrid2 # will change only if ipred2 is treated as discrete
    u2 <- ux.list[[ipred2]]
    nlevs2 <- length(u2)
    check.faclen(ipred2, nlevs2)
    if(is.factor(x[[ipred2]]) || nlevs2 <= ndiscrete) { # discrete?
        grid2 <- get.all.levs(x[[ipred2]], u2)
        n2 <- nlevs2
        xgrid <- xgrid[rep(seq_len(n1), times=n2), , drop=FALSE]
    } else
        grid2 <- seq(from=xranges[1,ipred2], to=xranges[2,ipred2], length=ngrid2)
    xgrid[[ipred2]] <- rep(grid2, each=n1)

    list(xframe=xgrid, grid1=as.numeric(grid1), grid2=as.numeric(grid2))
}
# we want to draw discrete variables in persp and contour plots using "blocks"

blockify.degree2.frame <- function(x, yhat, grid1, grid2,
                                   ipred1, ipred2, ux.list, ndiscrete)
{
    is.discrete2 <- function(ipred, grid1)
    {
        if(is.factor(x[[ipred]]))
            return(TRUE)
        u1 <- ux.list[[ipred]]
        # the integral check is necessary with the current
        # implementation which adds/subtracts a hardcoded .499
        # TODO make this like blockify.degree1.frame (which can handle non integers)
        length(u1) <= ndiscrete && is.integral(grid1)
    }
    if(is.discrete2(ipred1, grid1)) {
        yhat <- rep(yhat, each=2)                   # duplicate each elem in yhat
        grid1 <- rep(grid1, each=2)                 # duplicate each elem in grid1
        is.even <- (1:length(grid1)) %% 2 == 0
        grid1[!is.even] <- grid1[!is.even] - .499   # sub .5 from odd elems
        grid1[is.even]  <- grid1[is.even]  + .499   # add .5 to even elems
    }
    if(is.discrete2(ipred2, grid2)) {
        # duplicate each block in yhat (each block has n1 elements)
        y.old <- yhat
        yhat <- double(2 * length(yhat))
        n1 <- length(grid1)
        for(i in 1:length(grid2)) {
            start <- n1 * (i-1)
            end   <- n1 * i
            yhat[(2 * start + 1): (2 * end)] <- y.old[(start + 1): end]
        }
        grid2 <- rep(grid2, each=2)         # duplicate each elem in grid2
        is.even <- (1:length(grid2)) %% 2 == 0
        grid2[!is.even] <- grid2[!is.even] - .499   # sub .5 from odd elems
        grid2[is.even]  <- grid2[is.even]  + .499   # add .5 to even elems
    }
    list(yhat=yhat, grid1=grid1, grid2=grid2)
}
# Check grid.levels arg passed in by the user.  This checks that the names
# of the list elements are indeed predictor names.  The actual levels will
# be checked later in get.fixed.grid.lev.

check.grid.levels.arg <- function(x, grid.levels, pred.names)
{
    if(!is.null(grid.levels)) { # null is the default value
        if(!is.list(grid.levels))
            stop0("grid.levels must be a list.  ",
                  "Example: grid.levels=list(sex=\"male\")")
        for(name in names(grid.levels))
            if(!pmatch(name, pred.names, 0))
                stop0("illegal variable name '", name, "' in grid.levels")
    }
}
# x1 is a column in the data matrix x
get.fixed.grid.lev <- function(x1, pred.name, grid.func, grid.levels)
{
    lev <- NA
    ilev <- 1
    if(!is.null(grid.levels)) {
        # look for pred.name in the grid.levels list, if found use its value
        iname <- which(pmatch(names(grid.levels), pred.name, duplicates.ok=TRUE) == 1)
        if(length(iname) > 1)
            stop0("illegal grid.levels argument (\"",
                  names(grid.levels)[iname[1]], "\" and \"",
                  names(grid.levels)[iname[2]],
                  "\" both match \"", pred.name, "\")")
        if(length(iname)) {
            lev <- grid.levels[[iname]]
            if(length(lev) > 1 || is.na(lev))
                stop0("illegal value for ", pred.name, " in grid.levels")
            if(is.factor(x1)) {
                lev.name <- grid.levels[[iname]]
                if(!is.character(lev.name)    ||
                        length(lev.name) != 1 ||
                        !nzchar(lev.name))
                    stop0("illegal level for \"",
                          names(grid.levels)[iname], "\" in grid.levels ",
                          "(specify factor levels with a string)")
                lev.names <- levels(x1)
                ilev <- pmatch(lev.name, lev.names, 0)
                if(!ilev)
                    stop0("illegal level \"", lev.name, "\" for \"",
                          pred.name, "\" in grid.levels (allowed levels are ",
                          quotify(lev.names), ")")
            }
            else if(!(is.numeric(lev) || is.logical(lev)) || !is.finite(lev))
                stop0("illegal value for ", pred.name, " in grid.levels")
        }
    }
    if(is.factor(x1)) {         # use ilev
        lev.names <- levels(x1)
        lev <- if(is.ordered(x1))
                   ordered(lev.names, levels=lev.names)[ilev]
               else
                   factor(lev.names, levels=lev.names)[ilev]
    } else if(is.na(lev)) {     # use lev
        if(!is.function(grid.func))
            stop0("'grid.func' is not a function");
        formals <- names(formals(grid.func))
        # check grid.func signature, we allow argname "na.rm" for mean and median
        if(length(formals) < 2 || formals[1] != "x" ||
                (!any(formals == "na.rm") && formals[2] != "..."))
            stop0("the formal arguments of 'grid.func' should be 'x' and '...'\n",
                  "       Your 'grid.func' has ",
                  if(length(formals) == 0)      "no formal arguments"
                  else if(length(formals) == 1) "a single formal argument "
                  else                          "formal arguments ",
                  if(length(formals) > 0) paste0("'", formals, "'", collapse=" ")
                  else "")
        lev <- grid.func(x1, na.rm=TRUE)
        if(length(lev) != 1 || (!is.numeric(lev) && !is.numeric(lev)))
            stop0("'grid.func' should return a single number\n",
                  "       Your 'grid.func' returned ", as.char(lev))
    }
    lev
}
# if x is a factor
#   return a factor vector with nlevs elements, e.g. pclass1, pclass2, pclass3.
# else
#   return a vector with all unique values in x, e.g. 1,2,3 or FALSE, TRUE

get.all.levs <- function(x, levels)
{
    if(!is.factor(x))
        return(levels)

    # TODO Sanity check, quite expensive, make sure no gaps in factor coding
    #      Could remove this if convert levels to factors in a better way below?
    range <- range(as.numeric(x))
    if(range[1] < 1 || range[2] > length(levels))
        stop0("internal error: illegal factor range ", range[1], " ",  range[2],
              " for levels ", quotify(levels))

    if(is.ordered(x))
        ordered(1:length(levels), labels=levels)
    else
        factor(1:length(levels), labels=levels)
}
# Print the grid values, must do some finagling for a nice display

print.grid.values <- function(xgrid, trace)
{
    trace1(trace, "\n")
    row <- xgrid[1, , drop=FALSE]
    names(row) <- c(paste("grid:   ", names(row)[1]), names(row)[-1])
    rownames(row) <- ""
    print(row)
    trace1(trace, "\n")
}
