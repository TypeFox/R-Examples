`pairs.profileModel` <-
function (x, colours = 2:3, title = NULL, ...) 
{
    ## 'pairs.profileModel' is a minor modification of 'pairs.profile' in the
     # MASS lirary.
     #'pairs.profile' was modified by Ioannis Kosmidis under GPL 2 or greater
     # and after the permission of the authors, in order to comply with objects
     # of class  "profileModel".
     # Ioannis Kosmidis <I.Kosmidis@warwick.ac.uk> [15/02/2008]
     #
    ## Another plot method for profile objects showing pairwise traces.
     # Recommended only for diagnostic purposes
    ######### Begin
    ## Addition by Ioannis Kosmidis [15/02/2008]
    if (!attr(x, "includes.traces")) {
        cat("Updating to get profile traces...\n")
        x <- update(x, verbose = FALSE, profTraces = TRUE)
    }
    ######### End
    ######### Begin
    ## Modification by Ioannis Kosmidis [15/02/2008]
    isnotNA <- !x$isNA
    profNames <- names(x$profiles[isnotNA])
    parvals <- lapply(x$profiles[isnotNA], FUN = function(obj) obj[, -2])
    parvals <- lapply(parvals, FUN = function(obj) obj[, profNames])
    rng <- apply(do.call("rbind", parvals), 2, range, na.rm = TRUE)
    Pnames <- colnames(rng)
    npar <- length(Pnames)
    coefs <- coef(x$fit)[isnotNA]
    form <- paste(as.character(formula(x$fit))[c(2, 1, 3)], collapse = "")
    ######### End
    oldpar <- par(mar = c(0, 0, 0, 0), mfrow = c(1, 1), oma = c(3, 
        3, 6, 3), las = 1)
    on.exit(par(oldpar))
    ##
    ## The following dodge ensures that the plot region is square
    ##
    fin <- par("fin")
    dif <- (fin[2] - fin[1])/2
    if (dif > 0) 
        adj <- c(dif, 0, dif, 0)
    else adj <- c(0, -dif, 0, -dif)
    par(omi = par("omi") + adj)
    ##
    ##
    cex <- 1 + 1/npar
    frame()
    mtext(form, side = 3, line = 3, cex = 1.5, outer = TRUE)
    del <- 1/npar
    for (i in 1:npar) {
        ci <- npar - i
        pi <- Pnames[i]
        for (j in 1:npar) {
            pj <- Pnames[j]
            par(fig = del * c(j - 1, j, ci, ci + 1))
            if (i == j) {
                par(new = TRUE)
                plot(rng[, pj], rng[, pi], axes = FALSE, xlab = "", 
                  ylab = "", type = "n")
                op <- par(usr = c(-1, 1, -1, 1))
                text(0, 0, pi, cex = cex, adj = 0.5)
                par(op)
            }
            else {
                col <- colours
                if (i < j) 
                  col <- col[2:1]
                if (!is.null(parvals[[pj]])) {
                  par(new = TRUE)
                  plot(spline(x <- parvals[[pj]][, pj], y <- parvals[[pj]][, 
                    pi]), type = "l", xlim = rng[, pj], ylim = rng[, 
                    pi], axes = FALSE, xlab = "", ylab = "", 
                    col = col[2])
                  pu <- par("usr")
                  smidge <- 2/100 * (pu[4] - pu[3])
                  segments(x, pmax(pu[3], y - smidge), x, pmin(pu[4], 
                    y + smidge))
                }
                else plot(rng[, pj], rng[, pi], axes = FALSE, 
                  xlab = "", ylab = "", type = "n")
                if (!is.null(parvals[[pi]])) {
                  lines(x <- parvals[[pi]][, pj], y <- parvals[[pi]][, 
                    pi], type = "l", col = col[1])
                  pu <- par("usr")
                  smidge <- 2/100 * (pu[2] - pu[1])
                  segments(pmax(pu[1], x - smidge), y, pmin(pu[2], 
                    x + smidge), y)
                }
                points(coefs[pj], coefs[pi], pch = 3, cex = 3)
            }
            if (i == npar) 
                axis(1)
            if (j == 1) 
                axis(2)
            if (i == 1) 
                axis(3)
            if (j == npar) 
                axis(4)
        }
    }
    par(fig = c(0, 1, 0, 1))
    if (!is.null(title)) {
        par(oma = c(0, 0, 2, 0))
        title(title, outer = TRUE)
    }
    invisible(x)
}
`plot.profileModel` <-
function (x, cis = NULL, signed = FALSE, interpolate = TRUE, 
    n.interpolations = 100, print.grid.points = FALSE, title = NULL, 
    ...) 
{
    fitted <- x$fit
    if (!is.null(cis)) {
        fitted.name <- x$call[["fitted"]]
        prof.name <- match.call()[["x"]]
        fitted.attr <- attr(cis, "fitted object")
        prof.attr <- attr(cis, "profileModel object")
        if (is.null(fitted.attr)) 
            fitted.attr <- 1
        if (is.null(prof.attr)) 
            prof.attr <- 1
        if (fitted.name == fitted.attr | prof.name == prof.attr) {
        }
        else stop("Invalid confidence intervals were supplied.")
    }
    if (!(agreement <- x$agreement) & signed) 
        stop("The objective and the fitting procedure ", fitted$call[[1]], 
            " do not agree. Signed square roots cannot be calculated.")
    op <- par(no.readonly = TRUE)
    if (is.null(x$quantile)) {
        if (signed) {
            x <- signedSquareRoots.profileModel(x)
            temp.plot <- function(mat, nam) {
                plot(mat[, 1], mat[, 2], type = "l", xlab = nam, 
                  ylab = "Signed sqrt of objective")
            }
        }
        else {
            temp.plot <- function(mat, nam) {
                plot(mat[, 1], mat[, 2], type = "l", xlab = nam, 
                  ylab = "Profiled objective")
            }
        }
    }
    else {
        if (signed) {
            x <- signedSquareRoots.profileModel(x)
            temp.plot <- function(mat, nam) {
                plot(mat[, 1], mat[, 2], type = "l", xlab = nam, 
                  ylab = "Signed sqrt of objective")
                points(x = c(min(mat[, 1]), max(mat[, 1])), y = rep(-sqrt(x$quantile), 
                  2), type = "l", lty = 2)
                points(x = c(min(mat[, 1]), max(mat[, 1])), y = rep(sqrt(x$quantile), 
                  2), type = "l", lty = 2)
            }
        }
        else temp.plot <- function(mat, nam) {
            plot(mat[, 1], mat[, 2], type = "l", xlab = nam, 
                ylab = "Profiled objective")
            points(x = c(min(mat[, 1]), max(mat[, 1])), y = rep(x$quantile, 
                2), type = "l", lty = 2)
        }
    }
    profRes.or <- profRes <- x$profiles
    isNA <- x$isNA
    p <- length(profRes)
    profNames <- names(profRes)
    which <- x$profiled.parameters
    scale <- !is.null(Xmax <- fitted$X.max.scaleFit)
    Betas <- coef(fitted)[which]/(if (scale) 
        Xmax[which]
    else 1)
    if (agreement) {
        res.at.betas <- as.list(rep(NA, p))
        names(res.at.betas) <- profNames
        for (i in 1:p) {
            if (isNA[i]) 
                next
            res.at.betas[[profNames[i]]] <- matrix(c(Betas[i], 
                0), 1, 2)
        }
    }
    else suppressWarnings(res.at.betas <- update(x, grid.bounds = cbind(Betas, 
        Betas), gridsize = 1, quantile = NULL, verbose = FALSE, 
        profTraces = FALSE)$profiles)
    if (interpolate) {
        for (i in 1:p) {
            if (isNA[i]) 
                next
            profRes.i <- profRes[[i]][, 1:2]
            ### construct some information for the spline to use
            lin <- approx(profRes.i, n = 2 * nrow(profRes.i))
            smoothed <- spline(lin, n = n.interpolations)
            profRes[[i]] <- cbind(smoothed$x, smoothed$y)
        }
    }
    intersects <- x$intersects
    has.prelim <- attr(x$grid.bounds, "from.prelim")
    par(mfrow = c(ceiling(sqrt(p)), ceiling(sqrt(p))))
    for (i in 1:p) {
        if (isNA[i]) 
            next
        profRes.i <- profRes[[i]]
        profNames.i <- profNames[i]
        temp.plot(profRes.i, profNames.i)
        min.i <- min(profRes.i[, 2])
        max.i <- max(profRes.i[, 2])
        if (has.prelim) {
            intersects.i <- intersects[i, ]
            # draw cis
            if (!is.null(cis)) {
                cis.i <- cis[i, ]
                if (all(intersects.i)) {
                  points(x = rep(cis.i[1], 2), y = c(min.i, max.i), 
                    type = "l", lty = 3)
                  points(x = cis.i[1], y = min.i, pch = 6)
                  points(x = rep(cis.i[2], 2), y = c(min.i, max.i), 
                    type = "l", lty = 3)
                  points(x = cis.i[2], y = min.i, pch = 6)
                }
                if (sum(intersects.i) == 1) {
                  which.intersects.i <- which(intersects.i)
                  points(x = rep(cis.i[which.intersects.i], 2), 
                    y = c(min.i, max.i), type = "l", lty = 3)
                  points(x = cis.i[which.intersects.i], y = min.i, 
                    pch = 6)
                }
            }
            # draw estimates
            if (all(intersects.i) | all(!intersects.i)) {
                points(res.at.betas[[i]], pch = 4)
            }
            else {
                which.intersects.i <- which(intersects.i)
                if (which.intersects.i == 1) {
                  if (agreement) 
                    text(x = max(profRes.i[, 1]), y = 0, labels = expression(infinity))
                  else text(x = max(profRes.i[, 1]), y = min.i, 
                    labels = expression(infinity))
                }
                if (which.intersects.i == 2) {
                  if (agreement) 
                    text(x = min(profRes.i[, 1]), y = 0, labels = expression(-infinity))
                  else text(x = min(profRes.i[, 1]), y = min.i, 
                    labels = expression(-infinity))
                }
                if (!agreement) 
                  points(res.at.betas[[i]], pch = 4)
            }
        }
        if (print.grid.points) 
            points(x = profRes.or[[i]][, 1], y = profRes.or[[i]][, 
                2], pch = 16, cex = 0.6)
        title(profNames.i)
    }
    ######### Begin
    ## Addition by Ioannis Kosmidis [15/02/2008]
    if (!is.null(title)) {
        par(oma = c(0, 0, 2, 0))
        title(title, outer = TRUE)
    }
    ######### End
    par(op)
}
