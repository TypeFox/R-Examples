#' Description: Plot Potential
#'
#' Visualization of the potential function from the movpotential function
#'
#' Arguments:
#'    @param res output from movpotential function
#'    @param title title text
#'    @param xlab.text xlab text
#'    @param ylab.text ylab text
#'    @param cutoff parameter determining the upper limit of potential for visualizations
#'    @param plot.contours Plot contour lines.
#'    @param binwidth binwidth for contour plot
#'    @param bins bins for contour plot. Overrides binwidth if given
#' 
#' @importFrom tgp interp.loess
#' @importFrom ggplot2 ggplot
#' @return \item{ggplot2}{potential plotted}
#'
#' @export
#'
#' @references Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples # X = c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), rnorm(1000, mean = 2)); param = seq(0,5,length=3000); res <- movpotential_ews(X, param); PlotPotential(res$res, title = '', xlab.text = '', ylab.text = '', cutoff = 0.5, plot.contours = TRUE, binwidth = 0.2)
#'
#' @keywords early-warning

PlotPotential <- function(res, title = "", xlab.text, ylab.text, cutoff = 0.5, plot.contours = TRUE, 
    binwidth = 0.2, bins = NULL) {
    
    cut.potential <- max(apply(res$pots, 1, min)) + cutoff * abs(max(apply(res$pots, 
        1, min)))  # Ensure all minima are visualized
    pots <- res$pots
    pots[pots > cut.potential] <- cut.potential
    
    # Static contour Interpolate potential grid
    intp <- tgp::interp.loess(as.vector(res$pars), as.vector(res$xis), as.vector(pots), 
        gridlen = 2 * dim(pots))
    xy <- expand.grid(intp$x, intp$y)
    z <- as.vector(intp$z)
    z[is.na(z)] <- max(na.omit(z))
    
    potential <- NULL
    df <- data.frame(list(bg.var = xy[, 1], phylotype = xy[, 2], potential = z))
    bg.var <- NULL
    phylotype <- NULL
    
    p <- ggplot2::ggplot(df, aes(bg.var, phylotype, z = potential)) + geom_tile(aes(fill = potential)) + 
        scale_fill_gradient(low = "black", high = "white")
    
    if (plot.contours) {
        if (!is.null(bins)) {
            warning("bins argument is overriding the binwidth argument!")
            p <- p + stat_contour(bins = bins)
        } else {
            p <- p + stat_contour(binwidth = binwidth)
        }
        # p <- p + stat_contour(geom='polygon', aes(fill=..level..))
    }
    
    p <- p + xlab(xlab.text) + ylab(ylab.text) + labs(title = title)
    
    p
    
}







#' Description: Potential Analysis for univariate data
#'
#' \code{livpotential_ews} performs one-dimensional potential estimation derived from a uni-variate timeseries
#'
#' Arguments:
#'    @param x Univariate data (vector) for which the potentials shall be estimated
#'    @param std Standard deviation of the noise (defaults to 1; this will set scaled potentials)
#'    @param bw bandwidth for kernel estimation
#'    @param weights optional weights in ksdensity (used by movpotentials).
#'    @param grid.size Grid size for potential estimation.
#'    @param detection.threshold maximum detection threshold as fraction of density kernel height dnorm(0, sd = bandwidth)/N
#'    @param bw.adjust The real bandwidth will be bw.adjust*bw; defaults to 1
#'    @param density.smoothing Add a small constant density across the whole observation range to regularize density estimation (and to avoid zero probabilities within the observation range). This parameter adds uniform density across the observation range, scaled by density.smoothing.
#'    @param detection.limit ignore maxima that are below detection.limit * maximum density
#' 
# Returns:
#'   @return \code{livpotential} returns a list with the following elements:
#'   @return \item{xi}{the grid of points on which the potential is estimated}
#'   @return \item{pot}{The estimated potential: -log(f)*std^2/2, where f is the density.}
#'   @return \item{density}{Density estimate corresponding to the potential.}
#'   @return \item{min.inds}{indices of the grid points at which the density has minimum values; (-potentials; neglecting local optima)}
#'   @return \item{max.inds}{indices the grid points at which the density has maximum values; (-potentials; neglecting local optima)}
#'   @return \item{bw}{bandwidth of kernel used}
#'   @return \item{min.points}{grid point values at which the density has minimum values; (-potentials; neglecting local optima)}
#'   @return \item{max.points}{grid point values at which the density has maximum values; (-potentials; neglecting local optima)}
#'
#' @export
#'
#' @references Livina, VN, F Kwasniok, and TM Lenton, 2010. Potential analysis reveals changing number of climate states during the last 60 kyr . \emph{Climate of the Past}, 6, 77-82.
#' 
#' Dakos, V., et al (2012).'Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data.' \emph{PLoS ONE} 7(7): e41010. doi:10.1371/journal.pone.0041010
#' @author Based on Matlab code from Egbert van Nes modified by Leo Lahti. Implemented in early warnings package by V. Dakos.
#' @seealso 
#' \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}}; \code{\link{ch_ews}};\code{\link{movpotential_ews}}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples 
#' data(foldbif)
#' res <- livpotential_ews(foldbif)
#' @keywords early-warning

livpotential_ews <- function(x, std = 1, bw = "nrd", weights = c(), grid.size = NULL, 
    detection.threshold = 0.01, bw.adjust = 1, density.smoothing = 0, detection.limit = 0.1) {
    
    # std <- 1; bw <- 'nrd'; weights <- c(); grid.size = floor(.2*ncol(mydat));
    # detection.threshold = det.th; density.smoothing = 0
    
    x <- data.frame(x)
    
    if (is.null(grid.size)) {
        grid.size <- floor(0.2 * nrow(x))
    }
    
    if (is.null(bw)) {
        # following Silverman, B. W.: Density estimation for statistics and data
        # analysis,Chapman and Hall, 1986.
        bw <- 1.06 * sapply(x, sd)/nrow(x)^(1/5)
    }
    
    # Density estimation
    de <- density(ts(x), bw = bw, adjust = bw.adjust, kernel = "gaussian", weights = weights, 
        window = kernel, n = grid.size, from = min(x), to = max(x), cut = 3, na.rm = FALSE)
    
    # Estimated density f <- de$y
    
    # Smooth the density by adding a small probability across the whole observation
    # range (to avoid zero probabilities for points in the observation range)
    f <- de$y + density.smoothing * 1/diff(range(de$x))  # *max(de$y)
    
    # Normalize the density such that it integrates to unity
    f <- f/sum(diff(de$x[1:2]) * f)
    
    # Final grid points and bandwidth
    grid.points <- de$x
    bw <- de$bw
    
    # Compute potential
    U <- -log(f) * std^2/2
    # f <- exp(-2*U/std^2) # backtransform to density distribution Ignore very local
    # optima Note mins and maxs for density given here (not for potential, which has
    # the opposite signs)
    
    ops <- find.optima(f, detection.threshold = detection.threshold, bw = bw, x = x, 
        detection.limit = detection.limit)
    min.points <- grid.points[ops$min]
    max.points <- grid.points[ops$max]
    det.th <- ops$detection.threshold
    
    list(grid.points = grid.points, pot = U, density = f, min.inds = ops$min, max.inds = ops$max, 
        bw = bw, min.points = min.points, max.points = max.points, detection.threshold = det.th)
    
}


#' Description: find.optima
#'
#' Detect optima, excluding very local optima below detection.threshold
#'
#'  Arguments:
#'    @param f density
#'    @param detection.threshold detection threshold for peaks
#'    @param bw bandwidth
#'    @param x original data
#'    @param detection.limit ignore maxima that are below this value
#'
#' Returns:
#'    @return A list with the following elements:
#'      min minima
#'      max maxima
#'
#' @export
#'
#' @references See citation('TBA') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples #
#'
#' @keywords utilities

find.optima <- function(f, detection.threshold = 0, bw, x, detection.limit = 0) {
    
    # detection.multiplier = detection.threshold Set detection threshold as function
    # of kernel height kernel.height <- dnorm(0, sd = bw) / nrow(x)
    # detection.threshold <- detection.multiplier * kernel.height detlim <-
    # detection.limit * kernel.height detlim <- detection.limit*1/diff(range(x))
    detlim <- detection.limit
    
    # Detect minima and maxima of the density (see Livina et al.)  these correspond
    # to maxima and minima of the potential, respectively including end points of the
    # vector
    maxima <- which(diff(sign(diff(c(-Inf, f, -Inf)))) == -2)
    minima <- which(diff(sign(diff(c(Inf, -f, Inf)))) == -2)
    # minima <- which(extract(turnpoints(ts(f))) == -1) maxima <-
    # which(extract(turnpoints(ts(f))) == 1)
    
    # remove maxima that are below detection threshold
    maxima <- maxima[f[maxima] >= detlim]
    
    # remove minima that now became obsolete If there are multiple minima between two
    # consecutive maxima after removing the maxima that did not pass the threshold,
    # take the average of the minima;return the list of indices such that between
    # each pair of consecutive maxima, there is exactly one minimum
    
    if (length(maxima) > 1) {
        minima <- sapply(2:length(maxima), function(i) {
            
            mins <- minima[minima >= maxima[[i - 1]] & minima <= maxima[[i]]]
            if (length(mins) > 0) {
                round(mean(mins[which(f[mins] == min(f[mins]))]))
            } else {
                NULL
            }
        })
        
    } else {
        minima <- NULL
    }
    
    minima <- unlist(minima)
    maxima <- maxima[f[maxima] >= detlim]
    
    # Combine maxima that do not have minima in between
    if (length(maxima) > 1) {
        maxima2 <- c()
        
        for (i in 1:(length(maxima) - 1)) {
            nominima <- TRUE
            cnt <- 0
            while (nominima & (i + cnt) < length(maxima)) {
                cnt <- cnt + 1
                nominima <- sum(minima > maxima[[i]] & minima < maxima[[i + cnt]]) == 
                  0
                # if (is.na(nominima)) {nominima <- TRUE}
            }
            
            maxs <- maxima[i:(i + cnt - 1)]
            maxima2 <- c(maxima2, round(mean(maxs[which(f[maxs] == max(f[maxs]))])))
            
        }
        
        if (!maxima[[length(maxima)]] %in% maxima2) {
            maxima2 <- c(maxima2, maxima[[length(maxima)]])
        }
        
        
        maxima <- maxima2
        
    }
    
    # Remove minima that are outside the most extreme maxima
    minima <- minima[minima > min(maxima) & minima < max(maxima)]
    
    # Remove minima and maxima that are too shallow
    delmini <- logical(length(minima))
    delmaxi <- logical(length(maxima))
    
    for (j in 1:length(maxima)) {
        
        # Calculate distance of this maximum to all minima
        s <- minima - maxima[[j]]
        
        # Set distances to deleted minima to zero
        s[delmini] <- 0
        
        # identify the closest remaining minima
        i1 <- i2 <- NULL
        if (length(s) > 0) {
            
            minima.spos <- minima[s > 0]
            minima.sneg <- minima[s < 0]
            
            if (length(minima.spos) > 0) {
                i1 <- min(minima.spos)
            }
            if (length(minima.sneg) > 0) {
                i2 <- max(minima.sneg)
            }
            
        }
        
        # if no positive differences available, set it to same value with i2
        if ((is.null(i1) && !is.null(i2))) {
            i1 <- i2
        } else if ((is.null(i2) && !is.null(i1))) {
            # if no negative differences available, set it to same value with i1
            i2 <- i1
        }
        
        if (!is.null(i1) && is.na(i1)) {
            i1 <- NULL
        }
        if (!is.null(i2) && is.na(i2)) {
            i2 <- NULL
        }
        
        # If a closest minimum exists, check differences and remove if difference is
        # under threshold
        if (!is.null(i1)) {
            
            # Smallest difference between this maximum and the closest minima
            diff <- min(c((f[maxima[[j]]] - f[i1]), (f[maxima[[j]]] - f[i2])))
            
            if (diff < detection.threshold) {
                
                # If difference is below threshold, delete this maximum
                delmaxi[[j]] <- TRUE
                
                # Delete the larger of the two neighboring minima
                if (f[[i1]] > f[[i2]]) {
                  delmini[minima == i1] <- TRUE
                } else {
                  delmini[minima == i2] <- TRUE
                }
            }
            
        } else {
            # if both i1 and i2 are NULL, do nothing
        }
    }
    
    # Delete the shallow minima and maxima
    if (length(minima) > 0 && sum(delmini) > 0) {
        minima <- minima[!delmini]
    }
    
    if (length(maxima) > 0 && sum(delmaxi) > 0) {
        maxima <- maxima[!delmaxi]
    }
    
    list(min = minima, max = maxima, detection.threshold = detection.threshold, detection.limit = detlim)
    
}


#' Moving Average Potential
#'
#' This function reconstructs a potential derived from data along a gradient of a given parameter.
#'
#' Arguments:
#'  @param X a vector of the X observations of the state variable of interest
#'  @param param parameter values corresponding to the observations in X 
#'  @param bw Bandwidth for smoothing kernels. Automatically determined by default.
#'  @param bw.adjust Bandwidth adjustment constant
#'  @param detection.threshold Threshold for local optima to be discarded.
#'  @param std Standard deviation.
#'  @param grid.size number of evaluation points; number of steps between min and max potential; also used as kernel window size
#'  @param plot.cutoff cuttoff for potential minima and maxima in visualization
#'  @param plot.contours Plot contours on the landscape visualization
#'  @param binwidth binwidth for contour plot
#'  @param bins bins for contour plot. Overrides binwidth if given
#'
#'  @return A list with the following elements:
#'     pars values of the covariate parameter as matrix;
#'     xis values of the x as matrix;
#'     pots smoothed potentials;
#'     mins minima in the densities (-potentials; neglecting local optima);
#'     maxs maxima in densities (-potentials; neglecting local optima);
#'     plot an object that displays the potential estimated in 2D
#' 
#' @export
#'
#' @references Hirota, M., Holmgren, M., van Nes, E.H. & Scheffer, M. (2011). Global resilience of tropical forest and savanna to critical transitions. \emph{Science}, 334, 232-235.
#' @author L. Lahti, E. van Nes, V. Dakos.
#' @seealso \code{\link{generic_ews}}; \code{\link{ddjnonparam_ews}}; \code{\link{bdstest_ews}}; \code{\link{sensitivity_ews}};\code{\link{surrogates_ews}}; \code{\link{ch_ews}}; \code{livpotential_ews}
# ; \code{\link{timeVAR_ews}}; \code{\link{thresholdAR_ews}}
#' @examples X = c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), rnorm(1000, mean = 2)); param = seq(0,5,length=3000); res <- movpotential_ews(X, param)
#' @keywords early-warning

movpotential_ews <- function(X, param = NULL, bw = "nrd", bw.adjust = 1, detection.threshold = 0.1, 
    std = 1, grid.size = 50, plot.cutoff = 0.5, plot.contours = TRUE, binwidth = 0.2, 
    bins = NULL) {
    
    if (is.null(param)) {
        param <- seq(1, length(X), 1)
    }
    
    nas <- is.na(param) | is.na(X)
    if (sum(nas) > 0) {
        warning("The data contains NAs, removing the associated samples from X and param input arguments.")
        X <- X[!nas]
        param <- param[!nas]
    }
    
    minparam <- min(param)
    maxparam <- max(param)
    
    # Determine step size
    sdwindow <- step <- (maxparam - minparam)/grid.size
    
    # Place evaluation points evenly across data range
    xi <- seq(min(X), max(X), length = grid.size)
    
    # Initialize
    xis <- pars <- pots <- matrix(0, nrow = grid.size, ncol = length(xi))
    maxs <- mins <- matrix(0, nrow = grid.size, ncol = length(xi))
    
    for (i in 1:grid.size) {
        
        # Increase the parameter at each step
        par <- minparam + (i - 0.5) * step
        
        # Check which elements in evaluation range (param) are within 2*sd of par
        weights <- exp(-0.5 * (abs(par - param)/sdwindow)^2)
        
        # LL: Normalization was added in the R implementation 16.5.2012
        weights <- weights/sum(weights)
        
        # Calculate the potential
        tmp <- livpotential_ews(x = X, std = std, bw = bw, bw.adjust = bw.adjust, 
            weights = weights, grid.size = grid.size)
        
        # Store variables
        pots[i, ] <- tmp$pot
        xis[i, ] <- tmp$grid.points
        pars[i, ] <- par + rep(0, length(tmp$grid.points))
        mins[i, tmp$min.inds] <- 1
        maxs[i, tmp$max.inds] <- 1
        
    }
    
    res <- list(pars = pars, xis = xis, pots = pots, mins = mins, maxs = maxs, std = std)
    
    p <- PlotPotential(res, title = "Moving Average Potential", "parameter/time", 
        "state variable", cutoff = plot.cutoff, plot.contours = plot.contours, binwidth = binwidth, 
        bins = bins)
    
    list(res = res, plot = p)
    
}


 
