## coi_um.R

#' Estimate the coincidence as a function of micron distance
#'
#' Estimate the coincidence as a function of micron distance, with
#' data on XO locations in microns plus SC length in microns.
#'
#' The coincidence function is the probability of a recombination
#' event in both of two intervals, divided by the product of the two
#' intensity function for the two intervals.
#'
#' We estimate this as a function of the distance between the two
#' intervals in microns, taking account of varying SC lengths,.
#'
#' @param xoloc list of crossover locations (in microns) for each of several oocytes or spermatocytes.
#' @param sclength vector of SC lengths (in microns).
#' @param centromeres vector of centromere locations (in microns). If missing, taken to be \code{sclength/2}.
#' @param group nominal vector of groups; the intensity function of
#' the crossover process will be estimated separately for each group,
#' but a joint coincidence function will be estimated.
#' @param intwindow Window size used to smooth the estimated intensity
#' function.
#' @param coiwindow Window size used to smooth the estimated
#' coincidence function.
#' @param intloc Locations at which to estimate the intensity
#' function, in the interval [0,1]
#' @param coiloc Values at which the coincidence function is to be
#' estimated, in microns, less than \code{max(sclength)}
#' @return A list containing the estimated coincidence (as a matrix
#' with two columns, micron distance and corresponding estimated
#' coincidence) and the estimated intensity functions (as a matrix
#' with \code{length(group)+1} columns (the locations at which the
#' intensity functions were estimated followed by the group-specific estimates).
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{gammacoi}}, \code{\link{stahlcoi}},
#' \code{\link{kfunc}}, \code{\link{est.coi}}
#' @keywords models
#' @examples
#' # simple example using data simulated with no crossover interference
#' ncells <- 1000
#' L <- 2                      # chr lengths in Morgans (constant here)
#' nchi <- rpois(ncells, 2*L)  # number of chiasmata
#' xoloc <- lapply(nchi, function(a) runif(a, 0, L)) # chi locations
#' coi <- est.coi.um(xoloc, rep(L, ncells))
#'
#' # plot estimated coincidence and intensity
#' #    (intensity is after scaling chromosome to length 1)
#' par(mfrow=c(2,1), las=1)
#' plot(coi$coincidence, type="l", lwd=2, ylim=c(0, max(coi$coincidence[,2])))
#' plot(coi$intensity, type="l", lwd=2, ylim=c(0, max(coi$intensity[,2])))
#'
#' @useDynLib xoi
#' @export
est.coi.um <-
    function(xoloc, sclength, centromeres, group, intwindow=0.05, coiwindow,
             intloc, coiloc)
{
    # check inputs
    stopifnot(length(xoloc) == length(sclength))
    stopifnot(all(!is.na(unlist(xoloc))))
    stopifnot(all(!is.na(sclength) & sclength >= 0))
    for(i in seq(along=xoloc))
        stopifnot(all(xoloc[[i]] >= 0 & xoloc[[i]] <= sclength[i]))

    if(missing(centromeres)) {
        centromeres <- sclength/2
    } else {
        stopifnot(length(centromeres) == length(xoloc))
        stopifnot(all(!is.na(centromeres) & centromeres > 0 & centromeres < sclength))
    }

    if(missing(group)) {
        group <- rep(1, length(xoloc))
        ugroup <- "intensity"
    } else {
        stopifnot(length(group) == length(xoloc))
        ugroup <- unique(group)
        group <- match(group, ugroup) # turn into integers
    }

    if(min(table(group)) < 2)
        stop("At least one group with < 2 individuals")

    stopifnot(intwindow > 0 && intwindow < 1)
    if(missing(coiwindow)) coiwindow <- min(sclength)/10
    stopifnot(coiwindow > 0 && coiwindow < max(sclength))

    if(missing(intloc)) {
        intloc <- seq(0, 1, length=501)
    } else {
        intloc <- sort(intloc)
        stopifnot(length(intloc) > 0)
        stopifnot(all(!is.na(intloc) & intloc >= 0 & intloc <= 1))
    }

    if(missing(coiloc)) {
        coiloc <- seq(0, min(sclength)-coiwindow/2, length=501)
    } else {
        coiloc <- sort(coiloc)
        stopifnot(length(coiloc) > 0)
        stopifnot(all(!is.na(coiloc) & coiloc >= 0 & coiloc <= max(sclength)))
    }
    # end checks of inputs

    # do the estimation
    z <- .C("R_est_coi_um",
            as.integer(length(xoloc)),
            as.double(unlist(xoloc)),
            as.integer(sapply(xoloc, length)),
            as.double(sclength),
            as.double(centromeres),
            as.integer(group),
            as.integer(max(group)),
            as.double(intwindow),
            as.double(coiwindow),
            as.double(intloc),
            as.integer(length(intloc)),
            as.double(coiloc),
            as.integer(length(coiloc)),
            intensity=as.double(rep(0, length(intloc)*max(group))),
            coincidence=as.double(rep(0, length(coiloc))),
            PACKAGE="xoi")

    # reformat the results
    result <- list(coincidence = cbind(distance=coiloc, coincidence=z$coincidence),
                   intensity = cbind(position=intloc,
                   matrix(z$intensity, ncol=max(group))))
    colnames(result$intensity)[-1] <- ugroup

    result
}
