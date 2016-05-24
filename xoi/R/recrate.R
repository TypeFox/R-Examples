## recrate.R

#' Estimate recombination rate
#'
#' Obtain a smoothed estimate of the recombination rate along a chromosome,
#' using the cM and Mbp position of markers.
#'
#' We assume constant recombination rate within each marker interval.
#'
#' @param genmap Vector of cM positions of markers, or a list of such vectors.
#' @param phymap Vector of Mbp positions of markers, or a list of such vectors;
#' same length as \code{genmap}.
#' @param pos Vector of positions at which the recombination rate should be
#' estimated, or a list of such vectors.  If missing, we use the physical
#' marker positions plus a grid with 4 positions per Mbp.
#' @param window Length of sliding window (in Mbp).
#' @return A data.frame containing the positions and estimate recombination
#' rates.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{est.coi}}, \code{\link{intensity}}
#' @keywords models
#' @examples
#' # create equally-spaced map
#' pmap <- sim.map(100, n.mar=51, anchor=TRUE, include.x=FALSE, eq.spacing=TRUE)
#'
#' # simulate cross
#' x <- sim.cross(pmap, type="bc", n.ind=501)
#'
#' # estimate map for that cross
#' emap <- est.map(x)
#'
#' # empirical estimate of recombination rate
#' rr <- est.recrate(emap[[1]], pmap[[1]], window=5)
#' plot(rr, type="l", lwd=2)
#'
#' @useDynLib xoi
#' @export
est.recrate <-
    function(genmap, phymap, pos, window=5)
{
    # multiple-chromosome version
    if(is.list(genmap) && is.list(phymap)) {
        if(length(genmap) != length(phymap))
            stop("length(genmap) != length(phymap)")
        if(missing(pos)) {
            result <- apply(mapply(est.recrate, genmap, phymap, window=window), 2, as.data.frame)
        } else {
            if(length(pos) != length(genmap))
                stop("length(pos) != length(genmap)")
            result <- apply(mapply(est.recrate, genmap, phymap, pos, window=window), 2, as.data.frame)
        }
        names(result) <- names(genmap)
        return(result)
    }

    if(length(genmap) != length(phymap))
        stop("genmap and phymap should be the same length.")

    if(any(is.na(genmap) | is.na(phymap))) {
        drop <- is.na(genmap) | is.na(phymap)
        genmap <- genmap[!drop]
        phymap <- phymap[!drop]
        if(length(genmap) == 0)
            stop("Need multiple markers with known genetic and physical positions.")
    }

    if(length(unique(phymap)) != length(phymap))
        stop("Can't have markers right on top of each other (in physical distance).")

    if(window < 0)
        stop("window should be > 0")

    if(missing(pos)) {
        pos <- sort(c(phymap, seq(min(phymap), max(phymap), by=0.25)))
        pos <- unique(pos[pos >= min(phymap) & pos <= max(phymap)])
    }

    if(any(diff(genmap)<0))
        stop("genmap should be non-decreasing")
    if(any(diff(phymap)<0))
        stop("phymap should be non-decreasing")
    if(any(diff(pos)<0))
        stop("pos should be non-decreasing")

    if(any(pos < min(phymap) | pos > max(phymap))) {
        warning("pos should be within range of phymap")
        pos <- pos[pos >= min(phymap) & pos <= max(phymap)]
    }

    if(diff(range(phymap)) < window)
        stop("range of phymap should exceed the window length.")

    rate <- .C("R_est_recrate",
               as.integer(length(genmap)),
               as.double(genmap),
               as.double(phymap),
               as.integer(length(pos)),
               as.double(pos),
               rate=as.double(rep(0,length(pos))),
               as.double(window),
               as.double(rep(0, length(genmap)-1)),
               PACKAGE="xoi")$rate

    data.frame(pos=pos, rate=rate)
}


#' Convert recrate to scanone format
#'
#' Convert the result of \code{\link{est.recrate}} to the format
#' output by R/qtl's \code{\link[qtl]{scanone}} function.
#'
#' @param recrate A list of results from \code{\link{est.recrate}}
#' @param phymap A list of vectors of Mbp positions of markers
#' @return A data frame with class \code{"scanone"}, in the format output by \code{\link[qtl]{scanone}}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{est.recrate}}
#' @keywords models
#'
#' @examples
#' pmap <- sim.map(100, n.mar=51, anchor=TRUE, include.x=FALSE, eq.spacing=TRUE)
#'
#' # simulate cross
#' x <- sim.cross(pmap, type="bc", n.ind=501)
#'
#' # estimate map for that cross
#' emap <- est.map(x)
#'
#' # empirical estimate of recombination rate
#' rr <- est.recrate(emap[[1]], pmap[[1]], window=5)
#'
#' # make it a list (one component per chromosome, but here just the one chromosome)
#' rr <- list("1"=rr)
#'
#' # convert to scanone output and plot
#' rr_scanone <- recrate2scanone(rr)
#' plot(rr_scanone)
#'
#' @useDynLib xoi
#' @export
recrate2scanone <-
    function(recrate, phymap)
{
    if(is.data.frame(recrate) && names(recrate)[1]=="pos" && names(recrate)[2]=="rate")
        recrate <- list(recrate)
    else {
        if(is.null(names(recrate)))
            names(recrate) <- 1:length(recrate)
    }

    npos <- sapply(recrate, nrow)
    chr <- factor(rep(names(recrate), npos), levels=names(recrate))
    pos <- unlist(lapply(recrate, function(a) a$pos))
    rate <- unlist(lapply(recrate, function(a) a$rate))

    locnam <- vector("list", length(recrate))
    for(i in seq(along=recrate))
        locnam[[i]] <- paste("c", names(recrate)[i], ".loc", seq(along=recrate[[i]]$pos), sep="")

    if(!missing(phymap)) {
        for(i in seq(along=recrate)) {
            m <- match(phymap[[i]], recrate[[i]]$pos)
            locnam[[i]][m] <- names(phymap[[i]])
        }
    }

    result <- data.frame(chr=chr, pos=pos, "cM/Mbp"=rate)
    rownames(result) <- unlist(locnam)
    class(result) <- c("scanone", "data.frame")

    result
}
