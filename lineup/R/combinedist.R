## combinedist.R
## Karl W Broman

# combinedist
#
#' Combine distance matrices into a single such
#'
#' Combine multiple distance matrices into a single distance matrix providing
#' an overall summary
#'
#' The row and column names of the input distance matrices define the
#' individual IDs.
#'
#' If the input distance matrices all have an attribute \code{"denom"} (for
#' denominator) and \code{method="mean"}, we use a weighted mean, weighted by
#' the denominators.  This could be used to calculate an overall proportion.
#'
#' @param \dots Set of distance matrices, as calculated by \code{\link{distee}}
#' or \code{\link{disteg}}.
#' @param method Indicates whether to summarize using the median or the mean.
#' @return A distance matrix, with class \code{"lineupdist"}.  The individual
#' IDs are in the row and column names.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{distee}}, \code{\link{disteg}},
#' \code{\link{summary.lineupdist}}
#' @keywords utilities
#' @examples
#' library(qtl)
#'
#' # load example data
#' data(f2cross, expr1, expr2, pmap, genepos)
#' \dontshow{
#' keep <- c(1:20, 197, 553, 573, 740, 794, 822, 1474, 1522,
#'           1591, 1645, 2080, 2643, 2984, 3089, 3672, 4010, 4039,
#'           4159, 4191, 4198, 4213, 4401, 4544, 4593, 4925)
#' expr1 <- expr1[,keep]
#' expr2 <- expr2[,keep]
#' genepos <- genepos[keep,]}
#'
#' # calculate QTL genotype probabilities
#' f2cross <- calc.genoprob(f2cross, step=1)
#'
#' # find nearest pseudomarkers
#' pmark <- find.gene.pseudomarker(f2cross, pmap, genepos)
#'
#' # line up individuals
#' id1 <- findCommonID(f2cross, expr1)
#' id2 <- findCommonID(f2cross, expr2)
#'
#' # calculate LOD score for local eQTL
#' locallod1 <- calc.locallod(f2cross[,id1$first], expr1[id1$second,], pmark)
#' locallod2 <- calc.locallod(f2cross[,id2$first], expr2[id2$second,], pmark)
#'
#' # take those with LOD > 25
#' expr1s <- expr1[,locallod1>25,drop=FALSE]
#' expr2s <- expr2[,locallod2>25,drop=FALSE]
#'
#' # calculate distance between individuals
#' #     (prop'n mismatches between obs and inferred eQTL geno)
#' d1 <- disteg(f2cross, expr1s, pmark)
#' d2 <- disteg(f2cross, expr2s, pmark)
#'
#' # combine distances
#' d <- combinedist(d1, d2)
#'
#' # summary of problem samples
#' summary(d)
#'
#' @importFrom stats median
#' @export
combinedist <-
    function(..., method=c("median", "mean"))
{
    v <- list(...)

    # input is already a list?
    if(length(v) == 1 && is.list(v[[1]])) v <- v[[1]]

    if(!all(sapply(v, function(a) "lineupdist" %in% class(a))))
        stop("Input distance matrices must each be of class \"lineupdist\".")

    if(length(unique(sapply(v, function(a) class(a)[1]))) > 1)
        stop("Need all of the distance matrices to be the same type.")

    rn <- unique(unlist(lapply(v, rownames)))
    cn <- unique(unlist(lapply(v, colnames)))

    method <- match.arg(method)

    # combine into one big matrix
    d <- array(dim=c(length(rn), length(cn), length(v)))
    dimnames(d) <- list(rn, cn, names(v))
    for(i in seq(along=v))
        d[rownames(v[[i]]),colnames(v[[i]]),i] <- v[[i]]

    if(method=="mean" && all(sapply(v, function(a) !is.null(attr(a, "denom"))))) {
        use.denom <- TRUE
        denom <- array(dim=c(length(rn), length(cn), length(v)))
        dimnames(denom) <- list(rn, cn, names(v))
        for(i in seq(along=v))
            denom[rownames(v[[i]]), colnames(v[[i]]), i] <- attr(v[[i]], "denom")
    }
    else use.denom <- FALSE

    # summarize
    if(method=="median")
        ds <- apply(d, 1:2, median, na.rm=TRUE)
    else if(use.denom) {
        denom.sum <- apply(denom, 1:2, sum, na.rm=TRUE)
        ds <- apply(d*denom, 1:2, sum, na.rm=TRUE)/denom.sum
        attr(ds, "denom") <- denom.sum
    }
    else
        ds <- apply(d, 1:2, mean, na.rm=TRUE)

    class(ds) <- class(v[[1]])

    possible.attributes <- c("d.method", "compareWithin")
    for(i in possible.attributes[possible.attributes %in% names(attributes(v[[1]]))])
        attr(ds, i) <- attr(v[[1]], i)

    ds
}
