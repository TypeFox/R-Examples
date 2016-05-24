## findCommonID.R
## Karl W Broman

# findCommonID
#
#' Find individuals in common between a cross and a phenotype matrix
#'
#' Identify which individuals are in common between a QTL mapping data set and
#' a matrix of phenotypes, series of genes.
#'
#'
#' @param id1 A character vector of individual IDs.  This can also be a QTL
#' cross object (see \code{\link[qtl]{read.cross}}), in which case
#' \code{\link[qtl]{getid}} is used to grab individual IDs, or a matrix or data
#' frame, in which case the rownames are taken to be IDs.
#' @param id2 Like \code{id1}, can be a character vector, a cross or a
#' matrix/data frame.
#' @return A list with three components:
#'
#' First, a data frame with rows corresponding to all individuals (across the
#' two sets of individual IDs) and three columns: \code{indexInFirst} and
#' \code{indexInSecond} contain numeric indices to the locations of the
#' individuals within \code{cross} and \code{pheno}, and \code{inBoth} is a
#' logical vector to indicate which individuals appear in both crosses.  The
#' row names are the individual identifiers.
#'
#' The second and third components are vectors of indices in \code{id1} and
#' \code{id2}, respectively, indicating the paired locations of the individuals
#' that are in common.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{calc.locallod}}, \code{\link{corbetw2mat}}
#' @keywords utilities
#' @examples
#'
#' data(f2cross, expr1)
#'
#' # align IDs
#' id <- findCommonID(f2cross, expr1)
#'
#' # aligned data
#' f2cross_aligned <- f2cross[,id$first]
#' expr1_aligned <- expr1[id$second,]
#'
#' @export
findCommonID <-
    function(id1, id2)
{
    if("cross" %in% class(id1))
        id1 <- qtl::getid(id1)
    else if(!is.null(rownames(id1)))
        id1 <- rownames(id1)

    if("cross" %in% class(id2))
        id2 <- qtl::getid(id2)
    else if(!is.null(rownames(id2)))
        id2 <- rownames(id2)

    if(is.null(id1))
        stop("Can't find IDs in id1")
    if(is.null(id2))
        stop("Can't find IDs in id2")

    if(is.factor(id1)) id1 <- as.character(id1)
    if(is.factor(id2)) id2 <- as.character(id2)

    m1 <- match(id2, id1)
    m2 <- match(id1, id2)

    totind <- length(id1) + sum(is.na(m1))
    allID <- list(mat=data.frame(indexInFirst=rep(0, totind),
                  indexInSecond=rep(0, totind),
                  inBoth=rep(FALSE,totind)),
                  first=rep(NA, totind),
                  second=rep(NA, totind))

    rownames(allID$mat) <- c(id1, id2[is.na(m1)])

    allID$mat[,1] <- match(rownames(allID$mat), id1)
    allID$mat[,2] <- match(rownames(allID$mat), id2)
    allID$mat[,3] <- !is.na(allID$mat[,1]) & !is.na(allID$mat[,2])

    allID$first <- allID$mat[allID$mat[,3],1]
    allID$second <- allID$mat[allID$mat[,3],2]

    allID
}
