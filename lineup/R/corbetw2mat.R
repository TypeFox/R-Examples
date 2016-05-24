## corbetw2mat.R
## Karl W Broman

# corbetw2mat
#
#' Calculate correlations between columns of two matrices
#'
#' For matrices x and y, calculate the correlation between columns of x and
#' columns of y.
#'
#' Missing values (\code{NA}) are ignored, and we calculate the correlation
#' using all complete pairs, as in \code{\link[stats]{cor}} with
#' \code{use="pairwise.complete.obs"}.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix with the same number of rows as \code{x}.
#' @param what Indicates which correlations to calculate and return.  See
#' value, below.
#' @param corthresh Threshold on correlations if \code{what="bestpairs"}.
#' @return If \code{what="paired"}, the return value is a vector of
#' correlations, between columns of \code{x} and the corresponding column of
#' \code{y}.  \code{x} and \code{y} must have the same number of columns.
#'
#' If \code{what="bestright"}, we return a data frame of size \code{ncol(x)} by
#' \code{3}, with the \eqn{i}th row being the maximum correlation between
#' column \eqn{i} of \code{x} and a column of \code{y}, and then the
#' \code{y}-column index and \code{y}-column name with that correlation.  (In
#' case of ties, we give the first one.)
#'
#' If \code{what="bestpairs"}, we return a data frame with five columns,
#' containing all pairs of columns (with one in \code{x} and one in \code{y})
#' with correlation \eqn{\ge} \code{corthresh}.  Each row corresponds to a
#' column pair, and contains the correlation and then the \code{x}- and
#' \code{y}-column indices followed by the \code{x}- and \code{y}-column names.
#'
#' If \code{what="all"}, the output is a matrix of size \code{ncol(x)} by
#' \code{ncol(y)}, with all correlations between columns of \code{x} and
#' columns of \code{y}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{distee}}, \code{\link{findCommonID}}
#' @keywords array univar multivariate
#' @examples
#'
#' data(expr1, expr2)
#' \dontshow{expr1 <- expr1[,1:100]
#' expr2 <- expr2[,1:100]}
#'
#' # correlations with paired columns
#' r <- corbetw2mat(expr1, expr2)
#' # top 10, by absolute value
#' r[order(abs(r), decreasing=TRUE)[1:10]]
#'
#' # all pairs of columns with correlation >= 0.8
#' r_allpairs <- corbetw2mat(expr1, expr2, what="bestpairs", corthresh=0.6)
#'
#' # for each column in left matrix, most-correlated column in right matrix
#' r_bestright <- corbetw2mat(expr1, expr2, what="bestright")
#'
#' @useDynLib lineup
#' @export
corbetw2mat <-
    function(x, y, what=c("paired", "bestright", "bestpairs", "all"),
             corthresh=0.9)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    if(!is.matrix(y)) y <- as.matrix(y)

    n <- nrow(x)
    if(nrow(y) != n)
        stop("nrow(x)=", n, ", which is not equal to nrow(y)=", nrow(y))

    px <- ncol(x)
    py <- ncol(y)
    what <- match.arg(what)

    if(is.null(colnames(x))) colnames(x) <- paste("V", 1:ncol(x), sep="")
    if(is.null(colnames(y))) colnames(y) <- paste("V", 1:ncol(y), sep="")

    if(what=="paired" && py != px)
        stop("what=\"paired\", but ncol(x)=", px, ", which is not equal to ncol(y)=", py)

    if(what=="paired") {
        res <- .C("R_corbetw2mat_paired",
                  as.integer(n),
                  as.integer(px),
                  as.double(x),
                  as.double(y),
                  cor=as.double(rep(NA, px)),
                  PACKAGE="lineup",
                  NAOK=TRUE)$cor
        names(res) <- colnames(x)
    }

    else if(what=="bestright") {
        res <- .C("R_corbetw2mat_unpaired_lr",
                  as.integer(n),
                  as.integer(px),
                  as.double(x),
                  as.integer(py),
                  as.double(y),
                  cor=as.double(rep(NA, px)),
                  index=as.integer(rep(NA, px)),
                  PACKAGE="lineup",
                  NAOK=TRUE)
        res <- data.frame(cor=res$cor, yindex=res$index)
        rownames(res) <- colnames(x)
        res <- cbind(res, ycol=colnames(y)[res[,2]])
    }
    else if(what=="bestpairs") {
        res <- .C("R_corbetw2mat_unpaired_best",
                  as.integer(n),
                  as.integer(px),
                  as.double(x),
                  as.integer(py),
                  as.double(y),
                  cor=as.double(rep(NA, px*py)),
                  xindex=as.integer(rep(NA, px*py)),
                  yindex=as.integer(rep(NA, px*py)),
                  numpairs=as.integer(0),
                  as.double(corthresh),
                  PACKAGE="lineup",
                  NAOK=TRUE)
        res <- data.frame(cor=res$cor[1:res$numpairs],
                          xindex=res$xindex[1:res$numpairs],
                          yindex=res$yindex[1:res$numpairs])
        res <- cbind(res, xcol=colnames(x)[res[,2]], ycol=colnames(y)[res[,3]])

    }
    else {
        res <- .C("R_corbetw2mat_unpaired_all",
                  as.integer(n),
                  as.integer(px),
                  as.double(x),
                  as.integer(py),
                  as.double(y),
                  cor=as.double(rep(NA, px*py)),
                  PACKAGE="lineup",
                  NAOK=TRUE)$cor
        res <- matrix(res, nrow=px, ncol=py)
        dimnames(res) <- list(colnames(x), colnames(y))
    }

    res
}
