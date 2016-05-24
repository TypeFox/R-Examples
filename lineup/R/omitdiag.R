## omitdiag.R
## Karl W Broman

#' Replace the diagonal in a distance matrix with missing values
#'
#' Replace the diagonal (that is, self-self distances) from a distance matrix
#' calculated by \code{\link{distee}} or \code{\link{disteg}} with missing
#' values (so that only self-nonself distances are left).
#'
#' We use the row and column names to identify which entries are self-self.
#'
#' @param d A distance matrix calculated by \code{\link{distee}} or
#' \code{\link{disteg}}.
#' @return A matrix of the same form as the input, but with self-self distances
#' replaced with \code{NA}.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{pulldiag}}, \code{\link{distee}}, \code{\link{disteg}},
#' \code{\link{summary.lineupdist}}, \code{\link{plot2dist}},
#' \code{\link{plot.lineupdist}}
#' @keywords array
#' @examples
#' data(expr1, expr2)
#'
#' \dontshow{expr1 <- expr1[,1:500]
#' expr2 <- expr2[,1:500]}
#'
#' # distance as RMS difference
#' d <- distee(expr1, expr2)
#'
#' # focus on the self-nonself distances
#' # (replace self-self distances with NA)
#' d_selfnonself <- omitdiag(d)
#'
#' @export
omitdiag <-
    function(d)
{
    rn <- rownames(d)
    cn <- colnames(d)
    m <- match(rn, cn)
    wh <- which(!is.na(m))
    m <- m[!is.na(m)]
    for(i in seq(along=wh))
        d[wh[i],m[i]] <- NA

    d
}
