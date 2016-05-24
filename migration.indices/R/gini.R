#' Check Migration Matrix
#'
#' Checks if provided R object looks like a migration matrix.
#'
#' A migration matrix is a symmetric matrix with the exact same row and column names. The diagonal equals to zero. The upper triangle shows the in- and the lower triangle shows the out-migration.
#' @param m R object to check
#' @return (invisibly) TRUE
#' @keywords internal
check.migration.matrix <- function(m) {

    ## dummy checks on provided matrix
    if (missing(m))
        stop('No data provided!')
    if (!is.matrix(m))
        stop('Wrong data type (!matrix) provided!')
    if (nrow(m) != ncol(m))
        stop('Wrong data tpye (!symmetrical matrix) provided!')
    if (!is.numeric(m))
        stop('Wrong data tpye (!numeric) provided!')
    if (length(which(is.na(m[xor(upper.tri(m), lower.tri(m))]))) > 0)
        stop('Missing values (outside of diagonal) found in provided matrix!')
    if (!any(all(is.na(diag(m))), all(diag(m) == 0)))
        stop('Diagonal should be zero or missing.')

    return(invisible(TRUE))

}


#' Total Flows Gini Index
#'
#' The Total Gini Index shows the overall concentration of migration with a simple number computed by comparing each cell of the migration matrix with every other cell except for the diagonal:
#' \deqn{G^T = \frac{\sum_i \sum_{j \neq i} \sum_k \sum_{l \neq k} | M_{ij} - M_{kl} | }{ (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}}
#' This implementation solves the above formula by a simple loop for performance issues to compare all values to the others at one go, although smaller migration matrices could also be addressed by a much faster \code{dist} method. Please see the sources for more details.
#' @param m migration matrix
#' @param corrected Bell et al. (2002) updated the formula of Plane and Mulligan (1997) to have \eqn{2{n(n-1)-1}} instead of \eqn{2n(n-1)} in the denominator to "ensure that the index can assume the upper limit of 1".
#' @return A number between 0 and 1 where 0 means no spatial focusing and 1 shows that all migrants are found in one single flow.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.total(migration.hyp)           # 0.2666667
#' migration.gini.total(migration.hyp2)          # 0.225
#' migration.gini.total(migration.hyp, FALSE)    # 0.2222222
#' migration.gini.total(migration.hyp2, FALSE)   # 0.1875
#' @export
#' @seealso \code{\link{migration.gini.col}} \code{\link{migration.gini.row}} \code{\link{migration.gini.exchange}} \code{\link{migration.gini.in}} \code{\link{migration.gini.out}}
migration.gini.total <- function(m, corrected = TRUE) {

    check.migration.matrix(m)

    n           <- nrow(m)
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]

    if (corrected)
        denominator <- 2*(n*(n-1)-1)
    else
        denominator <- 2*n*(n-1)

    return(sum(apply(as.data.frame(m.val), 1, function(x) sum(abs(m.val-x))), na.rm = TRUE)/(denominator * sum(m, na.rm = TRUE)))

    ## faster method (fails with "low memory")
    diag(m)     <- NA
    sum(dist(as.vector(m)), na.rm = TRUE)*2/(denominator * sum(m, na.rm = TRUE))

}


#' Rows Gini Index
#'
#' The Rows Gini index concentrates on the "relative extent to which the destination selections of out-migrations are spatially focused":
#' \deqn{G^T_R = \frac{\sum_i \sum_{j \neq i} \sum_{h \neq i,j} | M_{ij} - M_{ih} | }{ (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}}
#' This implementation solves the above formula by computing the \code{dist} matrix for each row.
#' @param m migration matrix
#' @return A number between 0 and 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.row(migration.hyp)  # 0
#' migration.gini.row(migration.hyp2) # 0.02083333
#' @export
#' @seealso \code{\link{migration.gini.col}} \code{\link{migration.gini.row.standardized}}
migration.gini.row <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    sum(apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))/(2*n*(n-1)*sum(m, na.rm = TRUE))

}


#' Standardized Rows Gini Index
#'
#' The standardized version of the Rows Gini Index (\code{\link{migration.gini.row}}) by dividing that with the Total Flows Gini Index (\code{\link{migration.gini.total}}):
#' \deqn{G^{T*}_R = 100\frac{G^T_R}{G^T}}
#' As this index is standardized, it "facilitate comparisons from one period to the next of the rows" indices.
#' @param m migration matrix
#' @param gini.total optionally pass the pre-computed Total Flows Gini Index to save computational resources
#' @return A percentage range from 0\% to 100\% where 0\% means that the migration flows are uniform, while a higher value indicates spatial focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.row.standardized(migration.hyp)     # 0
#' migration.gini.row.standardized(migration.hyp2)    # 11.11111
#' @export
#' @seealso \code{\link{migration.gini.row}} \code{\link{migration.gini.col.standardized}}
migration.gini.row.standardized <- function(m, gini.total = migration.gini.total(m, FALSE)) {

    100 * migration.gini.row(m) / gini.total
}


#' Columns Gini Index
#'
#' The Columns Gini index concentrates on the "relative extent to which the destination selections of in-migrations are spatially focused":
#' \deqn{G^T_R = \frac{\sum_j \sum_{i \neq j} \sum_{g \neq i,j} | M_{ij} - M_{gj} | }{ (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}}
#' This implementation solves the above formula by computing the \code{dist} matrix for each columns.
#' @param m migration matrix
#' @return A number between 0 and 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.col(migration.hyp)  # 0.05555556
#' migration.gini.col(migration.hyp2) # 0.04166667
#' @export
#' @seealso \code{\link{migration.gini.row}} \code{\link{migration.gini.col.standardized}}
migration.gini.col <- function(m) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    sum(apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE)*2))/(2*n*(n-1)*sum(m, na.rm = TRUE))

}


#' Standardized Columns Gini Index
#'
#' The standardized version of the Columns Gini Index (\code{\link{migration.gini.col}}) by dividing that with the Total Flows Gini Index (\code{\link{migration.gini.total}}):
#' \deqn{G^{T*}_C = 100\frac{G^T_C}{G^T}}
#' As this index is standardized, it "facilitate comparisons from one period to the next" of the columns indices.
#' @param m migration matrix
#' @param gini.total optionally pass the pre-computed Total Flows Gini Index to save computational resources
#' @return A percentage range from 0\% to 100\% where 0\% means that the migration flows are uniform, while a higher value indicates spatial focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.col.standardized(migration.hyp)     # 25
#' migration.gini.col.standardized(migration.hyp2)    # 22.22222
#' @export
#' @seealso \code{\link{migration.gini.col}} \code{\link{migration.gini.row.standardized}}
migration.gini.col.standardized <- function(m, gini.total = migration.gini.total(m, FALSE)) {

    100 * migration.gini.col(m) / gini.total
}


#' Exchange Gini Index
#'
#' The Exchange Gini Index "indicates the contribution to spatial focusing represented by the \eqn{n(n-q)} net interchanges in the system":
#' \deqn{G^T_{RC, CR} = \frac{\sum_i \sum_{j \neq i} | M_{ij} - M_{ji} | }{ (2n(n-1)-1) \sum_i \sum_{j \neq i} M_{ij}}}
#' This implementation solves the above formula by simply substracting the transposed matrix's values from the original one at one go.
#' @param m migration matrix
#' @return A number between 0 and 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @export
#' @examples
#' data(migration.hyp)
#' migration.gini.exchange(migration.hyp)     # 0.05555556
#' migration.gini.exchange(migration.hyp2)    # 0.04166667
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.exchange.standardized}}
migration.gini.exchange <- function(m) {

    check.migration.matrix(m)

    n           <- nrow(m)
    m.t         <- t(m)
    m.val       <- m[xor(upper.tri(m), lower.tri(m))]
    m.t.val     <- m.t[xor(upper.tri(m.t), lower.tri(m.t))]

    sum(abs(m.val - m.t.val))/(2*n*(n-1)*sum(m))

}


#' Standardized Exchange Gini Index
#'
#' The standardized version of the Exchange Gini Index (\code{\link{migration.gini.exchange}}) by dividing that with the Total Flows Gini Index (\code{\link{migration.gini.total}}):
#' \deqn{G^{T*}_{RC, CR} = 100\frac{G^T_{RC, CR}}{G^T}}
#' As this index is standardized, it "facilitate comparisons from one period to the next" of the exchange indices.
#' @param m migration matrix
#' @param gini.total optionally pass the pre-computed Total Flows Gini Index to save resources
#' @return A percentage range from 0\% to 100\% where 0\% means that the migration flows are uniform, while a higher value indicates spatial focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.exchange.standardized(migration.hyp)  # 25
#' migration.gini.exchange.standardized(migration.hyp2) # 22.22222
#' @export
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.exchange}}
migration.gini.exchange.standardized <- function(m, gini.total = migration.gini.total(m, FALSE)) {

    100 * migration.gini.exchange(m) / gini.total

}


#' Out-migration Field Gini Index
#'
#' The Out-migration Field Gini Index is a decomposed version of the Rows Gini Index (\code{\link{migration.gini.row}}) representing "the contribution of each region's row to the total index" () (\code{\link{migration.gini.total}}):
#' \deqn{G^O_i = \frac{\sum_{j \neq i} \sum_{l \neq i,j} | M_{ij} - M_{il} | }{ 2(n-2) \sum_{j \neq k} M_{ij}}}
#' These Gini indices facilitates the direct comparison of different territories without further standardization.
#' @param m migration matrix
#' @param corrected Bell et al. (2002) updated the formula of Plane and Mulligan (1997) to be \eqn{2(n-2)} instead of \eqn{2(n-1)} because "the number of comparisons should exclude the diagonal cell in each row and column, and the comparison of each cell with itself".
#' @return A numeric vector with the range of 0 to 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.out(migration.hyp)         # 0 0 0
#' migration.gini.out(migration.hyp2)        # 0.000 0.25 0.000
#' migration.gini.out(migration.hyp, FALSE)  # 0 0 0
#' migration.gini.out(migration.hyp2, FALSE) # 0.000 0.125 0.000
#' @export
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.in}} \code{\link{migration.weighted.gini.out}}
migration.gini.out <- function(m, corrected = TRUE) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    if (corrected)
        denominator <- 2 * (n - 2)
    else
        denominator <- 2 * (n - 1)

    apply(m, 1, function(m.row) sum(dist(m.row), na.rm = TRUE) * 2) / (denominator * rowSums(m, na.rm = TRUE))

}


#' Migration-weighted Out-migration Gini Index
#'
#' The Migration-weighted Out-migration Gini Index is a weighted version of the Out-migration Field Gini Index (\code{\link{migration.gini.out}}) "according to the zone of destination's share of total migration and the mean of the weighted values is computed as":
#' \deqn{MWG^O = \frac{ \sum_i G^O_i \frac{\sum_j M_{ij}}{\sum_{ij} M_{ij}}}{n}}
#' @param m migration matrix
#' @param mgo optionally passed (precomputed) Migration In-migration Gini Index
#' @references \itemize{
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.weighted.gini.out(migration.hyp)   # 0
#' migration.weighted.gini.out(migration.hyp2)  # 0.02083333
#' @seealso \code{\link{migration.weighted.gini.in}} \code{\link{migration.weighted.gini.mean}}
#' @export
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.out}} \code{\link{migration.weighted.gini.in}} \code{\link{migration.weighted.gini.mean}}
migration.weighted.gini.out <- function(m, mgo = migration.gini.out(m)) {

    diag(m)     <- NA
    n           <- nrow(m)
    m.sum       <- sum(m, na.rm = TRUE)

    sum(mgo * colSums(m, na.rm = TRUE) / m.sum) / n

}


#' In-migration Field Gini Index
#'
#' The In-migration Field Gini Index is a decomposed version of the Columns Gini Index (\code{\link{migration.gini.col}}) representing "the contribution of each region's columns to the total index" () (\code{\link{migration.gini.total}}):
#' \deqn{G^I_j = \frac{\sum_{i \neq j} \sum_{k \neq j,i} | M_{ij} - M_{kj} | }{ 2(n-2) \sum_{i \neq j} M_{ij}}}
#' These Gini indices facilitates the direct comparison of different territories without further standardization.
#' @param m migration matrix
#' @param corrected Bell et al. (2002) updated the formula of Plane and Mulligan (1997) to be \eqn{2(n-2)} instead of \eqn{2(n-1)} because "the number of comparisons should exclude the diagonal cell in each row and column, and the comparison of each cell with itself".
#' @return A numeric vector with the range of 0 to 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.gini.in(migration.hyp)         # 0.2000000 0.5000000 0.3333333
#' migration.gini.in(migration.hyp2)        # 0.2000000 0.0000000 0.4285714
#' migration.gini.in(migration.hyp, FALSE)  # 0.1000000 0.2500000 0.1666667
#' migration.gini.in(migration.hyp2, FALSE) # 0.1000000 0.0000000 0.2142857
#' @export
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.out}} \code{\link{migration.weighted.gini.in}}
migration.gini.in <- function(m, corrected = TRUE) {

    check.migration.matrix(m)

    diag(m)     <- NA
    n           <- nrow(m)

    if (corrected)
        denominator <- 2 * (n - 2)
    else
        denominator <- 2 * (n - 1)

    apply(m, 2, function(m.row) sum(dist(m.row), na.rm = TRUE) * 2) / (denominator * colSums(m, na.rm = TRUE))

}


#' Migration-weighted In-migration Gini Index
#'
#' The Migration-weighted In-migration Gini Index is a weighted version of the In-migration Field Gini Index (\code{\link{migration.gini.in}}) "according to the zone of destination's share of total migration and the mean of the weighted values is computed as":
#' \deqn{MWG^I = \frac{ \sum_j G^I_j \frac{\sum_j M_{ij}}{\sum_{ij} M_{ij}}}{n}}
#' @param m migration matrix
#' @param mgi optionally passed (precomputed) Migration In-migration Gini Index
#' @references \itemize{
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.weighted.gini.in(migration.hyp)   # 0.1222222
#' migration.weighted.gini.in(migration.hyp2)  # 0.05238095
#' @seealso \code{\link{migration.gini}} \code{\link{migration.gini.in}} \code{\link{migration.weighted.gini.out}} \code{\link{migration.weighted.gini.mean}}
#' @export
migration.weighted.gini.in <- function(m, mgi = migration.gini.in(m)) {

    diag(m)     <- NA
    n           <- nrow(m)
    m.sum       <- sum(m, na.rm = TRUE)

    sum(mgi * rowSums(m, na.rm = TRUE) / m.sum) / n

}

#' Migration-weighted Mean Gini Index
#'
#' The Migration-weighted Mean Gini Index is simply the average of the  Migration-weighted In-migration (\code{\link{migration.weighted.gini.in}}) and the Migration-weighted Out-migration (\code{\link{migration.weighted.gini.out}}) Gini Indices:
#' \deqn{MWG^A = \frac{MWG^O + MWG^I}{2}}
#' @param m migration matrix
#' @param mwgi optionally passed (precomputed) Migration-weighted In-migration Gini Index
#' @param mwgo optionally passed (precomputed) Migration-weighted Out-migration Gini Index
#' @return This combined index results in a number between 0 and 1 where 0 means no spatial focusing and 1 shows maximum focusing.
#' @references \itemize{
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @examples
#' data(migration.hyp)
#' migration.weighted.gini.mean(migration.hyp)  # 0.06111111
#' migration.weighted.gini.mean(migration.hyp2) # 0.03660714
#' @seealso \code{\link{migration.weighted.gini.in}} \code{\link{migration.weighted.gini.out}}
#' @export
migration.weighted.gini.mean <- function(m, mwgi, mwgo) {

    if (missing(mwgi))
        mwgi <- migration.weighted.gini.in(m)
    if (missing(mwgo))
        mwgo <- migration.weighted.gini.out(m)

    (mwgi + mwgo) / 2

}


#' Spatial Gini Indexes
#'
#' This is a wrapper function computing all the following Gini indices:
#' \itemize{
#'   \item Total Flows Gini Index (\code{\link{migration.gini.total}})
#'   \item Rows Gini Index (\code{\link{migration.gini.row}})
#'   \item Standardized Rows Gini Index (\code{\link{migration.gini.row.standardized}})
#'   \item Columns Gini Index (\code{\link{migration.gini.col}})
#'   \item Standardized Columns Gini Index (\code{\link{migration.gini.col.standardized}})
#'   \item Exchange Gini Index (\code{\link{migration.gini.exchange}})
#'   \item Standardized Exchange Gini Index (\code{\link{migration.gini.exchange.standardized}})
#'   \item Out-migration Field Gini Index (\code{\link{migration.gini.out}})
#'   \item Migration-weighted Out-migration Gini Index (\code{\link{migration.weighted.gini.out}})
#'   \item In-migration Field Gini Index (\code{\link{migration.gini.in}})
#'   \item Migration-weighted In-migration Gini Index (\code{\link{migration.weighted.gini.in}})
#'   \item Migration-weighted Mean Gini Index (\code{\link{migration.weighted.gini.mean}})
#' }
#' @return List of all Gini indices.
#' @param m migration matrix
#' @param corrected to use Bell et al. (2002) updated formulas instead of Plane and Mulligan (1997)
#' @examples
#' data(migration.hyp)
#' migration.gini(migration.hyp)
#' migration.gini(migration.hyp2)
#' @export
#' @references \itemize{
#'   \item David A. Plane and Gordon F. Mulligan (1997) Measuring Spatial Focusing in a Migration System. \emph{Demography} \bold{34}, 251--262
#'   \item M. Bell, M. Blake, P. Boyle, O. Duke-Williams, P. Rees, J. Stillwell and G. Hugo (2002) Cross-National Comparison of Internal Migration. Issues and Measures. \emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)} \bold{165}, 435--464
#' }
#' @seealso \code{\link{migration.gini.col}} \code{\link{migration.gini.row}} \code{\link{migration.gini.exchange}} \code{\link{migration.gini.in}} \code{\link{migration.gini.out}}
migration.gini <- function(m, corrected = TRUE) {

    res <- list(
             migration.gini.total         = migration.gini.total(m, corrected),
             migration.gini.exchange      = migration.gini.exchange(m),
             migration.gini.row           = migration.gini.row(m),
             migration.gini.col           = migration.gini.col(m),
             migration.gini.in            = migration.gini.in(m, corrected),
             migration.gini.out           = migration.gini.out(m, corrected)
      )

    res$migration.gini.row.standardized           <- migration.gini.row.standardized(m, res$migration.gini.total)
    res$migration.gini.col.standardized           <- migration.gini.col.standardized(m, res$migration.gini.total)
    res$migration.gini.exchange.standardized      <- migration.gini.exchange.standardized(m, res$migration.gini.total)
    res$migration.gini.in.weighted                <- migration.weighted.gini.in(m, res$migration.gini.in)
    res$migration.gini.out.weighted               <- migration.weighted.gini.out(m, res$migration.gini.out)
    res$migration.gini.mean.weighted              <- migration.weighted.gini.mean(m, res$migration.gini.in.weighted, res$migration.gini.out.weighted)

    class(res) <- 'migration.gini'
    return(res)

}


#' @method print migration.gini
#' @S3method print migration.gini
print.migration.gini <- function(x, ...) {

    cat('\n')
    cat('Total Flows Gini Index:             ', x$migration.gini.total, '\n')
    cat('Rows Gini Index:                    ', x$migration.gini.row, '\n')
    cat('Standardized Rows Gini Index:       ', x$migration.gini.row.standardized, '\n')
    cat('Columns Gini Index:                 ', x$migration.gini.col, '\n')
    cat('Standardized Columns Gini Index:    ', x$migration.gini.col.standardized, '\n')
    cat('Exchange Gini Index:                ', x$migration.gini.exchange, '\n')
    cat('Standardized Exchange Gini Index:   ', x$migration.gini.exchange.standardized, '\n')
    cat('In-migration Field Gini Index:      ', 'vector', '\n')
    cat('Weighted In-migration Gini Index:   ', x$migration.gini.in.weighted, '\n')
    cat('Out-migration Field Gini Index:     ', 'vector', '\n')
    cat('Weighted Out-migration Gini Index:  ', x$migration.gini.out.weighted, '\n')
    cat('Migration-weighted Mean Gini Index: ', x$migration.gini.mean.weighted, '\n')
    cat('\n')

}
