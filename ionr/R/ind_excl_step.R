#' One step in indicator exclusion procedure
#'
#' @description See \code{\link{ind_excl}} for details.
#' @inheritParams ind_excl
#' @param exclude Exclude an item excluded at previous step, e.g., as decided by \code{\link{ind_excl_inc}}
#' @param round Allows rounding of values in returned matrix.
#' @return Provides the results of a single step in indicator exclusion procedure. See example for
#' details
#' @encoding utf-8
#'
#' @examples
#' ## Create a scale-outcome set that violates ION. Only 2 indicators out of 8 relate to
#' ## the outcome, the others just relate to the 2 indicators
#' set.seed(466)
#' a<-scale_sim(n=2500, to_n=2, tn_n=6)
#' res=ind_excl_step(a[[1]],a[[2]])
#' print(res)
#'
#' # note that the p-values for upper items (7 & 8 ) are much smaller than for the rest
#'
#' #row number   indicator number
#' #r.test.t     t value of the r.test.
#' #t.test.p     p value of the r.test.
#' #cor_excl     correlation between outcome and sum-score when an item is excluded.
#' #cor_all      correlation between outcome and sum-score when all items are included
#' # (i.e., full scale).
#' #cor.excl_all correlation between two sum-scores.
#' @export


ind_excl_step <- function(indicators, outcome, indicatornames = 1:ncol(indicators), exclude = vector(), coruse = "everything", round = F) {
    # drop indicators from scale, test them with outcome. Commented version supports multiple drop.
    # However, current algorithm makes no use of multiple drop.

    # foo.outcome <- function(i, o, drop = 1, coruse =
    # 'everything') { combs <- combn(ncol(i), drop) res <- numeric(length = ncol(combs)) #better than res = numeric(0)
    # for (j in 1:ncol(combs)) { res[j] <- cor(o, rowMeans(i[, -c(combs[, j,drop = FALSE])]), use = coruse) }
    # return(res) }


    # new version, drops just one
    foo_outcome_one <- function(i, o, coruse = "everything") {
        # assign temp variables and preallocate for speed
        nindic = ncol(i)
        res <- numeric(length = nindic)
        for (j in 1:nindic) {
            res[j] <- cor(o, rowMeans(i[, -j, drop = FALSE], na.rm = T), use = coruse)
        }
        return(res)
    }
    # test if some indicators are excluded from the test. if yes, then remove these indicators
    if (length(exclude) != 0) {
        # get indicators that are excluded
        excl_indicators <- grep(paste0(exclude, collapse = "|"), indicatornames)  #http://stackoverflow.com/questions/7597559/grep-in-r-with-a-list-of-patterns
        # exclude them from current setup
        indicators <- indicators[, -excl_indicators]  # before 160108: indicators=indicators[-excl_indicators]
        indicatornames <- indicatornames[-excl_indicators]
    }

    # calculate correlations with outcome, when each indicator of the scale is excluded
    dat <- foo_outcome_one(i = indicators, o = outcome, coruse = coruse)

    # correlation between all indicators and outcome
    cor_all <- cor(rowMeans(indicators,na.rm=T), outcome, use = coruse)

    # preallocate
    stats <- matrix(nrow = ncol(indicators), ncol = 5)

    # test the correlation difference between outcome and sumscore that has all indicators, vs. 1 indicator excluded
    for (j in 1:length(dat)) {
        cor_excl <- dat[j]
        cor_excl_all <- cor(rowMeans(indicators,na.rm=F), rowMeans(indicators[, -j,drop=FALSE],na.rm=T), use = coruse)
        res <- psych::r.test(n = nrow(indicators), r12 = cor_excl, r13 = cor_all, r23 = cor_excl_all)

        numbers <- (c(res$t, res$p, cor_excl, cor_all, cor_excl_all))

        stats[j, ] <- numbers
    }
    rownames(stats) <- indicatornames

    colnames(stats) <- c("r.test.t", "t.test.p", "cor_excl", "cor_all", "cor_excl_all")

    if (round == T) {

        stats[, 2] <- round(stats[, 2], 3)
        stats[, c(1, 3:5)] <- round(stats[, c(1, 3:5)], 2)
    }

    stats <- stats[order(abs(stats[, 1]), decreasing = T), ]
    return(stats)
}
