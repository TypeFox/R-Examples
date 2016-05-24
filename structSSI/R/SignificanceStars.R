SignificanceStars <- function(alpha, pvalues) {
        result <- rep(NA, length(pvalues))
        result[pvalues < alpha / 50] <- '***'
        result[pvalues >= alpha / 50 & pvalues < alpha / 5] <- '**'
        result[pvalues >= alpha / 5 & pvalues < alpha] <- '*'
        result[pvalues >= alpha & pvalues < 2 * alpha] <- '.'
        result[pvalues > 2 * alpha] <- '-'
        return(result)
    }
