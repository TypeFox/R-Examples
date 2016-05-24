##
## INPUT:
## x: output from profile.game
##
## RETURN:
## vector of parameters giving the highest profiled log-likelihood
##
svalsFromProfile <- function(x)
{
    x <- do.call(rbind, x)
    bestrow <- which.max(x[, 1])  ## highest log-likelihood
    xn <- names(x)[-1]
    ans <- as.numeric(x[bestrow, ][-1])
    names(ans) <- xn
    return(ans)
}
