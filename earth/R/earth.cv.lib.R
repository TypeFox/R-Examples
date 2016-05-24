# earth.cv.lib.R: library functions for cross validation of earth models

# The following functions were lifted from the Elith Leathwick code.
# In that code, the AUC calculation was adapted from Ferrier, Pearce and Watson.
# See Pearce and Ferrier (2000) Evaluating the predictive performance
# of habitat models developed using logistic regression.

get.binomial.deviance <- function(yhat, y) # yhat is predicted, y is observed
{
    deviance.contribs <- y * log(yhat) + (1 - y) * log(1 - yhat)
    mean(-2 * sum(deviance.contribs)) / length(y) # TODO length(y) ok, should use dof?
}
get.poisson.deviance <- function(yhat, y)
{
    deviance.contribs <- ifelse(y == 0, 0, y * log(y/yhat)) - (y - yhat)
    2 * sum(deviance.contribs) / length(y)
}
get.auc <- function(yhat, y) # area under ROC curve
{
    nx <- length(y[y == 0])
    ny <- length(y[y == 1])
    xy <- c(yhat[y == 0], yhat[y == 1])
    wilc <- nx * ny   +    (nx * (nx + 1)) / 2   -   sum(rank(xy)[seq_len(nx)])
    wilc / (nx * ny)
}
get.binomial.calib <- function(yhat, y) # returns c(intercept, slope)
{
    yhat <- yhat + 1e-005 # prevents log(0)
    yhat[yhat >= 1] <- .99999
    glm(y ~ log(yhat / (1 - yhat)), family = binomial)$coefficients
}
get.poisson.calib <- function(yhat, y) # returns c(intercept, slope)
{
    glm(y ~ log(yhat), family = poisson)$coefficients
}
get.maxerr <- function(errs) # get signed max absolute err; if matrix then of each col
{
    if(NCOL(errs) == 1)
        errs[which.max(abs(errs))]
    else
        apply(errs, 2, function(col)  col[which.max(abs(col))])
}
# get fraction correctly classified (one row of class.rate.tab)

get.class.rate <- function(yhat, y, ylevels)
{
    stopifnot(ncol(yhat) == ncol(y))
    n <- nrow(y)
    if(ncol(y) == 1)    # single response model?
        per.class.correct <- overall.correct <- sum((yhat > .5) == (y > .5))
    else {
        # multiple response model
        # y and yhat are indicator columns, convert back to levels
        y    <- ylevels[apply(y,    1, which.max)]
        yhat <- ylevels[apply(yhat, 1, which.max)]
        per.class.correct <- repl(0, length(ylevels))
        for(i in seq_along(ylevels)) {
            level <- ylevels[i]
            per.class.correct[i] <- sum(y == level & yhat == level) +
                                    sum(y != level & yhat != level)
        }
        overall.correct <- sum(y == yhat)
    }
    c(per.class.correct, overall.correct) / n
}
