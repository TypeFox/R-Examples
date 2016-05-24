############################################################################################
## package 'secr'
## LRtest.R
## likelihood ratio test for two models (assumed to be nested)
## last changed 2011 12 13 to include openCR
## 2013-10-29 fixed name of p.value; use logLik
############################################################################################

LR.test <- function (model1, model2) {
    if (is.null(getS3method('logLik', class(model1), TRUE)))
        stop ("no logLik method for model")
    if (class(model1) != class(model2))
        stop ("models must have same class")
    statistic <- as.numeric(2 * abs(logLik(model1) - logLik(model2)))
    if (length(statistic) != 1)
        stop ("problem with 'model1' or 'model2'")
    parameter <- abs(length(model1$fit$par) - length(model2$fit$par))
    p.value <- 1 - pchisq(statistic, parameter)
    names(statistic) <- 'X-square'
    names(parameter) <- 'df'
    names(p.value) <- NULL
    s1name <- deparse(substitute(model1))
    s2name <- deparse(substitute(model2))
    temp <- list (
        statistic = statistic,
        parameter = parameter,
        p.value = p.value,
        method = 'Likelihood ratio test for two models',
        data.name = paste (s1name, 'vs', s2name)
    )
    class(temp) <- 'htest'
    temp
}
############################################################################################

