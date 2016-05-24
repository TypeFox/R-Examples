noEffect <- function(object)
{
    if (identical(object$"type", "binomial"))
    {
        respVec <- object$"dataList"$resp
        weiVec <- object$"dataList"$weights
        llNull <- logLik(glm(cbind(respVec*weiVec, (1-respVec)*weiVec) ~ 1, family = binomial))
    }
    
    if (identical(object$"type", "Poisson"))
    {
        llNull <- logLik(glm(resp ~ 1, family = poisson))
    }

    if (identical(object$"type", "continuous"))
    {
        llNull <- logLik(lm(object$dataList$resp ~ 1))
    }
    lldrc <- logLik(object)
    lrt <- -2*(llNull - lldrc)
    dfDiff <- attr(lldrc, "df") - attr(llNull, "df")

    retVec <- c(lrt, dfDiff, 1 - pchisq(lrt, dfDiff))
    names(retVec) <- c("Chi-square test", "Df", "p-value")
    retVec
}