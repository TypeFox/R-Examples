poi.ci <- function (x, conf.level = 0.95) 
{
    nn <- length(x)
    LCI <- qchisq((1 - conf.level)/2, 2 * sum(x))/2/nn
    UCI <- qchisq(1 - (1 - conf.level)/2, 2 * (sum(x) + 1))/2/nn
    res <- cbind(mean(x), LCI, UCI)
    ci.prefix <- paste(round(100 * conf.level, 1), "%", sep = "")
    colnames(res) <- c("PointEst", paste(ci.prefix, "LCI"), paste(ci.prefix, 
        "UCI"))
    res
}                 