predictProb.coxph <- function (object, response, x, times, ...)
{
    require(survival)
    newdata <- as.data.frame(x)
    newdata$time <- response[, "time"]
    newdata$status <- response[, "status"]
    survival.survfit.coxph <- getFromNamespace("survfit.coxph",
        ns = "survival")
    survival.summary.survfit <- getFromNamespace("summary.survfit",
        ns = "survival")
    survfit.object <- survival.survfit.coxph(object, newdata = newdata,
        se.fit = FALSE, conf.int = FALSE)
    inflated.pred <- survival.summary.survfit(survfit.object,
        times = times)
    p <- t(inflated.pred$surv)
    if ((miss.time <- (length(times) - NCOL(p))) > 0)
        p <- cbind(p, matrix(rep(NA, miss.time * NROW(p)), nrow = NROW(p)))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop("Prediction failed")
    p
}