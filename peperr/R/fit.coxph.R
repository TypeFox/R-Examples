fit.coxph <- function(response, x, cplx,...)
{
require(survival)
actual.data <- as.data.frame(x)
actual.data$time <- response[,"time"]
actual.data$status <- response[,"status"]

res <- coxph(Surv(time, status)~., data=actual.data, ...)
res
}

