refit.simex <-
function (object, fitting.method = "quadratic",
jackknife.estimation = "quadratic", asymptotic = TRUE, ...)
.refit(object, fitting.method = fitting.method,
jackknife.estimation = jackknife.estimation,
asymptotic = asymptotic,
allowed.fitting = c("quad", "line", "nonl"),
allowed.jackknife = c("quad", "line", "nonl", FALSE), ...)

