plot.cv.customizedGlmnet <-
function(x, ...)
{
    plot(x$fit, x$lambda.min)
}
