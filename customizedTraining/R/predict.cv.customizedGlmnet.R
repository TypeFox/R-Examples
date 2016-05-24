predict.cv.customizedGlmnet <-
function(object, ...)
{
    c(predict.customizedGlmnet(object$fit, lambda = object$lambda.min))
}
