predict.mcsimex <-
function (object, newdata, ...)
{
    new.object <- object$model
    new.object$coefficients <- object$coefficients
if (missing(newdata)) {
predict(new.object, ...)
} else {
predict(new.object, newdata = data.frame(newdata), ...)
}
}

