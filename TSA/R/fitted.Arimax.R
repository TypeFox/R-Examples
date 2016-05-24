`fitted.Arima` <-
function (object,...) 
{
fitted=eval(object$call$x)-residuals(object)
fitted
}

