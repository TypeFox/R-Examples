"coef.HOF" <-
    function (
			object, 
			model, 
			...
) {
    maxNrofParameters <- 5  
    out <- sapply(object$models, function(x) c(x$par, rep(NA, maxNrofParameters - length(x$par))))
    rownames(out) <- letters[1:maxNrofParameters]
    if (!missing(model)) {
        out <- out[,model]
    }
    out
}
"deviance.HOF" <-
    function (object, model, ...) 
{
    out <- sapply(object$models, function(x) x$deviance)
    if (!missing(model))
        out <- out[model]
    out
}
"fitted.HOF" <-
    function (object, model, ...) 
{
    out <- sapply(object$models, function(x) x$fitted)  
    if(!missing(model)) out <- out[,model]
    out
}
"predict.HOF" <-
    function (object, model, newdata, ...) {
    if(missing(model)) model <- pick.model(object, ...)
    p <- coef(object, model, ...)
    xrange <- object$range
    if (missing(newdata)) x <- object$x else x <- newdata
    fv <- HOF.fun(x=x, model=as.character(model), p=as.numeric(p), M=1, xrange)
    fv
}
