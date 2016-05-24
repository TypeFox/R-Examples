predict.coxme <- function(object, newdata, 
                       type=c("lp", "risk")) {
    # This is an early skeleton of the function
    type <-match.arg(type)
    n <- object$n
    Terms <-  object$terms

    if (!missing(newdata)) stop("newdata argument not yet supported")
    
    out <- object$linear.predictor
    if (type=="risk") out <- exp(out)
    if (!is.null(object$na.action))
        napredict(object$na.action, out)
    else out
}
