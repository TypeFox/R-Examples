coef.cv.gcdnet <- function(object, s = c("lambda.1se", 
    "lambda.min"), ...) {
    if (is.numeric(s)) 
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    coef(object$gcdnet.fit, s = lambda, ...)
} 
