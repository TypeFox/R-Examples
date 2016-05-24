derivative <-
function (fcn, x, epsilon = .Machine$double.eps, curve.scale.min = 0.1, 
    ...) 
{
    if (any(!is.finite(x))) 
        stop("NA/NaN/Inf values not allowed in 'x'")
    xc <- ifelse(abs(x) < curve.scale.min, curve.scale.min, x)
    h <- (epsilon^(1/3)) * xc
    (fcn(x + h, ...) - fcn(x - h, ...))/(2 * h)
}
