simp_int <-
function (x, fx, n.pts = 256, ret = FALSE) 
{
       if (class(fx) == "function") 
        fx = fx(x)
    n.x = length(x)
    if (n.x != length(fx)) 
        stop("Unequal input vector lengths")
    if (n.pts < 64) 
        n.pts = 64
    ap = approx(x, fx, n = 2 * n.pts + 1)
    h = diff(ap$x)[1]
    integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + 
        ap$y[2 * (1:n.pts) + 1])/3
    invisible(list(value = sum(integral), cdf = list(x = ap$x[2 * 
        (1:n.pts)], y = cumsum(integral))))
}
