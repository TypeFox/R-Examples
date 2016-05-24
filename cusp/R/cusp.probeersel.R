`cusp.probeersel` <-
function (alpha, beta) 
{
    w1 = cusp.extrema(alpha, beta)
    w2 = cusp.curve(alpha, beta)
    x = sort(c(w1, range(w2)))
    y = dcusp.unnorm(x, alpha, beta)
    dx = diff(x)
    dy = filter(y, c(1/2, 1/2))[-5]
    0.5 * y[1] * dx[1] + sum(dy * dx) + 0.5 * y[5] * dx[4]
}

