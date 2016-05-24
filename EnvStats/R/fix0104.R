fix0104 <-
function (f, d = -1, alpha = 2, beta = 2, mu = 0, sigma = 1, 
    eta = 0, kappa = 1) 
{
    a <- mu - sigma
    b <- mu + sigma
    fa <- f(d, alpha, beta, mu, sigma, eta = a, kappa)
    fb <- f(d, alpha, beta, mu, sigma, eta = b, kappa)
    for (i in 1:1000) {
        if (fb != fa) 
            x <- a - (fa * (b - a))/(fb - fa)
        else x <- (a + b)/2
        fx <- f(d, alpha, beta, mu, sigma, eta = x, kappa)
        if (abs(fx) < 0.001) 
            break
        if ((sign(fb) != sign(fa) && sign(fx) != sign(fa)) || 
            (sign(fb) == sign(fa) && abs(x - a) < abs(x - b))) {
            b <- x
            fb <- fx
        }
        else {
            a <- x
            fa <- fx
        }
    }
    c(x = x, f = fx)
}
