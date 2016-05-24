fix0111 <-
function (converge = 0.001, d = -1, alpha = 2, beta = 2, mu = 0, 
    sigma = 1, eta = 0, kappa = 1) 
{
    a <- 0.33/sigma^2
    b <- 3/sigma^2
    fa <- fix0103(d, alpha, beta, mu, sigma, eta, kappa = a^(-0.5))
    fb <- fix0103(d, alpha, beta, mu, sigma, eta, kappa = b^(-0.5))
    for (i in 1:1000) {
        if (fb != fa) 
            x <- a - (fa * (b - a))/(fb - fa)
        else x <- (a + b)/2
        if (x <= 0) 
            x <- i/(100 * sigma^2)
        vec.eta <- fix0104(fix0102, d, alpha, beta, mu, sigma, 
            eta, kappa = x^(-0.5))
        fx <- fix0103(d, alpha, beta, mu, sigma, vec.eta[1], 
            kappa = x^(-0.5))
        if (abs(fx) < converge) 
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
    c(eta = vec.eta[1], exp.eta = vec.eta[2], kappa = x^(-0.5), 
        exp.kappa = fx)
}
