fix0103 <-
function (d = -1, alpha = 2, beta = 2, mu = 0, sigma = 1, eta = 0, 
    kappa = 1) 
{
    r <- (d - mu)/sigma
    rt <- (d - eta)/kappa
    g20 <- g2.m.singly.censored(r, 0, alpha, beta)
    g1.5 <- g1.m.singly.censored(r, 0.5, alpha, beta)
    s.r <- g1.5/g20 * (g20 > 0)
    if (g20 == 0) 
        s.r <- 0
    exp.chi <- -2 * alpha + pnorm(rt) * (r * s.r + 2 * alpha) + 
        fix0107(fix0109, d, eta + 10 * kappa, 1000, alpha, beta, 
            mu, sigma, eta, kappa) * (2 * alpha + 1) * 2 * beta
    exp.chi
}
