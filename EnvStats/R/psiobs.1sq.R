psiobs.1sq <-
function (x, mu = 0, sig = 1, eta = 0, kappa = 1, alf = 2, bet = 2) 
{
    r <- (x - mu)/sig
    (1/sig * s.2.m.singly.censored(r, alf, bet))^2 * dnorm(x, 
        mean = eta, sd = kappa)
}
