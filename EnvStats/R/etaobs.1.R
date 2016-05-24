etaobs.1 <-
function (x, mu = 0, sig = 1, eta = 0, kappa = 1, alf = 2, bet = 1) 
{
    z <- (x - eta)/kappa
    eps <- (x - mu)/sig
    s.2.m.singly.censored(eps, alf, bet) * z * dnorm(x, mean = eta, 
        sd = kappa)
}
