gamcen.11 <-
function (x, mu = 0, sig = 1, eta = 0, kappa = 1, alf = 2, bet = 2) 
{
    r <- (x - mu)/sig
    1/sig^2 * (s.1.m.singly.censored(r, alf, bet) - s.0.m.singly.censored(r, 
        alf, bet)^2)
}
