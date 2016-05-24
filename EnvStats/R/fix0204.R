fix0204 <-
function (func = 1, a = 0, b = 1, n, mu, sig, eta, kappa, alf, 
    bet) 
{
    (sum(func(a + ((1:n - 0.5) * (b - a))/n, mu, sig, eta, kappa, 
        alf, bet)) * (b - a))/n
}
