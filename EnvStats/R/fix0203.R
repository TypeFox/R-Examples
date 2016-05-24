fix0203 <-
function (func = 1, a = 0, b = 1, mu, sig, eta, kappa, alf = 2, 
    bet = 2, conv = 1e-05, n = 100) 
{
    i0 <- fix0204(func, a, b, n, mu, sig, eta, kappa, alf, bet)
    for (i in (10:20) * n) {
        i1 <- fix0204(func, a, b, i, mu, sig, eta, kappa, alf, 
            bet)
        if (abs((i1 - i0)/i1) < conv) 
            break
        i0 <- i1
    }
    i1
}
