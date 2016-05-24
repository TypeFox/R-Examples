multrnks <-
function (rk, m1 = 2, m2 = 2, m = 2) 
    {
        n <- length(rk)
        pk <- rk/n
        q <- rep(0, n)
        for (l in m1:m2) {
            q <- q + (l * choose(m, l) * (pk^(l - 1)) * ((1 - pk)^(m - 
                                                                       l)))
        }
        q
    }
