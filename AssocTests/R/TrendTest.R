TrendTest <- function(ri, si, varphi = 0.5)
{
    phi <- c(0,varphi,1)
    ni <- ri + si + 1e-10
    r <- sum(ri)
    s1 <- sum(si)
    n <- r + s1

    numerator <- sqrt(n) * sum(phi*(s1*ri-r*si))
    denominator <- sqrt(r*s1*(n*sum(phi^2*ni)-(sum(phi*ni))^2))

    z <- numerator/denominator

    pv = pnorm(-abs(z)) * 2

    list(test.stat=z, p.val=pv)
}
