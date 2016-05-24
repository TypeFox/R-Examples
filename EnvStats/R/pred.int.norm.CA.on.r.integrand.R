pred.int.norm.CA.on.r.integrand <-
function (v.weird, n.weird, df.weird = n.weird - 1, n.mean = 1, 
    K.weird, delta.over.sigma, m.weird, r.weird) 
{
    term.1 <- pT(q = sqrt(n.weird) * K.weird, df = df.weird, 
        ncp = sqrt(n.weird/n.mean) * (qnorm(v.weird) + sqrt(n.mean) * 
            delta.over.sigma))
    a <- v.weird^(m.weird - 2)
    term.2 <- r.weird * (v.weird * (1 + a * (1 - v.weird)))^(r.weird - 
        1)
    term.3 <- 1 + a * (m.weird - 1 - m.weird * v.weird)
    term.1 * term.2 * term.3
}
