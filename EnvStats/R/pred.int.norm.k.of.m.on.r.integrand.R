pred.int.norm.k.of.m.on.r.integrand <-
function (v.weird, n.weird, df.weird = n.weird - 1, n.mean = 1, 
    K.weird, delta.over.sigma, k.weird, m.weird, r.weird) 
{
    term.1 <- pT(q = sqrt(n.weird) * K.weird, df = df.weird, 
        ncp = sqrt(n.weird/n.mean) * (qnorm(v.weird) + sqrt(n.mean) * 
            delta.over.sigma))
    term.2 <- r.weird * pbeta(q = v.weird, shape1 = k.weird, 
        shape2 = m.weird + 1 - k.weird)^(r.weird - 1)
    term.3 <- (v.weird^(k.weird - 1) * (1 - v.weird)^(m.weird - 
        k.weird))/beta(k.weird, m.weird + 1 - k.weird)
    term.1 * term.2 * term.3
}
