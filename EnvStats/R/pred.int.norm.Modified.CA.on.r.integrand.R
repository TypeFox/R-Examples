pred.int.norm.Modified.CA.on.r.integrand <-
function (v.weird, n.weird, df.weird = n.weird - 1, n.mean = 1, 
    K.weird, delta.over.sigma, r.weird) 
{
    term.1 <- pT(q = sqrt(n.weird) * K.weird, df = df.weird, 
        ncp = sqrt(n.weird/n.mean) * (qnorm(v.weird) + sqrt(n.mean) * 
            delta.over.sigma))
    term.2 <- r.weird * (v.weird * (1 + v.weird * (3 - v.weird * 
        (5 - 2 * v.weird))))^(r.weird - 1)
    term.3 <- 1 + v.weird * (6 - v.weird * (15 - 8 * v.weird))
    term.1 * term.2 * term.3
}
