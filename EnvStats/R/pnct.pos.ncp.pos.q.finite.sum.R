pnct.pos.ncp.pos.q.finite.sum <-
function (q, df, ncp, start, end) 
{
    j <- start:end
    ncp.sq.o.2 <- (ncp^2)/2
    jo2 <- j/2
    q.sq <- q^2
    summand.vec <- exp((-ncp.sq.o.2) + jo2 * log(ncp.sq.o.2) - 
        lgamma(jo2 + 1) + log(pbeta(q.sq/(df + q.sq), jo2 + 0.5, 
        df/2)))/2
    sum(summand.vec)
}
