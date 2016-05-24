`igp` <-
function (t, y, params) 
{
    B <- y[1]
    N <- y[2]
    P <- y[3]
    with(as.list(params), {
        dPdt <- bpb * abp * B * P + bpn * anp * N * P - mp * 
            P
        dNdt <- bnb * abn * B * N - mn * N - anp * N * P
        dBdt <- r * B * (1 - abb * B) - abn * B * N - abp * B * 
            P
        return(list(c(dBdt, dNdt, dPdt)))
    })
}
