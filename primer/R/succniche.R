`succniche` <-
function (t, y, params) 
{
    S <- y[1]
    E <- y[2]
    M <- y[3]
    R <- y[4]
    F <- max(c(0, 1 - params["D"] - S - E - M - R))
    with(as.list(params), {
        dS = c * (S + R + M) * F - a * c * (M + E) * S - g * 
            S - m * S
        dR = g * (S + M) - m * R
        dM = a * c * (M + E) * S + c * (S + R + M) * E - g * 
            M - m * M
        dE = a * c * (M + E) * F - c * (S + R + M) * E - m * 
            E
        return(list(c(dS, dE, dM, dR)))
    })
}
