`scheffer` <-
function (t, y, p) 
{
    F <- y[1]
    S <- y[2]
    with(as.list(p), {
        n <- N/(1 + qs * S + qf * F)
        dF <- rf * F * (n/(n + hf)) * (1/(1 + af * F)) - lf * 
            F
        dS <- rs * S * (n/(n + hs)) * (1/(1 + as * S + b * F + 
            W)) - ls * S
        return(list(c(dF, dS)))
    })
}
