"simpPol" <-
function (sumT, par, lambda, lambdac) 
{
    for (i in 1:length(par)) {
        sumT <- sumT + (par[i] * ((lambda - lambdac)/100)^i)
    }
    sumT
}

