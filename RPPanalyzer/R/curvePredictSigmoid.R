`curvePredictSigmoid` <-
function(x,params)  {

    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]

    v <- alpha + beta*(2^(x*gamma))/(1+2^(x*gamma))

    return(v)

}

