IyenGreenWeight <- function(x, b, q, type = 1, alpha = 0.05){

    ## Weight function w_2 in Iyengar & Greenhouse (1988, p. 112)
    ## note that the parametrization here is that x lives on the 
    ## t-statistic scale (not on p-value scale)

    t <- qt(1 - alpha / 2, df = q)
    ind <- abs(x) <= t
    res <- rep(1, length(x))

    if (type == 1){res[ind] <- (abs(x[ind]) / t) ^ b}
    if (type == 2){res[ind] <- rep(exp(- b), sum(ind))}
    return(res)
    }
