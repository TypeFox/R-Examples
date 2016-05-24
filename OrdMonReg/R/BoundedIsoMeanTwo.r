BoundedIsoMeanTwo <- function(g1, w1, g2, w2, K1 = 1000, K2 = 400, delta = 10 ^ (-4), errorPrec = 10, output = TRUE){

tmp <- BoundedAntiMeanTwo(-g2, w2, -g1, w1, K1, K2, delta, errorPrec, output)

res <- list("g1" = -tmp$g2, "g2" = -tmp$g1, "L" = tmp$L, "error" = tmp$error, "k" = tmp$k, "tau" = tmp$t_old)
return(res)
}

