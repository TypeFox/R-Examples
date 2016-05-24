coxLoglik <- function(beta, ttf, tf, Z){

pro <- Z %*% beta
S1 <- sum(tf * pro)
S2 <- 0

event <- ttf[tf == 1]

for (j in 1:sum(tf)){
     risk <- ttf >= event[j]
     risk <- (1:length(tf)) * risk
     risk <- risk[risk > 0]
     S2 <- S2 + log(sum(exp(pro[risk])))
}

L <- S1 - S2
return(list("L" = L))
}

