coxDeriv <- function(beta, ttf, tf, Z){

p <- dim(Z)[2]
prod <- Z %*% beta
alpha <- exp(prod)
G1 <- apply(matrix(rbind(Z[tf == 1, ], rep(0, p)), ncol = p), 2, sum)
G2 <- 0

event <- ttf[tf == 1]

# log-likelihood
for (j in 1 : sum(tf)){
     risk <- ttf >= event[j]
     risk <- (1:length(tf)) * risk
     risk <- risk[risk > 0]
     
     G3 <- apply(rbind(matrix(apply(Z, 2,  function(vec, alpha){vec * alpha}, alpha)[risk, ], ncol = p), rep(0, length(beta))), 2, sum) / sum(alpha[risk])
     G2 <- G2 + G3
}

G <- G1 - G2

return(list("dL" = G))
}


