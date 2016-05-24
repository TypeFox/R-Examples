`Q_function_addendo3` <-
function(Sigmaeta, G, n, S11, S00, S10) {
     Q_addendo3 <- n*log(det(Sigmaeta)) + sum(diag( solve(Sigmaeta) %*%  (S11 - S10%*%t(G) - G%*%t(S10) + G%*%S00%*%t(G)) ))
     return(Q_addendo3)
}

