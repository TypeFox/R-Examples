krugman.spec <-
function (e_ij, e_il, e_j, e_l) {

industries <- length(e_ij)
s_ij <- vector()
s_il <- vector()

for (i in 1:industries) {
s_ij[i] <- e_ij[i]/e_j
s_il[i] <- e_il[i]/e_l
}

K_jl <- sum(abs(s_ij-s_il))
return(K_jl)
}
