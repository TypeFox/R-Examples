gini.conc <-
function (e_ij, e_j) {

e_i <- sum(e_ij)
e <- sum(e_j)
if (e_i > e) { return (NA) }

regions <- length(e_ij)
C_j <- vector()

for (r in 1:regions) { 
C_j[r] <- (e_ij[r]/e_i)/(e_j[r]/e)
}

C_j_sort <- sort(C_j)
lambda <- 1:(nrow(as.data.frame(C_j_sort)))
C_mean <- mean(C_j)
C_j_minus_C_mean <- C_j_sort-C_mean
sum_C <- sum(lambda*C_j_minus_C_mean)
G_i <- (2/((regions^2)*C_mean))*sum_C
return(G_i)

}
