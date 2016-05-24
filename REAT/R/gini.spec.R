gini.spec <-
function (e_ij, e_i) {

e_j <- sum(e_ij)
e <- sum(e_i)
if (e_j > e) { return (NA) }

industries <- length(e_ij)

R_i <- vector()

for (i in 1:industries) { 
R_i[i] <- (e_ij[i]/e_j)/(e_i[i]/e)
}

R_i_sort <- sort(R_i)
lambda <- 1:(nrow(as.data.frame(R_i_sort)))
R_mean <- mean(R_i)
R_i_minus_R_mean <- R_i_sort-R_mean
sum_R <- sum(lambda*R_i_minus_R_mean)
G_j <- (2/((industries^2)*R_mean))*sum_R
return(G_j)

}
