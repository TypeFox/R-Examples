h_conditional <- function(t, z, beta, lambda, a){
     result <- lambda * a * (t ^ (a - 1)) * (exp(z %*% beta))
     return(result)
}









#