maxStep <- function(beta, beta_new, V, eps = 0){

# function t on page 6
prod <- t(V) %*% beta    
prod_new <- t(V) %*% beta_new    
    
ratio <- - prod / (prod_new - prod)
ratio <- ratio[prod_new > eps]

t <- min(1, ratio)
    
return(t)}
