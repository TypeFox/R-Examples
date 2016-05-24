##
## Function for creating 'SOCC' objects 
socc <- function(F, g, d, f){
    G <- matrix(0, nrow = nrow(F) + 1, ncol = ncol(F))
    G[1, ] <- -d
    G[-1, ] <- -F
    h <- matrix(c(f, g), ncol = 1)   
    return(list(conType = "SOCC", G = G, h = h, dims = nrow(G)))    
}
