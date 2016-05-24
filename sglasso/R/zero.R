zero <- function(theta, truncate){
    id <- abs(theta) <= truncate
    theta[id] <- 0
    theta
}
