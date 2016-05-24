qloglin <- function(u, t){
    if (abs(t) > 1e-6){
        z <- log(1 + ((exp(t) - 1) * u)) / t} else {
        z <- u + t * u * (1 - u) / 2}
    return(z)
}



# implemented before e-Mail of LD of August 2010
# qloglin <- function(u, t){
#     m <- length(u)
#     z <- 1:m * NA
#     II <- (1:m)[abs(t) > 10^(-6)]
#     z[II] <- log( 1 + ((exp(t) - 1) * u[II])) / t
#     II <- (1:m)[abs(t) <= 10^(-6)]
#     z[II] <- u[II] + t * u[II]*(1-u[II])/2
#     return(z)
# }

