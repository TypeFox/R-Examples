sdfPLS <- function(p, n){
    lams <- 2*pi*seq(from=1/n, to=1/2, by=1/n)
    0.5*(pi^(-p))*p*lams^(p-1)
    }
