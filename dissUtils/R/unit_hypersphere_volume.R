unit.hypersphere.volume <- function(D){
#D is the number of dimensions of the hypersphere
    return(exp((D / 2) * log(pi) - log(gamma(1 + D/2))));
}
