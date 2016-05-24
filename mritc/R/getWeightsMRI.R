getWeightsMRI <- function(nnei, sigma=1){
    if(! nnei %in% c(26))
        stop("'nnei' has to be among (26)")

    width <- ((nnei + 1)^(1/3) - 1) /  2
    expand <- t(data.matrix(expand.grid(-width:width,
                                     -width:width,
                                     -width:width)[-(nnei/2 +1),]))
 
    weights <- dnorm(sqrt((expand[1,]/width)^2 + (expand[2,]/width)^2 +
                          (expand[3,]/width)^2), mean=0, sd=sigma)
    weights /sum(weights)
}
