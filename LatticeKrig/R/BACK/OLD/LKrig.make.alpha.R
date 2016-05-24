LKrig.make.alpha<- function( alpha=NA, nu=NULL, nlevel){

    if( is.na(alpha[1]) ) {
        alpha<- rep( NA, nlevel)
    }
    scalar.alpha <- !is.list(alpha) | !is.null(nu)
# determine alpha if nu has been specified
    if (!is.null(nu)) {
        alpha <- exp(-2 * (1:nlevel) * nu)
        alpha <- alpha/sum(alpha)
    }
    if (scalar.alpha & (nlevel != 1) & (length(alpha) == 1)){
                stop( "Only one alpha specifed for multiple levels")}     
# coerce alpha to a list if it is passed as something else
    if (!is.list(alpha)) {
        alpha <- as.list(alpha)
    }
    return( list( alpha=alpha, scalar.alpha= scalar.alpha))
    
  }
