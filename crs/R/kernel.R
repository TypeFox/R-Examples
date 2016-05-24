## Kernel function used for the categorical variables z

kernel <- function(Z,
                   z,
                   lambda,
                   is.ordered.z=NULL) {

  if(is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda)) stop(" must provide is.ordered.z, Z, z, and lambda")

  if(!is.ordered.z) {
    return(ifelse(Z==z,1,lambda))
  } else {
    return(ifelse(Z==z,1,lambda^abs(Z-z)))
  }

}

## Product kernel function. Z is a vector/matrix, z is a scalar/vector

prod.kernel <- function(Z,
                        z,
                        lambda,
                        is.ordered.z=NULL) {

  Z <- as.matrix(Z)

  if(is.null(is.ordered.z) || missing(Z) || missing(z) || missing(lambda)) stop(" must provide is.ordered.z, Z, z, and lambda")
  if(length(is.ordered.z) != NCOL(Z)) stop(" is.ordered.z and Z incompatible")

  num.z <- NCOL(Z)

  if(num.z != NROW(z) || num.z != NROW(lambda)) stop(paste(" incompatible dimensions for Z, z, and lambda (",num.z,",",NROW(z),",",NROW(lambda),")",sep=""))

  prodker <- kernel(Z=Z[,1],z=z[1],lambda=lambda[1],is.ordered.z=is.ordered.z[1])
  if(num.z > 1) for(i in 2:num.z) prodker <- prodker * kernel(Z=Z[,i],z=z[i],lambda=lambda[i],is.ordered.z=is.ordered.z[i])

  return(prodker)
  
}
