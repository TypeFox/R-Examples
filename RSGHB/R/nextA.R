nextA <- function(b,d,env)
{
     SQNP <- sqrt(env$gNP)
     
     # draws from multivariate normals can be taken using 
     # draws = mu + l * eta where l*l' = d (l is the cholesky decomposition)      
     
     #accounting for non-diffuse priors
     return(colMeans(b) + t(chol(d))%*%matrix(rnorm(env$gNIV),nrow=env$gNIV,ncol=1)/SQNP)
}
