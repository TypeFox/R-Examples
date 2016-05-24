pca.covridge <- function( x , ridge = 10^(-10) ){
    cx <- cov(x)
    diag(cx) <- diag(cx) + ridge
    pcax <- stats::princomp( covmat=cx )
    L <- pcax$loadings
    sdev <- pcax$sdev
    D <- diag( pcax$sdev^2)
    scores <- t( t(L) %*%  t( x ) )  # = x %*% L
    res <- list( "loadings" = L , "scores" = scores , 
                "sdev" = sdev )
    return(res)            
        }
