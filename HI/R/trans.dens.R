####
#### Routines for hyperplane inflation simulation technique
#### (see Petris & Tardella)
####

lpi <- log(base::pi)
lsqPi <- log(sqrt(base::pi))

rballunif <- function(n, d) {
    ## generate a point uniformly in the n-dimensional ball
    ## centered at the origin and having radius `d'
    x <- rnorm(n)
    d * runif(1)^(1/n) * x / sqrt(crossprod(x))
}

"trans.dens" <-
function(y, ldens.list, which.models, ..., back.transform=F) {
    ##   y can be a vector
    ##   or a n by p matrix, whose rows are the points at which
    ##   to evaluate trans.dens
    ##   ldens.list is a list of logdensities
    ##   which.models is a sequence of model indices
    ##   For each point of `y', the function returns one of the following: 
    ##   1) if back.transform=F only the _log_ of the transformed density
    ##   2) if back.transform=T a vector of 
    ##      (n+2) elements where 
    ##           *the first one, named "trans.dens",
    ##                is  the _log_ of the transformed density
    ##          **the second one, named "model.index",
    ##                is the corresponding submodel region
    ##         ***the last n elements represent the corresponding 
    ##                x vector in the original submodel space
#     if ( is.null(ldens.list) )
#         stop("ldens.list empty")
#     if ( any(diff(which.models) <= 0) )
#         stop("model indices must be given in increasing order")
    if ( which.models[1] != 0 )
        which.models <- c(0, which.models)
#     if ( length(which.models) != length(ldens.list) )
#         stop("ldens.list and which.models must have the same length")
    if ( !is.matrix(y) ) 
        y <- t(as.matrix(y))
    n.points <- nrow(y)
    n <- ncol(y)
    ans <- matrix( 0, n.points, 2+n )

    ## n=dimension of the largest embedding space
    ## ??? shouldn't it be compatible with (that is the same as)
    ## which.models[length(wihich.models)] ????
    ## add a new check ???

    for ( loop in 1:n.points ) {
        z <- y[loop,]
        i <- length( which.models )
        h <- n - ( k <- which.models[i] )
        if ( h == 0 ){
            logd <- -( ldens.list[[1]](rep(0,n),...) + (n/2)*lpi -
                  (( ldens.list[[i]](...) + lgamma(k/2+1)) ))/k
         }
        else{
            logd <- -( ldens.list[[1]]( c(z[1:h], rep(0,k)), ... ) + (k/2)*lpi -
                  (( ldens.list[[i]](z[1:h],...) + lgamma(k/2+1))) )/k
            }

        ## rewrite the if below in a more efficient way...?
        if ( ( h == 0 && log(sum(z^2)) > 2*logd ) || ( h > 0 && log(sum(z[-(1:h)]^2)) > 2*logd ) ) {
            if ( h == 0 ){
                z <- ldeflate( z, logd )
            }
            else{
                z[-(1:h)] <- ldeflate( z[-(1:h)], logd )
                }
            i <- i - 1           
            h <- n - ( k <- which.models[i] )
            if ( i > 1 ) {
                logd <- -(( ldens.list[[1]]( c(z[1:h], rep(0,k)), ... )+(k/2)*lpi) -
                          (( ldens.list[[i]](z[1:h],...) + lgamma(k/2+1)) ))/k
                while ( (i > 1) && (log(sum(z[-(1:h)]^2)) > (2*logd) )) {
                    z[-(1:h)] <- ldeflate( z[-(1:h)], logd )
                    i <- i - 1
                    h <- n - ( k <- which.models[i] )
                    if ( k > 0 ) {
                        logd <- -( ldens.list[[1]]( c(z[1:h], rep(0,k)), ... ) + (k/2)*lpi -
                              (( ldens.list[[i]](z[1:h],...) + lgamma(k/2+1)) ))/k
       c(z[1:h], rep(0,k))             }
                }
            }
        }
        if ( k == n )
            x <- rep(0,k)
        else
            x <- c( z[1:h], rep(0,k) )
        ans[loop,] <- c(ldens.list[[1]](x,...), k, x )
    }
    if ( back.transform ) { 
        dimnames(ans) <- list(NULL,c("trans.dens", "model.index", paste("x",1:n,sep=".")))
        return(ans)
    }
    else
        return(ans[,1])
}



linflate <- function(y, logh) {
    norm.y <- sqrt( y %*% y )
    n <- length(y)
    y * ( norm.y^n + exp(n*logh) )^(1/n) / norm.y
}
   

ldeflate <- function(y, logh) {
    norm.y <- sqrt( y %*% y )
    n <- length(y)
    y * ( norm.y^n - exp(n*logh) )^(1/n) / norm.y
}

"trans.up" <-
function(x, ldens.list, which.models, ...) {
    ## aim: it maps back a point in the "original model scale"
    ##      to an appropriate point in the "inflated scale"
    ##      corresponding to the same submodel subspace

    ## x              can be a vector or a n by p matrix, 
    ##                whose rows are the points to be transformed 
    ## ldens.list      is a list of densities
    ## which.models   is a sequence of model indices
#     if ( is.null(ldens.list) )
#         stop("ldens.list empty")
#     if ( any(diff(which.models) <= 0) )
#         stop("model indices must be given in increasing order")
    if ( which.models[1] != 0 )
        which.models <- c(0, which.models)
#     if ( length(which.models) != length(ldens.list) )
#         stop("ldens.list and which.models must have the same length")
    if ( !is.matrix(x) ) 
        x <- t(as.matrix(x))
    n.points <- nrow(x)
    n <- ncol(x)
    ans <- matrix( 0, n.points, n )
    for ( loop in 1:n.points ) {
        z <- x[loop,]
        k <- 0
        while ( (z[n-k] == 0) && (k < n) ) k <- k+1
        if ( k > 0 ) # some components are zero
            if ( length(j <- which( which.models == k )) > 0 ) {
                m <- n-k
                if ( m == 0 ) {
                    logd <- ( ldens.list[[j]](...) - ldens.list[[1]]( rep(0,k), ... ) +
                           lgamma(k/2+1) ) / k  -  0.5*lpi
                    z <- rballunif( k, exp(logd) )
                }
                else {
                    logd <- ( ldens.list[[j]](z[1:m],...) -
                           ldens.list[[1]]( c(z[1:m], rep(0,k)), ... ) +
                           lgamma(k/2+1) ) / k  -  0.5*lpi
                    z[-(1:m)] <- rballunif( k, exp(logd) )
                }
                
            }
        if ( k < n ) { # not all the components are zero
            for ( i in (k+1):n )
                if ( length(j <- which( which.models == i )) > 0 ) {
                    m <- n-i
                    if ( m == 0 ) {
                        logd <- ( ldens.list[[j]](...) - ldens.list[[1]]( rep(0,i), ... ) +
                               lgamma(i/2+1) ) / i  -  0.5*lpi
                        z <- linflate( z, logd )
                    }
                    else {
                        logd <- ( ldens.list[[j]](z[1:m],...) -
                               ldens.list[[1]]( c(z[1:m], rep(0,i)), ... ) +
                               lgamma(i/2+1) ) / i  -  0.5*lpi
                        z[-(1:m)] <- linflate( z[-(1:m)], logd )
                    }
                }
        }
        ans[loop,] <- z
    }
    return(ans)
}

### transformed log density for the mixture of two components
"trans2" <-
function(y, ldens.list, k, ...) {
    ## 'y' is a vector
    ## returns also the model index of the vector 'y'
    if (length(ldens.list) != 2 || length(k) != 1)
        stop("ldens.list must have length 2\nand k must have length 1")
    n <- length(y)
    h <- n - k # dimension of submodel 'k'
    if (h==0)
    {
        ldk <- ldens.list[[2]](...)
        if ( is.infinite(ldk) && ldk < 0 )
            return(c(ldens.list[[1]](y,...), 0))
        ld0 <- ldens.list[[1]](rep(0,n),...)
    }
    else
    {
        ldk <- ldens.list[[2]](y[1:h],...)  
        if ( is.infinite(ldk) && ldk < 0 )
            return(c(ldens.list[[1]](y,...), 0))
        ld0 <- ldens.list[[1]](c(y[1:h],rep(0,k)),...)
    }
    if ( is.infinite(ld0) && ld0 < 0 )
        stop(paste("ldens.list[[1]] may not take the value",ld0))
    u <-  ldk - ld0 + lgamma(k/2 + 1) - k*lsqPi -
        0.5 * k * log(crossprod(y[(n-k+1):n]))
    if ( u > 0 )
        return(c(ld0, k))
    else 
        return(c(ldens.list[[1]](y * rep(c(1,(1-exp(u))^(1/k)),c(h,k)),...), 0))
}

### Map points in the auxiliary space back to the original one
"transBack2" <-
function(y, ldens.list, k, ...) {
    ## 'y' is a vector or a matrix
    ## back.transform not implemented yet
    ## returns also the model index of the vector 'y'
    if (length(ldens.list) != 2 || length(k) != 1)
        stop("ldens.list must have length 2\nand k must have length 1")
    if ( !is.null(dim(y)) )
    {
        ans <- matrix(0,NROW(y),NCOL(y))
        for (i in 1:NROW(ans))
            ans[i,] <- Recall(y[i,], ldens.list, k, ...)
        return(ans)
    }
    n <- length(y)
    h <- n - k # dimension of submodel 'k' 
    if (h==0)
    {
        ldk <- ldens.list[[2]](...)
        if ( is.infinite(ldk) && ldk < 0 )
            return(y)
        z <- rep(0,n)
    }
    else
    {
        ldk <- ldens.list[[2]](y[1:h],...)  
        if ( is.infinite(ldk) && ldk < 0 )
            return(y)
        z <- c(y[1:h],rep(0,k))
    }
    ld0 <- ldens.list[[1]](z, ...)
    if ( is.infinite(ld0) && ld0 < 0 )
        stop(paste("ldens.list[[1]] may not take the value",ld0))
    u <-  ldk - ld0 + lgamma(k/2 + 1) - k*lsqPi -
        0.5 * k * log(crossprod(y[(n-k+1):n]))
    if ( u > 0 )
        return(z)
    else 
        return(y * rep(c(1,(1-exp(u))^(1/k)),c(h,k)))
}

"transUp2" <-
function(y, ldens.list, k, ...) {
    ## 'y' is a vector
    n <- length(y)
    h <- n - k # dimension of submodel 'k' 
    ind.h <- seq(1,length.out=h)
    if (h==0)
    {
        ldk <- ldens.list[[2]](...)
        if ( is.infinite(ldk) && ldk < 0 ) return(y)
        ld0 <- ldens.list[[1]](rep(0,n),...)
    }
    else
    {
        ldk <- ldens.list[[2]](y[ind.h],...)  
        if ( is.infinite(ldk) && ldk < 0 ) return(y)
        ld0 <- ldens.list[[1]](c(y[ind.h],rep(0,k)),...)
    }
    if ( is.infinite(ld0) && ld0 < 0 )
        stop(paste("ldens.list[[1]] may not take the value",ld0))
    r <-  ldk - ld0 + lgamma(k/2 + 1) - k*lsqPi 
    u <- r - 0.5 * k * log(crossprod(y[(n-k+1):n]))
    if ( u > 0 )
    {
        ## model 'k'
        return( c(y[ind.h], rballunif(k,exp(r/k))) )
    }
    else
    {
        ## model '0'
        return( y * rep(c(1,(1-exp(u))^(1/k)),c(h,k)) )
        
    }
}    

