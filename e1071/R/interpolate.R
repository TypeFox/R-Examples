interpolate <- function(x, a, adims=lapply(dimnames(a), as.numeric),
                        method="linear"){

    if(is.vector(x)) x<- matrix(x, ncol=length(x))

    if(!is.array(a))
        stop("a is not an array")
    
    ad <- length(dim(a))
    
    method <- pmatch(method, c("linear", "constant"))
    if (is.na(method)) 
        stop("invalid interpolation method")
    
    if(any(unlist(lapply(adims, diff))<0))
        stop("dimensions of a not ordered")
    
    retval <- rep(0, nrow(x))
    bincombi <- bincombinations(ad)

    convexcoeff <- function(x, y) {
        ok <- y>0
        x[ok] <- y[ok]-x[ok]
        x
    }

    for(n in 1:nrow(x)){
        
        ## the "leftmost" corner of the enclosing hypercube 
        leftidx <- rep(0, ad)
        xabstand <- rep(0, ad)
        aabstand <- rep(0, ad)
        
        for(k in 1:ad){
            if(x[n,k] < min(adims[[k]]) || x[n,k] > max(adims[[k]]))
                stop("No extrapolation allowed")
            else{
                leftidx[k] <- max(seq(adims[[k]])[adims[[k]] <= x[n,k]])
                ## if at the right border, go one step to the left
                if(leftidx[k] == length(adims[[k]]))
                    leftidx[k] <- leftidx[k] - 1
                
                xabstand[k] <- x[n,k] - adims[[k]][leftidx[k]]
                aabstand[k] <- adims[[k]][leftidx[k]+1] - 
                    adims[[k]][leftidx[k]]
            }
        }

        coefs <- list()
        if(method==1){
            for(k in 1:(2^ad)){
                retval[n] <- retval[n] +
                    element(a, leftidx+bincombi[k,]) *
                        prod((aabstand-
                              convexcoeff(xabstand,
                                          aabstand*bincombi[k,]))/aabstand)
            }
        }
        else if(method==2){
            retval[n] <- element(a, leftidx)
        }
    }
    
    names(retval) <- rownames(x)
    retval
}

