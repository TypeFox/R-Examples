`lanczos` <- function(z, log=FALSE){
    stopifnot(all(Re(z)>=0.5))

    g <- 7
    p <- c(
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
        -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
        -176.61502916214059906584551354,
        12.507343278686904814458936853,
        -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
        )

     
    z <- z+0i-1  #NB coerces to complex
    
    x <- p[1]
    for(i in seq.int(from=2,to=g+2)){ x <- x + p[i]/(z+i-1) }
    tee <- z + g + 0.5

    if(log){
        return(log(2*pi)/2 +(z+0.5)*log(tee) -tee +log(x))
    } else {
        return(sqrt(2*pi) *tee^(z+0.5) *exp(-tee) *x)
    }
}

`complex_gamma` <- function(z,log=FALSE){
    out <- z*NaN +0i
    z <- z+0i;
    left <- Re(z)<0.5
    zl <- z[ left]
    zr <- z[!left]

    if(log){        
        out[ left] <- log(pi)-log(sin(pi*zl))-lanczos(1-zl,log=TRUE)
        out[!left] <- lanczos(zr,log=TRUE)
    } else {
        out[ left] <- pi/(sin(pi*zl)*lanczos(1-zl,log=FALSE))
        out[!left] <- lanczos(zr,log=FALSE)
    }
    return(out)
}

`complex_factorial` <- function(z,log=FALSE){complex_gamma(z+1,log=log)}
