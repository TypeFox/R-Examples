### Contains functions for the Gompertz distribution

pgompertz <- function(q, shape = 1, scale = 1,
                     lower.tail = TRUE, log.p = FALSE,
                      param = c("default", "canonical")){

    n <- length(q)
    if (any(scale == 0)){
        if (log.p) return(rep(-Inf, n))
        else return(rep(0, n))
    }
    ##if ( any(c(shape, scale) <= 0) ){
    if ( any(scale < 0) ){
        cat("scale = ", scale, "\n")
        warning("Negative scale")
        return(NaN)
    }

    if ( any(shape < 0) ){
        warning("Negative shape")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    y <- ifelse(q <= 0, 0, -shape * scale * expm1(q / scale))

    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          -Inf,
                          log(-expm1(y))
                          )
        }else{
            ret <- ifelse(q <= 0,
                          0,
                          y
                          )
        }
    }else{
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          0,
                          -expm1(y) # = 1 - exp(y)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          1,
                          exp(y)
                          )
        }
    }

    return ( ret )
}

dgompertz <- function(x, shape = 1, scale = 1, log = FALSE,
                      param = c("default", "canonical")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    
    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    y <- ifelse(x < 0, 0, x / scale)

    if (log){
        ret <- ifelse(x < 0,
                      0,
                      log(shape) + y - shape * scale * expm1(y))
    }else{
        ret <- ifelse(x < 0,
                      0,
    ##                  shape * exp(y) * exp(-shape * scale * expm1(y)))
                      shape * exp(y - shape * scale * expm1(y)))
    }

    return ( ret )
}

hgompertz <- function(x, shape = 1, scale = 1, log = FALSE,
                      param = c("default", "canonical")){

    if (any(scale == 0)) {
        return(rep(Inf, length(x)))
    }
    ##if ( any(c(shape, scale) <= 0) ){
    if ( any(scale <= 0) ){
        cat("scale = ", scale, "\n")
        warning("Non-positive scale")
        return(NaN)
    }

    if ( any(shape < 0) ){
        cat("shape = ", shape, "\n")
        warning("Negative shape")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    if (log) {
        ret <- ifelse(x < 0,
                      -Inf,
                      log(shape) + x / scale
                      )
    }else{
        ret <- ifelse(x < 0,
                      0,
                      shape * exp(x / scale)
                      )
    }

    return ( ret )
}

qgompertz <- function(p, shape = 1, scale = 1,
                     lower.tail = TRUE, log.p = FALSE,
                      param = c("default", "canonical")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##
    
    if (log.p) p <- exp(p)

    ok <- (p >= 0) & (p <= 1)


    ret <- ifelse(ok,
                  scale * log1p(-log1p(-p) / (shape * scale)),
                  NaN)

    if (!all(ok)) warning("qgompertz produced NaN's")

    return ( ret )
}

Hgompertz <- function(x, shape = 1, scale = 1, log.p = FALSE,
                      param = c("default", "canonical")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    y <- x / scale

    ret <- ifelse(x <= 0,
                  0,
                  shape * scale * expm1(y)
                  )
    if (log.p) ret <- log(ret)

    return ( ret )
}

rgompertz <- function(n, shape = 1, scale = 1,
                      param = c("default", "canonical")){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ## New in 2.1-3:
    param <- param[1]
    if (param == "canonical"){
        ## Transform to "default"
        shape <- shape / scale ## ?
    }else if (param != "default") stop("Illegal 'param'")
    ##

    
    y <- runif(n)

    return ( qgompertz(y, shape, scale, param = "default") )# Note!
}
