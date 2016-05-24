## initialize method 
setMethod("initialize", "Gumbel",
    function(.Object, loc = 0, scale = 1) {
        .Object@img <- Reals()
        .Object@param <- new("GumbelParameter", loc = loc, scale = scale, 
                             name = gettext("parameter of a Gumbel distribution"))
        .Object@r <- function(n){}
        body(.Object@r) <- substitute({ rgumbel(n, loc = loc1, scale = scale1) },
                                     list(loc1 = loc, scale1 = scale))
        .Object@d <- function(x, log = FALSE){}
        body(.Object@d) <- substitute({  dgumbel(x, loc = loc1, scale = scale1, log = log) },
                                     list(loc1 = loc, scale1 = scale))
        .Object@p <- function(q, lower.tail = TRUE, log.p = FALSE){}
        body(.Object@p) <- substitute({p1 <- pgumbel(q, loc = loc1, scale = scale1, lower.tail = lower.tail) 
                                       return(if(log.p) log(p1) else p1)},
                                     list(loc1 = loc, scale1 = scale))
        .Object@q <- function(p, loc = loc1, scale = scale1, lower.tail = TRUE, log.p = FALSE){}
            body(.Object@q) <- substitute({
                        ## P.R.: changed to vectorized form 
                        p1 <- if(log.p) exp(p) else p
                                                                        
                        in01 <- (p1>1 | p1<0)
                        i01 <- distr:::.isEqual01(p1)
                        i0 <- (i01 & p1<1)   
                        i1 <- (i01 & p1>0)
                        ii01 <- distr:::.isEqual01(p1) | in01
                                      
                        p0 <- p
                        p0[ii01] <- if(log.p) log(0.5) else 0.5
                                      
                        q1 <- qgumbel(p0, loc = loc1, scale = scale1, 
                                      lower.tail = lower.tail) 
                        q1[i0] <- if(lower.tail) -Inf else Inf
                        q1[i1] <- if(!lower.tail) -Inf else Inf
                        q1[in01] <- NaN
                        
                        return(q1)  
                     },  list(loc1 = loc, scale1 = scale))
        .Object@.withSim   <- FALSE
        .Object@.withArith <- FALSE
        .Object@.logExact <- FALSE
        .Object@.lowerExact <- TRUE
        .Object
    })

