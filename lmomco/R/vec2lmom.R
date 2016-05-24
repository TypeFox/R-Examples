"vec2lmom" <-
function(vec, lscale=TRUE, trim=NULL, leftrim=NULL, rightrim=NULL, checklmom=TRUE) {
    z <- list(lambdas=vector(mode="numeric", length=length(vec)),
              ratios=vector(mode="numeric", length=length(vec)),
              trim=trim,
              leftrim=leftrim,
              rightrim=rightrim,
              source="vec2lmom"
             )
    n <- length(vec)
    if(n < 2) {
        warning("function expects a minimum of first two L-moments")
        return(NA)
    }
    z$lambdas[1] <- vec[1]   # the mean
    z$lambdas[2] <- ifelse(lscale == TRUE, vec[2], vec[2]*z$lambdas[1])
    z$ratios[1]  <- NA
    z$ratios[2]  <- z$lambdas[2]/z$lambdas[1]

    if(n >= 3) {
       z$ratios[3:n]  <- vec[3:n] # ratios mandated
       z$lambdas[3:n] <- z$ratios[3:n]*z$lambdas[2]
       if(checklmom) {
          if(! are.lmom.valid(z)) {
             warning("L-moments are invalid, but still returning the values")
          }
       }
    }

    return(z)
}


