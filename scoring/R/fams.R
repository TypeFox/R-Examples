betafam <-
function(p, d, param, ...){
    ## Define integrand:
    ## Inf and NaN arise when 0 is raised to a negative power.
    ## This happens when d==p.

    ## Convert to vector form
    p <- p[2]
    d <- 1*d==2
    
    "betaint" <- function(p, d, param){
      ifelse(d==p, 0,
             (d*(1-p) + (1-d)*p)*(p^(param[1]-1))*((1-p)^(param[2]-1)))
    }

    lbnd <- d*p
    ubnd <- d + (1-d)*p

    res <- try(integrate(betaint, lbnd, ubnd, d=d, param=param, subdivisions=100)$value, silent=TRUE)

    if(inherits(res,"try-error")){
        warning("A score evaluates to infinity.")
        res <- Inf
    }

    res
}


powfam <-
function(p, d, param, ...){
    ## Obtain score from power family with baseline for a single
    ## forecast.
    gam <- param[1]
    if(length(param)==1){
      ## We are using the family without baselines:
      cpar <- rep(1, length(p))
    } else {
      ## We are using the family with baselines
      cpar <- param[2:length(param)]
    }

    ## JNW 2009, eq (11), numerator of first term
    num1 <- (p[d]/cpar[d])^(gam-1) - 1
    ## JNW 2009, eq (11), numerator of second term
    num2 <- sum(p^gam/cpar^(gam-1)) - 1

    ## Now combine them
    -(num1/(gam-1) - num2/gam)
}


sphfam <-
function(p, d, param, ...){
    ## Obtain score from pseudospherical family with baseline for a single
    ## forecast.
    gam <- param[1]
    if(length(param)==1){
      ## We are using the family without baselines:
      cpar <- rep(1, length(p))
    } else {
      ## We are using the family with baselines
      cpar <- param[2:length(param)]
    }

    ## JNW 2009, eq (12), numerator and denominator of term in square brackets
    num <- p[d]/cpar[d]
    denom <- (sum(p^gam/cpar^(gam-1)))^(1/gam)

    ## Now combine them
    -(1/(gam-1)) * ((num/denom)^(gam-1) - 1)
}


## NB Beta family param is just for a 2-alternative rule,
##    while pow and sph are for multi-alternatives.
ordwrap <- function(data, param, fam){
    ## Wrapper to get ordinal score out of any other family
    ## (See Jose, Nau, Winkler, 2009, Management Science, Eq (6) + (13))
    nalts <- ncol(data) - 1
    p <- data[,1:nalts]
    p1 <- t(apply(p, 1, cumsum))
    p1 <- p1[,1:(ncol(p1) - 1)]
    p2 <- 1 - p1
    d <- data[,(nalts + 1)]
    dmat <- matrix(2, nrow(p1), ncol(p1))
    tmpd <- cbind(d, dmat)
    dmat <- t(apply(tmpd, 1, function(x){
              ## Handle case where outcome is last
              ## alternative (so R doesn't add to
              ##              the vector)
              if(x[1] < length(x)){
                  x[(x[1]+1):length(x)] <- 1
              }
              x[2:length(x)]}))

    ## Cumulative baselines for pow and sph
    if(length(param) > 1){
        Q <- cumsum(param[2:length(param)])
    }

    ## FIXME: Modify calcscore to accept different param values
    ##        for each forecast row, then this loop can be removed.
    ## TODO: This assumes a single family parameter (ok for pow and sph),
    ##       followed by baselines.  Greater flexibility will be allowed
    ##       if baselines are separate argument from family parameters.
    runscore <- rep(0, nrow(data))
    for(i in 1:ncol(p1)){
        if(length(param) > 1){
            tmpparam <- c(param[1], Q[i], 1 - Q[i])
        } else {
            tmpparam <- param
        }
        tmpdat <- cbind(p1[,i], p2[,i], dmat[,i])

        tmpscore <- scoreitems(tmpparam, tmpdat, fam = fam,
                               ordered = FALSE)

        runscore <- runscore + tmpscore
    }

    runscore
}


logfam <- function(p, d, ...){
    -log(p[d])
}
