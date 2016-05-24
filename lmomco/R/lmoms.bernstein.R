"lmoms.bernstein" <-
function(x, bern.control=NULL,
            poly.type=c("Bernstein", "Kantorovich", "Cheng"),
            bound.type=c("none", "sd", "Carv", "either"),
            fix.lower=NULL, fix.upper=NULL, p=0.05) {

    if(! is.null(bern.control)) {
         poly.type  <- bern.control$poly.type
         bound.type <- bern.control$bound.type
         fix.lower  <- bern.control$fix.lower
         fix.upper  <- bern.control$fix.upper
         p          <- bern.control$p
    }

    n <- length(x)
    if (length(unique(x)) == 1)
        stop("all values are equal--Lmoments can not be computed")
    "afunc" <- function(f,PWMf=function(f) { 1 }, ...) {
                          dat2bernqua(f, ...)*PWMf(f)  }
    L1 <- integrate(afunc, 0,1, x=x, poly.type=poly.type,
             bound.type=bound.type, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
    L2 <- integrate(afunc, 0,1, x=x, poly.type=poly.type,
             PWMf=function(f) { (2*f - 1) },
             bound.type=bound.type, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
    L3 <- integrate(afunc, 0,1, x=x, poly.type=poly.type,
             PWMf=function(f) { (6*f^2 - 6*f + 1) },
             bound.type=bound.type, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
    L4 <- integrate(afunc, 0,1, x=x, poly.type=poly.type,
             PWMf=function(f) { (20*f^3 - 30*f^2 + 12*f - 1) },
             bound.type=bound.type, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
    L5 <- integrate(afunc, 0,1, x=x, poly.type=poly.type,
             PWMf=function(f) { (70*f^4 - 140*f^3 + 90*f^2 - 20*f + 1) },
             bound.type=bound.type, p=p, fix.lower=fix.lower, fix.upper=fix.upper)
    lambdas <- c(L1$value, L2$value, L3$value, L4$value, L5$value)
    ratios  <- rep(NA, length(lambdas))
    ratios[2] <- lambdas[2]/lambdas[1]
    for(i in 3:length(lambdas)) { ratios[i] <- lambdas[i]/lambdas[2] }
    z <- list(lambdas=lambdas, ratios=ratios)
    return(z)
}

