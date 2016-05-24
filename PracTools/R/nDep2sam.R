
nDep2sam <- function(S2x, S2y, g, r, rho, alt, del, sig.level=0.05, pow=0.80){
            # check type of test
    alt.ok <- alt %in% c("one.sided", "two.sided")
    if (!alt.ok)
        stop("alt must be either 'one.sided' or 'two.sided'.\n ")
    else {
        if (alt == "one.sided")
            za <- qnorm(1-sig.level)
        if (alt=="two.sided")
            za <- qnorm(1-sig.level/2)
    }
            # check ranges of parameters
    if (g<0 | g>1) stop("g must be in [0,1].\n")
    if (r<0) stop("r must be positive.\n")    
    if (rho<0 | rho>1) stop("rho must be in [0,1].\n")
    if (pow<0 | pow>1) stop("pow must be in [0,1].\n")
    
            # compute sample size in group 1
    zb <- qnorm(1-pow)
    n1 <- (S2x + r*S2y - 2*g*r*rho*sqrt(S2x*S2y)) / del^2 * (za-zb)^2
    
    METHOD <- "Two-sample comparison of means\n Sample size calculation for overlapping samples"

    structure(list(n1 = ceiling(n1),
         n2 = ceiling(n1/r),
         S2x.S2y = c(S2x, S2y),
         delta = del,
         gamma = g,
         r = r,
         rho = rho,
         alt = alt,
         sig.level = sig.level,
         power = pow,
         method = METHOD
    ), class="power.htest")    
}
