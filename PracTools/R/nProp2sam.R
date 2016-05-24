
nProp2sam <- function(px, py, pxy, g, r, alt, sig.level=0.05, pow=0.80){
                    # check for allowable alt values
    alt.ok <- alt %in% c("one.sided", "two.sided")
    if (!alt.ok)
        stop("alt must be either 'one.sided' or 'two.sided'.\n ")
    else {
        if (alt == "one.sided")
            za <- qnorm(1-sig.level)
        if (alt=="two.sided")
            za <- qnorm(1-sig.level/2)
    }
    
    zb <- qnorm(1-pow)
    del <- px - py
    s2x <- px*(1-px)
    s2y <- py*(1-py)
    sxy <- pxy - px*py
                    # check for allowable pxy
    if ( (pxy <= px*py + sqrt(s2x*s2y)) &
         (pxy >= px*py - sqrt(s2x*s2y)) &
         (pxy <= min(px,py))
       )
         { n1 <- (s2x + r*s2y - 2*g*r*sxy) / del^2 * (za-zb)^2 }
    else stop("pxy not in allowable range.")
    
    METHOD <- "Two-sample comparison of proportions\n Sample size calculation for overlapping samples"

    structure(list(n1 = ceiling(n1),
         n2 = ceiling(n1/r),
         px.py.pxy = c(px, py, pxy),
         gamma = g,
         r = r,
         alt = alt,
         sig.level = sig.level,
         power = pow,
         method = METHOD
    ), class="power.htest")    
}
