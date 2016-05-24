#### Tests for psi(), rho() etc

library(robust)

### Test internal consistency
## (1) rho'(x) = psi(x)
## (2) psi'(x) = psp(x)

## (3) chi(x) = pmin(R, rho(x))  ( *not* the definition for S-estimates !! )
##             but really chi() = 1/2 psi()^2  Huber; for the others chi() == rho()

##
x  <- seq(-4,4, length =  801)# large n
x. <- seq(-4,4, length = 2001)# larger n

verbose <- TRUE
for(ipsi in 1:3) { # later  1:4 --  ipsi = 4 nowhere documented
    cat("ipsi = ", ipsi, "\n")
    f.x <- cbind(psi = psi.weight(x, ips=ipsi), psp = psp.weight(x, ips=ipsi),
                 chi = chi.weight(x, ips=ipsi), rho = rho.weight(x, ips=ipsi))
    rhoF <- splinefun(x,f.x[,"rho"])
    psiF <- splinefun(x,f.x[,"psi"])
    pspF <- splinefun(x,f.x[,"psp"])
##     chiF <- splinefun(x,f.x[,"chi"])
    p1 <- psiF(x., deriv=1); p. <- pspF(x.)
    r1 <- rhoF(x., deriv=1); ps <- psiF(x.)
    if(verbose) {
        cat("psi'(.) = psp(.):", all.equal(p1, p., tol = 1e-6),"\n")
        cat("rho'(.) = psi(.):", all.equal(r1, ps, tol = 1e-6),"\n")
        ## TODO: chi ?
    }
    stopifnot(all.equal(p1, p., tol = if(ipsi == 3) .05 else 1e-3),
              all.equal(r1, ps, tol = 1e-4))
    if(verbose) cat("\n---------------------\n\n")
}


### Plots --> ../man/weight.funs.Rd
