################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Gaussian spatial interaction function for twinstim's epidemic component
###
### Copyright (C) 2009-2014 Sebastian Meyer
### $Revision: 771 $
### $Date: 2014-02-17 10:52:54 +0100 (Mon, 17 Feb 2014) $
################################################################################


## nTypes: determines the number of parameters=(log-)standard deviations of the
##   Gaussian kernel. In a multitype epidemic, the different types may share the
##   same spatial interaction function (type-invariant), in which case nTypes=1.
##   Otherwise nTypes should equal the number of event types of the epidemic, in
##   which case every type has its own variance parameter.
## logsd: logical indicating if the gaussian kernel should be reparametrized
##   such that the log-standard deviation is the parameter in question. This
##   avoids constrained optimisation (L-BFGS-B) or the use of 'validpars'.
## density: logical. If TRUE, the isotropic Gaussian density (on R^2) will not
##   be scaled to have maximum value of 1 at the mean c(0,0).
## effRangeMult: determines the effective range for numerical integration in
##   terms of multiples of the parameter, i.e. with effRangeMult=6 numerical
##   integration only considers the 6-sigma area around the event instead of the
##   whole observation region W.
## validpars: If logsd = FALSE, you should either use
##   constrained optimisation (L-BFGS-B) or set 'validpars' to function (pars)
##   pars > 0. 

siaf.gaussian <- function (nTypes = 1, logsd = TRUE, density = FALSE,
                           F.adaptive = TRUE, effRangeMult = 6, validpars = NULL)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)
    if (isScalar(F.adaptive)) {
        adapt <- F.adaptive
        F.adaptive <- TRUE
    } else adapt <- 0.1

    f <- function (s, pars, types) {}       # coordinate matrix s, length(types) = 1 or nrow(s)
    F <- if (F.adaptive) {
        as.function(c(alist(polydomain=, f=, pars=, type=),
                      list(adapt=adapt), quote({})))
    } else siaf.fallback.F
    Fcircle <- function (r, pars, type) {}  # single radius and type
    effRange <- function (pars) {}
    deriv <- function (s, pars, types) {}   # coordinate matrix s, length(types) = 1 or nrow(s)
    Deriv <- function (polydomain, deriv, pars, type, nGQ = 20L) {} # single "owin" and type
    simulate <- function (n, pars, type, ub) {} # n=size of the sample,
                                                # type=single type,
                                                # ub=upperbound (unused here)

    ## if there is only one type, we set the default type(s) argument to 1
    ## (it is actually unused inside the functions)
    if (nTypes == 1L) {
        formals(f)$types <- formals(Fcircle)$type <- formals(deriv)$types <-
            formals(Deriv)$type <- formals(simulate)$type <- 1L
    }
    
    # helper expressions
    tmp1 <- if (logsd) expression(sds <- exp(pars)) else expression(sds <- pars)
    tmp1.1 <- if (nTypes==1L) expression(sd <- sds) else expression(sd <- sds[type])
    tmp2 <- c(
        expression(sLengthSquared <- .rowSums(s^2, L <- nrow(s), 2L)),
        if (nTypes == 1L) expression(sdss <- sds) else expression(
            types <- rep_len(types, L),
            sdss <- sds[types]
            )
        )

    # spatial interaction function
    body(f) <- as.call(c(as.name("{"),
        tmp1, tmp2,
        expression(fvals <- exp(-sLengthSquared/2/sdss^2)),
        if (density) expression(fvals / (2*pi*sdss^2)) else expression(fvals)
    ))

    # numerically integrate f over a polygonal domain
    if (F.adaptive) {
        body(F) <- as.call(c(as.name("{"),
            tmp1, tmp1.1,
            expression(
                eps <- adapt * sd,
                intf <- polyCub.midpoint(polydomain, f, pars, type, eps=eps),
                intf
                )
        ))
    }
    
    # calculate the integral of f over a circular domain around 0
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp1, tmp1.1,
        expression(val <- pchisq((r/sd)^2, 2)), # cf. Abramowitz&Stegun formula 26.3.24
        if (!density) expression(val <- val * 2*pi*sd^2),
        expression(val)
    ))

    # effective integration range of f as a function of sd
    if (isScalar(effRangeMult)) {
        body(effRange) <- as.call(c(as.name("{"),
            tmp1,
            substitute(effRangeMult*sds)
        ))
    } else effRange <- NULL

    # derivative of f wrt pars
    derivexpr <- if (logsd) { # derive f wrt psi=log(sd) !!
        if (density) {
            quote(deriv[cbind(seq_len(L),colidx)] <- exp(-frac) / pi/sdss^2 * (frac-1))
        } else {
            quote(deriv[cbind(seq_len(L),colidx)] <- exp(-frac) * 2*frac)
        }
    } else { # derive f wrt sd !!
        if (density) {
            quote(deriv[cbind(seq_len(L),colidx)] <- exp(-frac) / pi/sdss^3 * (frac-1))
        } else {
            quote(deriv[cbind(seq_len(L),colidx)] <- exp(-frac) * 2*frac/sdss)
        }
    }
    derivexpr <- do.call("substitute", args=list(expr=derivexpr,
                         env=list(colidx=if (nTypes==1L) 1L else quote(types))))
    body(deriv) <- as.call(c(as.name("{"),
        tmp1, tmp2,
        expression(
            deriv <- matrix(0, L, length(pars)),
            frac <- sLengthSquared/2/sdss^2
            ),
        derivexpr,
        expression(deriv)
    ))

    # integrate 'deriv' over a polygonal domain
    body(Deriv) <- as.call(c(as.name("{"),
        ## Determine a = argmax(abs(deriv(c(x,0))))
        if (density) {
            expression(a <- 0)          # maximum absolute value is at 0
        } else {
            c(tmp1, tmp1.1,
              expression(
                  xrange <- polydomain$xrange,           # polydomain is a "owin"
                  a <- min(max(abs(xrange)), sqrt(2)*sd), # maximum absolute value
                  if (sum(xrange) < 0) a <- -a # is more of the domain left of 0?
                  )
              )
        },
        if (nTypes == 1L) {
            expression(deriv.type <- function (s) deriv(s, pars, 1L)[,1L,drop=TRUE])
        } else { # d f(s|type_i) / d sigma_{type_j} is 0 for i != j
            expression(deriv.type <- function (s) deriv(s, pars, type)[,type,drop=TRUE])
        },
        expression(int <- polyCub.SV(polydomain, deriv.type, nGQ=nGQ, alpha=a)),
        if (nTypes == 1L) expression(int) else expression(
            res <- numeric(length(pars)), # zeros
            res[type] <- int,
            res
            )
    ))

    ## sampler (does not obey the 'ub' argument!!)
    body(simulate) <- as.call(c(as.name("{"),
        tmp1, tmp1.1,
        expression(matrix(stats::rnorm(2*n, mean=0, sd=sd), nrow=n, ncol=2L))
    ))

    ## set function environments to the global environment
    environment(f) <- environment(F) <- environment(Fcircle) <-
        environment(deriv) <- environment(Deriv) <-
            environment(simulate) <- .GlobalEnv
    if (is.function(effRange)) environment(effRange) <- .GlobalEnv

    ## return the kernel specification
    list(f=f, F=F, Fcircle=Fcircle, effRange=effRange, deriv=deriv, Deriv=Deriv,
         simulate=simulate, npars=nTypes, validpars=validpars)
}
