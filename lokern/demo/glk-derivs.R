## Q: Can we have *same* kernel, *same* bandwidth  with different  'deriv'
##    similarly to smooth.spline() ?
##
## Answer: not really,  mainly because don't have enough choices
##     (nu, k_{ord}), i.e., because currently,  nu - k  must be even ...

## "dput(.) of simple integer vector to character :
myF <- function(d, A="(", O=")", B = length(d) > 1)
    paste(if(B) A, paste(d, collapse=", "), if(B) O, sep="")

library(lokern)

p.3glks <- function(x.dat, y.dat, korder, derivs = 0:2,
                    is.rand=FALSE, useBandwidth, bw.factor = 1.8,
                    col = 2, lwd = 1.5)
{
    ## Purpose: Plot  glkerns(*,  deriv = {0, 1, 2})
    ## ----------------------------------------------------------------------
    ## Arguments: (x.dat, y.dat): the numeric data vectors
    ##            korder : the kernel order -- automatically is diminuished by one
    ##                     if needed to keep  'korder - deriv' an even number
    ##            derivs : integer vectors of derivatives to compute
    ##            useBandwidth: possibly a user specified bandwidth
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  2 Jul 2009, 09:24

    if(!missing(useBandwidth) && is.numeric(useBandwidth) && useBandwidth > 0)
        bw <- useBandwidth
    else {
        ## Determine the fixed bandwidth :
        bw0 <- glkerns(x.dat, y.dat, korder=korder, is.rand=is.rand)$bandwidth
        bw <- bw0 * bw.factor     # more smoothing for the derivatives
    }

    stopifnot(derivs == (d <- as.integer(derivs)), length(derivs) >= 1)
    derivs <- d
    stopifnot(0 <= derivs, derivs <= 4)
    glist <- as.list(derivs)
    ## Estimates for   g, g' , g'' ... {well, depending on derivs} :
    for(i.d in seq_along(derivs)) {
        nu <- derivs[i.d]
        k0 <- korder
        if((korder - nu) %% 2 != 0) { ## 'korder - nu' must be even {theory; Fortran code}
            k0 <- korder - 1
            message(gettextf("deriv = %d: modifying korder from %d to %d",
                             nu, korder, k0))
        }
        glist[[i.d]] <-
            glkerns(x.dat, y.dat, korder=k0, is.rand=is.rand, bandw = bw, deriv = nu)
    }

    names(glist) <- paste("nu=", derivs, sep="")

    ## ---------  Plots ---------------

    op <- par(mfrow= c(length(derivs), 1), mgp = c(1.25, 0.6, 0),
              mar = c(3,3,2.5,1) + .1, oma = c(0,0, 2, 0))
    on.exit(par(op))

    for(i.d in seq_along(derivs)) {
        nu <- derivs[i.d]
        tit <-
            switch(nu + 1,
                   expression(widehat(g)(.)), # 0
                   expression(widehat(g * minute)(.)), # 1: g'
                   expression(widehat(g * second)(.)), # 2: g''
                   expression(widehat(g * minute*second)(.)),# 3: g'''
                   expression(widehat(g ^ {(4)})),
                   expression(widehat(g ^ {(5)})),
                   expression(widehat(g ^ {(6)})))

        with(glist[[i.d]], {
            plot(est ~ x.out, type = "l", main = tit, col=col, lwd=lwd)
            if(nu == 0) ## data
                points(y ~ x, cex = 0.5)
            else ## y = 0 line (to see zeros):
                abline(h = 0, col = "gray", lty=3)

            mtext(substitute(list(bw == B,k[ord] == K),
                             list(B = formatC(bandwidth), K = korder)),
                  adj = 1, line = .25) })
    }
    mtext(sprintf("glkerns(*, deriv = %s, bandwidth = <fixed>, korder = %d)",
                  myF(derivs, "{","}"), korder),
          line = 0.5, outer = TRUE, cex = par("cex.main"), font = par("font.main"))

    invisible(glist)
}

data(xSim)
n <- length(xSim)
tt <- ((1:n) - 1/2)/n # equidistant x

p.3glks(tt, xSim, kord = 4)

## Chose bandwidth yourself; see all available derivatives:
## Store results
r <- p.3glks(tt, xSim, kord = 6,
             derivs = 0:4, useBand = 0.15)
## and inspect them
str(r)
