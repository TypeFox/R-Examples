rbw.ics <- function(n, p, frac = 1/p,
                    sig1 = 0.05, sig2 = 1/10, ...)
{
    ## Purpose: rbwheel() +  3 ICS() methods ... mainly for nice plotting
    ## ICS: Compare M-estimate [Max.Lik. of t_{df = 2}] with high-breakdown
    ## ----------------------------------------------------------------------
    ## Arguments: as rbwheel(.)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 13 Jun 2009, 20:00

    call <- match.call()
    call[[1]] <- as.name("rbwheel")
    Lab <- paste("X <-", deparse(call, width.cutoff= 200))
    X <- rbwheel(n=n, p=p, frac=frac, sig1 = 0.05, sig2 = 1/10, ...)
    X.paM <- ics(X, S1=cov, S2= function(.) cov.trob(., nu=2)$cov, stdKurt=FALSE)
    X.paM.<- ics(X, S1=cov, S2= function(.) tM(., df=2)$V, stdKurt = FALSE)
    X.paR <- ics(X, S1=cov, S2= function(.) covMcd(.)$cov, stdKurt = FALSE)

    par.s <- list(plot.symbol = list(alpha=0.4, cex= .5, pch = 16))
    pX <- splom(~ X, xlab = Lab, pscales = 0,
                par.settings = list(plot.symbol = c(par.s$plot.symbol,
                                    col = "dark gray")))
    p1 <- splom(~ X.paM @ Scores, par.settings=par.s,
                pscales = 0, xlab = "ics(X, .. S2= cov.trob(., nu=2))")
    p2 <- splom(~ X.paM.@ Scores, par.settings=par.s,
                pscales = 0, xlab = "ics(X, .. S2= tM(., df=2))")
    pMcd <- splom(~ X.paR @ Scores, par.settings=par.s,
                pscales = 0, xlab = "ics(X, .. S2= covMcd(.))")

    structure(list(X=X, pX=pX, p1=p1, p2=p2, pMcd=pMcd), class = "rbwIcs4plot")
}

## the (lattice-like) print method for such an object:
print.rbwIcs4plot <- function(x, ...) {
    print(x$pX,  split=c(1,1, 2,2), more=TRUE)
    print(x$p1,  split=c(1,2, 2,2), more=TRUE)
    print(x$p2,  split=c(2,2, 2,2), more=TRUE)
    print(x$pMcd,split=c(2,1, 2,2))
    invisible(x)
}


library(lattice)

if(require("ICS")) {
  stopifnot(require("MASS"))

  print(rr <- rbw.ics(200, p=4))

  cat("Now you can experiment yourself!   E.g.,\n",
      "     rbw.ics(500, p=5, sig2= 0.3)",
      "     rbw.ics(1000, p=4, spherize=TRUE)",
      "     rbw.ics(100, p=10)\n",
      "etc ..\n", sep="\n")
} else {
    cat("You really should install the  'ICS'  package for this demo",
        "see  ?install.packages\n", sep="\n")
}
