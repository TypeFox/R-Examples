p.rbwheel <- function(n, p, frac = 1/p,
                    sig1 = 0.05, sig2 = 1/10, ...)
{
    ## Purpose: rbwheel() +  plot
    ## ----------------------------------------------------------------------
    ## Arguments: as rbwheel(.)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 14 Jun 2009, 16:48

    call <- match.call()
    call[[1]] <- as.name("rbwheel")
    Lab <- paste("X <-", deparse(call, width.cutoff= 200))

    rb <- rbwheel(n=n, p=p, frac=frac, sig1=sig1, sig2=sig2, ...,
                  fullResult = TRUE)
    ##            ^^^^^^^^^^^^^^^^^

    par.s <- list(plot.symbol = list(alpha=0.4, cex= .5, pch = 16))

    pX <- splom(~ rb$X, xlab = Lab, pscales = 0,
                par.settings = list(plot.symbol = c(par.s$plot.symbol)))
    pX0 <- splom(~ rb$X0, pscales = 0,
                par.settings = list(plot.symbol = c(par.s$plot.symbol,
                                    col = "dark gray")),
                xlab = "X0 *before* rotation (and scaling)")

    ## 2 x 1  plots :  [ ]   [ ]
    print(pX,  split=c(1,1, 2,1), more=TRUE)
    print(pX0, split=c(2,1, 2,1))
    invisible(rb)
}

library(robustX)
library(lattice)

p.rbwheel(500,4)

cat("Now you can experiment yourself!   ...\n",

      "etc ..\n", sep="\n")

if(FALSE) {
pdf("/u/maechler/R/Meetings-Kurse-etc/2009-Parma/p_rbwheel-500-6.pdf", height=4,width=8)
p.rbwheel(500,6)
dev.off()
}
