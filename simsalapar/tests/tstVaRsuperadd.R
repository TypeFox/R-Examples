#### A "fast" and simplified version of demo(VaRsuperadd)
#### see ../demo/VaRsuperadd.R
####     ~~~~~~~~~~~~~~~~~~~~~

require(simsalapar)
source(system.file("xtraR/assertErr-etc.R", package="simsalapar", mustWork=TRUE))

## Must be fast, rather than "interesting":
n.obs  <- 16
n.alpha <- 8
doExtras <- FALSE

demo(VaRsuperadd) # crucial part!

## The only difference to doOne() is the missing t(.) at the very end
do1.wrong <- function(n, d, family, tau, qmargin, alpha)
{
    cop <- switch(family,
                "normal" =
                  ellipCopula("normal", param=iTau(ellipCopula("normal"), tau=tau),
                              dim=d),
                "t" =
                  ellipCopula("t", param=iTau(ellipCopula("t"), tau=tau), dim=d),
                "Clayton" =
                  onacopulaL("Clayton", list(th=iTau(archmCopula("clayton"), tau),
                                             seq_len(d))),
                "Gumbel" =
                  onacopulaL("Gumbel", list(th=iTau(archmCopula("gumbel"), tau),
                                            seq_len(d))),
                stop("unsupported 'family'"))
    U <- rCopula(n, copula=cop)
    sapply(qmargin, function(FUN) ecdf(rowSums(FUN(U)))(d*FUN(alpha)) - alpha)
}

assertError( doCheck(do1.wrong, varList) )
