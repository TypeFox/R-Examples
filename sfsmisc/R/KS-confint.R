### Fixme: In the following, computing and plotting should be separated

###--> ./ecdf.R plot.ecdf() should get conf.type and conf.int argument!!

### Also, I've posted a pre-version of this:
## Date: Mon, 22 Oct 2001 19:15:35 +0200
## From: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: [R] Re: conf.int. for ecdfs {was "Two questions"}
## To: cblouin@is2.dal.ca
## Cc: Kjetil Halvorsen <kjetilh@umsanet.edu.bo>, r-help@stat.math.ethz.ch

### Note -- this is related to the pkstwo() function inside ks.test()
### ====    in stats : ~/R/D/r-devel/R/src/library/stats/R/ks.test.R

ecdf.ksCI <- function(x, main = NULL, sub = NULL,
                      xlab = deparse(substitute(x)), ci.col = "red", ...)
{
    force(xlab)
    if(is.null(main))
        main <- paste0("ecdf(",deparse(substitute(x)),") + 95% K.S. bands")
    n <- length(x)
    if(is.null(sub))
        sub <- paste("n = ", n)
    ec <- ecdf(x)
    xx <- get("x", envir=environment(ec))# = sort(x)
    yy <- get("y", envir=environment(ec))
    D <- KSd(n)
    yyu <- pmin(yy + D, 1)
    yyl <- pmax(yy - D, 0)
    ecu <- stepfun(xx, c(yyu, 1) )
    ecl <- stepfun(xx, c(yyl, yyl[n]) )

    ## Plots -- all calling  plot.stepfun

    plot(ec, main = main, sub = sub, xlab = xlab, ...)
    plot(ecu, add=TRUE, verticals=TRUE, do.points=FALSE,
         col.hor= ci.col, col.vert= ci.col, ...)
    plot(ecl, add=TRUE, verticals=TRUE, do.points=FALSE,
         col.hor= ci.col, col.vert= ci.col, ...)
}

KSd <- function(n)
{
    ## `approx.ksD()'
    ## approximations for the critical level for Kolmogorov-Smirnov statistic D,
    ## for confidence level 0.95. Taken from Bickel & Doksum, table IX, p.483
    ## and Lienert G.A.(1975) who attributes to Miller,L.H.(1956), JASA
    ifelse(n > 80,
           1.358/( sqrt(n) + .12 + .11/sqrt(n)),## Bickel&Doksum, table IX,p.483

           splinefun(c(1:9, 10, 15, 10 * 2:8),# from Lienert
                     c(.975,   .84189, .70760, .62394, .56328,# 1:5
                       .51926, .48342, .45427, .43001, .40925,# 6:10
                       .33760, .29408, .24170, .21012,# 15,20,30,40
                       .18841, .17231, .15975, .14960)) (n))
}
