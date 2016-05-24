pkgname <- "ExpDes"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ExpDes')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ccboot")
### * ccboot

flush(stderr()); flush(stdout())

### Name: ccboot
### Title: Multiple comparison: Bootstrap
### Aliases: ccboot

### ** Examples

data(ex1)
attach(ex1)
crd(trat, ig, quali = TRUE, mcomp='ccboot', sigF = 0.05)



cleanEx()
nameEx("crd")
### * crd

flush(stderr()); flush(stdout())

### Name: crd
### Title: One factor Completely Randomized Design
### Aliases: crd

### ** Examples

data(ex1)
attach(ex1)
crd(trat, ig, quali = FALSE, sigF = 0.05)



cleanEx()
nameEx("est21Ad")
### * est21Ad

flush(stderr()); flush(stdout())

### Name: est21Ad
### Title: Stink bugs in corn: additional treatment
### Aliases: est21Ad

### ** Examples

data(est21Ad)
## maybe str(est21Ad) ; plot(est21Ad) ...



cleanEx()
nameEx("ex")
### * ex

flush(stderr()); flush(stdout())

### Name: ex
### Title: Vines: Split-Plot in Randomized Blocks Design
### Aliases: ex

### ** Examples

data(ex)
## maybe str(ex) ; plot(ex) ...



cleanEx()
nameEx("ex1")
### * ex1

flush(stderr()); flush(stdout())

### Name: ex1
### Title: Yacon: CRD
### Aliases: ex1

### ** Examples

data(ex1)
## maybe str(ex1) ; plot(ex1) ...



cleanEx()
nameEx("ex2")
### * ex2

flush(stderr()); flush(stdout())

### Name: ex2
### Title: Food bars: RBD
### Aliases: ex2

### ** Examples

data(ex2)
## maybe str(ex2) ; plot(ex2) ...



cleanEx()
nameEx("ex3")
### * ex3

flush(stderr()); flush(stdout())

### Name: ex3
### Title: Forage: LSD
### Aliases: ex3

### ** Examples

data(ex3)
## maybe str(ex3) ; plot(ex3) ...



cleanEx()
nameEx("ex4")
### * ex4

flush(stderr()); flush(stdout())

### Name: ex4
### Title: Composting: Doble Factorial scheme in CRD
### Aliases: ex4

### ** Examples

data(ex4)
## maybe str(ex4) ; plot(ex4) ...



cleanEx()
nameEx("ex5")
### * ex5

flush(stderr()); flush(stdout())

### Name: ex5
### Title: Food bars: Double Factorial scheme in RBD
### Aliases: ex5

### ** Examples

data(ex5)
## maybe str(ex5) ; plot(ex5) ...



cleanEx()
nameEx("ex6")
### * ex6

flush(stderr()); flush(stdout())

### Name: ex6
### Title: Fictional data 1
### Aliases: ex6

### ** Examples

data(ex6)
## maybe str(ex6) ; plot(ex6) ...



cleanEx()
nameEx("ex7")
### * ex7

flush(stderr()); flush(stdout())

### Name: ex7
### Title: Height of corn plants 21 days after emergence.
### Aliases: ex7

### ** Examples

data(ex7)



cleanEx()
nameEx("ex8")
### * ex8

flush(stderr()); flush(stdout())

### Name: ex8
### Title: Composting: double factorial scheme plus one additional
###   treatment in CRD.
### Aliases: ex8

### ** Examples

data(ex8)
## maybe str(ex8) ; plot(ex8) ...



cleanEx()
nameEx("ex9")
### * ex9

flush(stderr()); flush(stdout())

### Name: ex9
### Title: Vegetated: Split-plot in CRD
### Aliases: ex9

### ** Examples

data(ex9)
## maybe str(ex9) ; plot(ex9) ...



cleanEx()
nameEx("fat2.ad.crd")
### * fat2.ad.crd

flush(stderr()); flush(stdout())

### Name: fat2.ad.crd
### Title: Double factorial scheme plus one additional treatment in CRD
### Aliases: fat2.ad.crd

### ** Examples

data(ex8)
attach(ex8)
data(secaAd)
fat2.ad.crd(inoculante, biodiesel, vaso, seca, secaAd, quali = c(TRUE,FALSE), mcomp = "tukey", fac.names = c("Inoculant", "Biodiesel"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.ad.rbd")
### * fat2.ad.rbd

flush(stderr()); flush(stdout())

### Name: fat2.ad.rbd
### Title: Double factorial scheme plus one additional treatment in RBD
### Aliases: fat2.ad.rbd

### ** Examples

data(ex7)
attach(ex7)
data(est21Ad)
fat2.ad.rbd(periodo, nivel, bloco, est21, est21Ad, quali = c(TRUE, FALSE), mcomp = "sk", fac.names = c("Period", "Level"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.crd")
### * fat2.crd

flush(stderr()); flush(stdout())

### Name: fat2.crd
### Title: Double factorial scheme in CRD
### Aliases: fat2.crd

### ** Examples

data(ex4)
attach(ex4)
fat2.crd(revol, esterco, zn, quali=c(FALSE,TRUE), mcomp="tukey", fac.names=c("Revolving","Manure"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat2.rbd")
### * fat2.rbd

flush(stderr()); flush(stdout())

### Name: fat2.rbd
### Title: Double factorial scheme in RBD
### Aliases: fat2.rbd

### ** Examples

data(ex5)
attach(ex5)
fat2.rbd(trat, genero, bloco, sabor ,quali=c(TRUE,TRUE), mcomp="lsd", fac.names=c("Samples","Gender"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.ad.crd")
### * fat3.ad.crd

flush(stderr()); flush(stdout())

### Name: fat3.ad.crd
### Title: Triple factorial scheme plus an additional treatment in CRD
### Aliases: fat3.ad.crd

### ** Examples

data(ex6)
attach(ex6)
data(respAd)
fat3.ad.crd(fatorA, fatorB, fatorC, rep, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "duncan", fac.names = c("Factor A", "Factor B", "Factor C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.ad.rbd")
### * fat3.ad.rbd

flush(stderr()); flush(stdout())

### Name: fat3.ad.rbd
### Title: Triple factorial scheme plus an additional treatment in RBD
### Aliases: fat3.ad.rbd

### ** Examples

data(ex6)
attach(ex6)
data(respAd)
fat3.ad.rbd(fatorA, fatorB, fatorC, rep, resp, respAd, quali = c(TRUE, TRUE, TRUE), mcomp = "snk", fac.names = c("Factor A", "Factor B", "Factor C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.crd")
### * fat3.crd

flush(stderr()); flush(stdout())

### Name: fat3.crd
### Title: Triple factorial scheme in CRD
### Aliases: fat3.crd

### ** Examples

data(ex6)
attach(ex6)
fat3.crd(fatorA, fatorB, fatorC, resp, quali = c(TRUE, TRUE, TRUE), mcomp = "lsdb", fac.names = c("Factor A", "Factor B", "Factor C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("fat3.rbd")
### * fat3.rbd

flush(stderr()); flush(stdout())

### Name: fat3.rbd
### Title: Triple factorial scheme in RBD
### Aliases: fat3.rbd

### ** Examples

data(ex6)
attach(ex6)
fat3.rbd(fatorA, fatorB, fatorC, rep, resp, quali = c(TRUE, TRUE, TRUE), mcomp = "tukey", fac.names = c("Factor A", "Factor B", "Factor C"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("ginv")
### * ginv

flush(stderr()); flush(stdout())

### Name: ginv
### Title: Generalized inverse
### Aliases: ginv

### ** Examples

## Not run: 
# The function is currently defined as
function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, dimnames = dnx[2:1])
}

## End(Not run)



cleanEx()
nameEx("lastC")
### * lastC

flush(stderr()); flush(stdout())

### Name: lastC
### Title: Setting the last character of a chain
### Aliases: lastC

### ** Examples

x<-c("a","ab","b","c","cd")
lastC(x)
# "a" "b" "b" "c" "d"



cleanEx()
nameEx("latsd")
### * latsd

flush(stderr()); flush(stdout())

### Name: latsd
### Title: Latin Square Design
### Aliases: latsd

### ** Examples

data(ex3)
attach(ex3)
latsd(trat, linha, coluna, resp, quali=TRUE, mcomp="snk", sigT=0.05, sigF=0.05)



cleanEx()
nameEx("rbd")
### * rbd

flush(stderr()); flush(stdout())

### Name: rbd
### Title: Randomized Blocks Design
### Aliases: rbd

### ** Examples

data(ex2)
attach(ex2)
rbd(trat, provador, aparencia, quali = TRUE, mcomp='lsd', sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("respAd")
### * respAd

flush(stderr()); flush(stdout())

### Name: respAd
### Title: Fictional data: additional treatment
### Aliases: respAd

### ** Examples

data(respAd)
## maybe str(respAd) ; plot(respAd) ...



cleanEx()
nameEx("secaAd")
### * secaAd

flush(stderr()); flush(stdout())

### Name: secaAd
### Title: Composting: additional treatment
### Aliases: secaAd

### ** Examples

data(secaAd)
## maybe str(secaAd) ; plot(secaAd) ...



cleanEx()
nameEx("split2.crd")
### * split2.crd

flush(stderr()); flush(stdout())

### Name: split2.crd
### Title: Split-plots in CRD
### Aliases: split2.crd

### ** Examples

data(ex9)
attach(ex9)
split2.crd(cobertura, prof, rep, pH, quali = c(TRUE, TRUE), mcomp = "lsd", fac.names = c("Cover", "Depth"), sigT = 0.05, sigF = 0.05)



cleanEx()
nameEx("split2.rbd")
### * split2.rbd

flush(stderr()); flush(stdout())

### Name: split2.rbd
### Title: Split-plots in RBD
### Aliases: split2.rbd

### ** Examples

data(ex)
attach(ex)
split2.rbd(trat, dose, rep, resp, quali = c(TRUE, FALSE), mcomp = "tukey", fac.names = c("Treatament", "Dose"), sigT = 0.05, sigF = 0.05)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
