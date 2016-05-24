### R code from vignette source 'rhoAMH-dilog.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
## Custom graphics device (for cropping .pdf):
pdfCrop <- function(name, width, height, ...) {
    f <- paste0(name, ".pdf")
    grDevices::pdf(f, width=width, height=height, onefile=FALSE)
    assign(".pdfCrop.name", f, envir=globalenv())
}
pdfCrop.off <- function() { # used automagically
    grDevices::dev.off() # closing the pdf device
    f <- get(".pdfCrop.name", envir=globalenv())
    system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"),
           intern=FALSE) # crop the file (relies on PATH)
}
op.orig <-
options(width = 70, useFancyQuotes = FALSE,
        ## SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        prompt="> ",  continue="   ")
## JSS: prompt="R> ", continue=">  ")
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")
## if(Sys.getenv("USER") == "maechler")# take CRAN's version, not development one:
##    require("copula", lib="~/R/Pkgs/CRAN_lib")
require("copula")


###################################################
### code chunk number 2: ak-frac1 (eval = FALSE)
###################################################
## require(sfsmisc) #--> mat2tex(), mult.fig(), eaxis()
## k <- 1:9; ak <- MASS::fractions(12/((k+1)*(k+2))^2)


###################################################
### code chunk number 3: ak-frac2 (eval = FALSE)
###################################################
## rbind(k = k, `$a_k$` = as.character(ak))


###################################################
### code chunk number 4: ak-frac-show (eval = FALSE)
###################################################
## require(sfsmisc) #--> mat2tex(), mult.fig(), eaxis()
## k <- 1:9; ak <- MASS::fractions(12/((k+1)*(k+2))^2)
## rbind(k = k, `$a_k$` = as.character(ak))


###################################################
### code chunk number 5: ak-frac-do
###################################################
require(sfsmisc) #--> mat2tex(), mult.fig(), eaxis()
k <- 1:9; ak <- MASS::fractions(12/((k+1)*(k+2))^2)
mat2tex(
rbind(k = k, `$a_k$` = as.character(ak))
      , stdout())


###################################################
### code chunk number 6: def-rhos
###################################################
##' Version 1:  Direct formula from Nelsen:
.rhoAmh.1 <- function(a) {
  Li2 <- gsl::dilog(a)
  12 * (1 + a) / a^2 * Li2 - 24 * (1 - a) / a^2 * log1p(- a) - 3 * (a + 12) / a
}
.rhoAmh.1b <- function(a) {
  Li2 <- gsl::dilog(a)
  ## factored out 3/a from version 1:
  3/a * (4 * (1 + a) / a * Li2 - 8 * (1 - a) / a * log1p(- a) - (a + 12))
}

##' Version 2:
.rhoAmh.2 <- function(a, e.sml = 1e-11) {
  stopifnot(length(a) <= 1)
  if(abs(a) < e.sml) { ## if |a| << 1, do better than the direct formula:
      a*(1/3 + a*(1/12 + a*(3/100 + a/75)))
  } else { ## regular a
      Li2 <- gsl::dilog(a)
      3/a * (4 * (1 + 1/a) * Li2 - 8 * (1/a - 1) * log1p(- a) - (a + 12))
  }
}

##' Series version with N terms:
rhoAmh.T <- function(a, N) {
    stopifnot(length(N) == 1, N == as.integer(N), N >= 1)
    if(N <= 4)
        switch(N,
               a/3,
               a/3*(1 + a/4),
               a*(1/3 + a*(1/12 + a* 3/100)),
               a*(1/3 + a*(1/12 + a*(3/100 + a/75))))
    else { ## N >= 5
        n <- N:1 #--> sum smallest to largest
        if(is(a, "mpfr")) ## so all computations work in high precision
            n <- mpfr(n, precBits=max(.getPrec(a)))
        cf <- ## 3/choose(n+2, 2)^2
            3/((n+1)*(n+2)/2)^2
        a2n <- outer(n,a, function(x,y) y^x) ## a2n[i,j] := a[j] ^ n[i]
        colSums(cf * a2n)
    }
}


###################################################
### code chunk number 7: first-curve
###################################################
r1  <- curve( .rhoAmh.1 (x),  1e-20, .1, log="x", n=1025)
r1b <- curve( .rhoAmh.1b(x), n=1025, add=TRUE, col=2)
r2 <-  curve( Vectorize(.rhoAmh.2)(x), n=1025, add=TRUE,
             col=adjustcolor("blue4",1/4), lwd = 5)
tab <- cbind(as.data.frame(r1), y.b = r1b$y, y2 = r2$y)


###################################################
### code chunk number 8: more-curves
###################################################
if(require("sfsmisc")) {
    myAxes <- function(sides) for(s in sides) eaxis(s)
} else {
    myAxes <- function(sides) for(s in sides)  axis(s)
}
rhoAcurve <- function(k, ..., log = "",
                      ylab = substitute({rho^"*"}[2](x,KK), list(KK=k)))
    curve(Vectorize(.rhoAmh.2)(x, k), n=1025, ylab=ylab, log=log,
          xaxt = if(grepl("x", log, fixed=TRUE)) "n" else "s",
          yaxt = if(grepl("y", log, fixed=TRUE)) "n" else "s", ...)

e.s <- eval(formals(.rhoAmh.2)$e.sml); t0 <- e.s * .99999
op <- sfsmisc::mult.fig(2, marP = -c(1.4,1,1,1))$old.par
rhoAcurve(e.s, 1e-18, 1e-1, log = "xy", ylab=""); myAxes(1:2)
lines(t0, .rhoAmh.2(t0), type="h", lty=3, lwd = 3/4)
rhoAcurve(1e-6, add=TRUE, col=adjustcolor(2, 1/3), lwd=4)
rhoAcurve(1e-6, 1e-18, 1, log="x", col="tomato"); myAxes(1)
par(op)


###################################################
### code chunk number 9: curve4
###################################################
rhoAcurve(1e-6, 1e-7, 1e-5, log = "y", col="tomato"); myAxes(2)
abline(v=1e-6, lty=3, lwd=1/2)


###################################################
### code chunk number 10: curve5
###################################################
cc <- 1e-4 ; op <- mult.fig(2, marP= -c(1,0,1,1))$old.par
rhoAcurve(cc, 1e-6, 1e-3, log = "xy", col="tomato",ylab=""); myAxes(1:2)
abline(v=cc, lty=3, lwd=1/2)
## zoom in extremely:
rhoAcurve(cc, cc*(1-1e-4), cc*(1+1e-4), col="tomato")
abline(v=cc, lty=3, lwd=1/2);          par(op)


###################################################
### code chunk number 11: curve7
###################################################
cc <- 1e-3
rhoAcurve(cc, cc*(1-2^-20), cc*(1+2^-20), log="y",yaxt="s", col="tomato")
abline(v=cc, lty=3, lwd=1/2)
rhoAcurve(cc*10, add=TRUE, col=adjustcolor(1,.25), lwd=3)


###################################################
### code chunk number 12: curve8
###################################################
cc <- 0.01
rhoAcurve(cc, cc*(1-2^-20), cc*(1+2^-20), log="y",yaxt="s", col="tomato")
abline(v=cc, lty=3, lwd=1/2)
rhoAcurve(cc*10, add=TRUE, col=adjustcolor(1,.25),lwd=5)


###################################################
### code chunk number 13: Taylor-curves-1
###################################################
a <- 2^seq(-30,-1, by = 1/32)# 0 < a <= 0.5
rhoA.T <-  vapply(1:6, rhoAmh.T, a=a, numeric(length(a)))
op <- mult.fig(mfcol=c(1,3), mgp=c(2.5,.8,0))$old.par
matplot(a, rhoA.T, type="l")
matplot(a, rhoA.T, type="l", log="y", yaxt="n")   ; myAxes(2)
matplot(a, rhoA.T, type="l", log="xy", axes=FALSE); myAxes(1:2);box()
par(op)


###################################################
### code chunk number 14: Taylor-curves-2
###################################################
rhoA.true <- rhoAmh.T(a,50)
chk.w.mpfr <- FALSE ## Sys.info()[["user"]] == "maechler"
if(chk.w.mpfr) {
    require(Rmpfr)## get the "really" "true" values:
    print(system.time(rhA.mp  <- rhoAmh.T(mpfr(a, prec=256), 50))) ## 3.95 sec (lynne)
    print(system.time(rhA.mp1 <- rhoAmh.T(mpfr(a, prec=256), 60))) ## 4.54 sec
    stopifnot(all.equal(rhA.mp, rhoA.true, tol = 1e-15))
        print(all.equal(rhA.mp, rhoA.true, tol = 1e-20)) ## 6.99415....e-17 [64bit, lynne]
    ## see if the 50 terms have converged:
    print( all.equal(rhA.mp, rhA.mp1, tol = 1e-30) )
    ## "Mean relative difference: 2.4958....e-22"
    ## ==> 50 terms seem way enough for double prec
}
matplot(a, 1 - rhoA.T / rhoA.true, type="l", log="y")


###################################################
### code chunk number 15: def-plot-relE
###################################################
pl.relE.rhoAMH <- function(N.max, N.inf = 50, N.min = 1, l2a.lim = c(-30, -1),
                           n.p.u = 2^round(log2(1000 / diff(l2a.lim))),
                           cut.rA2 = 1e-7,
                           colX = adjustcolor("midnightblue", 0.5), ...)
{
    stopifnot(length(l2a.lim) >= 2, l2a.lim < 0, n.p.u >= 1,
              N.max >= N.min, N.min >= 1, N.inf > N.max + 4,
              (N3 <- c(N.min, N.max, N.inf)) == as.integer(N3))
    a <- 2^seq(l2a.lim[1], l2a.lim[2], by = 1/n.p.u)
    N.s <- N.min:N.max
    rhoA.true <- rhoAmh.T(a, N.inf)
    rhoA.T <- vapply(N.s, rhoAmh.T, a=a, numeric(length(a))) # matrix
    rhoA.v2 <- Vectorize(.rhoAmh.2)(a, cut.rA2) # "Li2()+direct" below

    ## matplot() compatible colors and lty's
    cols <- palette()[1 + (N.s-1) %% 6]
    ltys <- (1:5)    [1 + (N.s-1) %% 5]
    matplot(a, 1 - rhoA.T / rhoA.true, type="l", log="xy",
            col=cols, lty=ltys, axes=FALSE, frame=TRUE, ...)
    myAxes(1:2)
    lines(a, 1 - rhoA.v2 / rhoA.true, col= colX, lwd=3)
    legend("topleft", c(paste0("N=",N.s), "Li2()+direct"),
           col=c(cols, colX), lty=c(ltys, 1), lwd=c(rep(1,length(N.s)), 3),
           cex=.75, bty="n")
    invisible(list(a=a, rhoA.T=rhoA.T, rhoA.v2 = rhoA.v2))
}


###################################################
### code chunk number 16: relE-curves-1-2
###################################################
op <- mult.fig(2, marP=-c(1.5,1.5,2,1))$old.par
pl.relE.rhoAMH(4, l2a=c(-53,-1), ylab="")
pl.relE.rhoAMH(6,                ylab="")


###################################################
### code chunk number 17: relE-curves-3-4
###################################################
mult.fig(2, marP=-c(1.5,1.5,2,1))
pl.relE.rhoAMH(12, l2a=c(-12, -3),            ylab="")
pl.relE.rhoAMH(20, l2a=c(-8, -.5), N.min = 4, ylab="")


###################################################
### code chunk number 18: relE-curves-5
###################################################
par(op); pl.relE.rhoAMH(40, l2a=c(-5, -.5), N.min = 10)


###################################################
### code chunk number 19: relE-C-1
###################################################
pl.relE.rhoAMH(6)
abline(v=1e-4, col="gray", lty=2)#-> N=2 cutoff
abline(v=2e-3, col="gray", lty=2)#-> N=3 cutoff


###################################################
### code chunk number 20: relE-C-2
###################################################
pl.relE.rhoAMH(12, l2a=c(-12, -3))
abline(v= 2e-3, col="gray", lty=2)#-> N=3 cutoff
abline(v= 7e-3, col="gray", lty=2)#-> N=4 cutoff
abline(v=16e-3, col="gray", lty=2)#-> N=5 cutoff


###################################################
### code chunk number 21: rhoAmh-def-show (eval = FALSE)
###################################################
## copula ::: .rhoAmhCopula


###################################################
### code chunk number 22: rhoAmh-def-do
###################################################
rhoA <- copula ::: .rhoAmhCopula; environment(rhoA) <- environment()
rhoA


###################################################
### code chunk number 23: remnant-fig
###################################################
rhoAMH <- Vectorize(copula:::.rhoAmhCopula)
curve(rhoAMH, n=1025, -1, 1, ylim= c(-1,1), xlab = quote(theta),
      ylab="", col="tomato", lwd=2, las=1)
abline(0, 1/3, lty=2, col=(adjustcolor(c2 <- "orange2", 2/3)))
curve(x/3*(1+x/4), lty=2, col=(adjustcolor(c3 <- "blue", 1/2)),
      -1.1,1.1, add=TRUE); x. <- .65
text(.4 , .3 , quote(rho[plain(AMH)](theta)),col="tomato")
text(.88, .23, quote(y == theta/3), col=c2)
text(.7,  .05, quote(y == theta/3*(1+theta/4)), col=adjustcolor(c3, 1/2))
segments(.55, .10, x., x./3*(1+x./4), lty="82", col=adjustcolor(c3, 1/2))
abline(h=0,v=0, lty=3); rect(-1,-1,1,1, lty=3)


###################################################
### code chunk number 24: regression-tests
###################################################
t0 <- seq(-1,1, by=2^-8)[1:512]
t1 <- seq(-1/2, 1/2, by = 2^-8)
th <- 10^-(6:99); i <- -(1:9)
rth <- rhoAMH(th)
stopifnot(all.equal(rhoAMH(1), 4*pi^2 - 39, tol = 8e-15),# <- gave NaN
          all.equal(rhoAMH(t0), t0/3 * (1 + t0/4), tol = 0.06),
          all.equal(rhoAMH(t1), t1/3 * (1 + t1/4), tol = 1/85),
          all.equal(rth,      th / 3 * (1 + th/4), tol = 1e-15),
          all.equal(rth,      th / 3, tol = 1e-6),
          all.equal(rth[i], th[i]/ 3, tol = 6e-16))
th <- 10^-(16:307)
stopifnot(all.equal(th/3, rhoAMH(th), tol=4e-16),
          rho(amhCopula(0, use.indepC="FALSE")) == 0)


###################################################
### code chunk number 25: sessionInfo
###################################################
toLatex(sessionInfo(), locale=FALSE)


###################################################
### code chunk number 26: copula-version (eval = FALSE)
###################################################
## my.strsplit(  packageDescription("copula")[["Date"]]  )


###################################################
### code chunk number 27: copula-version
###################################################
pd <- gsub("\\$", '', packageDescription("copula")[["Date"]])
p2 <- strsplit(sub("^Date: +", '', pd),"\n")[[1]][2:1]
cat(p2[1], "-- ", sub(" *, *$", '', p2[2]), "\n", sep="")


###################################################
### code chunk number 28: finalizing
###################################################
options(op.orig)


