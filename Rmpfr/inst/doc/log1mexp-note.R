### R code from vignette source 'log1mexp-note.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
## Our custom graphics device:
pdfaCrop <- function(name, width, height, ...) {
    fn <- paste(name, "pdf", sep = ".")
    if(FALSE)## debug
        cat("pdfaCrop: fn = ",fn,"; call:\n\t",deparse(match.call()),"\n")
    grDevices::pdf(fn, width = width, height = height, onefile=FALSE)# ...)
    assign(".pdfaCrop.name", fn, envir = globalenv())
}
## This is used automagically :
pdfaCrop.off <- function() {
    dev.off()# for the pdf
    f <- get(".pdfaCrop.name", envir = globalenv())
    ## and now crop that file:
    pdfcrop <- "pdfcrop" # relying on PATH - fix if needed
    pdftex  <- "pdftex"  # relying on PATH - fix if needed
    system(paste(pdfcrop, "--pdftexcmd", pdftex, f, f, "1>/dev/null 2>&1"),
           intern=FALSE)
}
op.orig <-
options(width = 75,
	SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
	digits = 5,
	useFancyQuotes = "TeX",
	## for JSS, but otherwise MM does not like it:
	## prompt="R> ",
	continue="  ")# 2 (or 3) blanks: use same length as 'prompt'

if((p <- "package:fortunes") %in% search())
    try(detach(p, unload=TRUE, char=TRUE))
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")
if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')
library("sfsmisc")# e.g., for eaxis()
library("Rmpfr")
.plot.BC <- FALSE # no Box-Cox plot


###################################################
### code chunk number 2: def-t3
###################################################
library(Rmpfr)

t3.l1e <- function(a)
{
    c(def   = log(1 - exp(-a)),
      expm1 = log( -expm1(-a)),
      log1p = log1p(-exp(-a)))
}


###################################################
### code chunk number 3: def-leg
###################################################
leg <- local({ r <- body(t3.l1e)[[2]]; r[[1]] <- `expression`; eval(r) })
## will be used below


###################################################
### code chunk number 4: def-test-2 (eval = FALSE)
###################################################
## ##' The relative Error of log1mexp computations:
## relE.l1e <- function(a, precBits = 1024) {
##     stopifnot(is.numeric(a), length(a) == 1, precBits > 50)
##     da <- t3.l1e(a)  ## double precision
##     a. <- mpfr(a, precBits=precBits)
##     ## high precision *and* using the correct case:
##     mMa <- if(a <= log(2)) log(-expm1(-a.)) else log1p(-exp(-a.))
##     structure(as.numeric(1 - da/mMa), names = names(da))
## }


###################################################
### code chunk number 5: def-test-funs
###################################################
library(Rmpfr)

t3.l1e <- function(a)
{
    c(def   = log(1 - exp(-a)),
      expm1 = log( -expm1(-a)),
      log1p = log1p(-exp(-a)))
}
##' The relative Error of log1mexp computations:
relE.l1e <- function(a, precBits = 1024) {
    stopifnot(is.numeric(a), length(a) == 1, precBits > 50)
    da <- t3.l1e(a)  ## double precision
    a. <- mpfr(a, precBits=precBits)
    ## high precision *and* using the correct case:
    mMa <- if(a <= log(2)) log(-expm1(-a.)) else log1p(-exp(-a.))
    structure(as.numeric(1 - da/mMa), names = names(da))
}


###################################################
### code chunk number 6: comp-big (eval = FALSE)
###################################################
## a.s <- 2^seq(-55, 10, length = 256)
## ra.s <- t(sapply(a.s, relE.l1e))


###################################################
### code chunk number 7: bigpic-show (eval = FALSE)
###################################################
## a.s <- 2^seq(-55, 10, length = 256)
## ra.s <- t(sapply(a.s, relE.l1e))
## cbind(a.s, ra.s) # comparison of the three approaches


###################################################
### code chunk number 8: bigpic-do
###################################################
a.s <- 2^seq(-55, 10, length = 256)
ra.s <- t(sapply(a.s, relE.l1e))
capture.and.write(cbind(a.s, ra.s), 8, last = 6)


###################################################
### code chunk number 9: drop-large-a
###################################################
ii <- a.s < 710
a.s  <-  a.s[ii]
ra.s <- ra.s[ii, ]


###################################################
### code chunk number 10: a.small
###################################################
t3.l1e(1e-20)
as.numeric(t3.l1e(mpfr(1e-20, 256)))


###################################################
### code chunk number 11: bigpict-setup (eval = FALSE)
###################################################
## par(mar = c(4.1,4.1,0.6,1.6))
## cc <- adjustcolor(c(4,1,2),.8, red.f=.7)
## lt <- c("solid","33","3262")
## ll <- c(.7, 1.5, 2)


###################################################
### code chunk number 12: bigpict-def (eval = FALSE)
###################################################
## matplot(a.s, abs(ra.s), type = "l", log = "xy",
##         col=cc, lty=lt, lwd=ll, xlab = "a", ylab = "", axes=FALSE)
## legend("top", leg, col=cc, lty=lt, lwd=ll, bty="n")
## draw.machEps <- function(alpha.f = 1/3, col = adjustcolor("black", alpha.f)) {
##     abline(h = .Machine$double.eps, col=col, lty=3)
##     axis(4, at=.Machine$double.eps, label=quote(epsilon[c]), las=1, col.axis=col)
## }
## eaxis(1); eaxis(2); draw.machEps(0.4)


###################################################
### code chunk number 13: zoomin-comp
###################################################
a. <- (1:400)/256
ra <- t(sapply(a., relE.l1e))
ra2 <- ra[,-1]


###################################################
### code chunk number 14: bigpict-fig
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar = c(4.1,4.1,0.6,1.6))
cc <- adjustcolor(c(4,1,2),.8, red.f=.7)
lt <- c("solid","33","3262")
ll <- c(.7, 1.5, 2)
matplot(a.s, abs(ra.s), type = "l", log = "xy",
        col=cc, lty=lt, lwd=ll, xlab = "a", ylab = "", axes=FALSE)
legend("top", leg, col=cc, lty=lt, lwd=ll, bty="n")
draw.machEps <- function(alpha.f = 1/3, col = adjustcolor("black", alpha.f)) {
    abline(h = .Machine$double.eps, col=col, lty=3)
    axis(4, at=.Machine$double.eps, label=quote(epsilon[c]), las=1, col.axis=col)
}
eaxis(1); eaxis(2); draw.machEps(0.4)
## draw the zoom-in region into the plot:
yl <- range(pmax(1e-18, abs(ra2)))
rect(min(a.), yl[1], max(a.), yl[2],
     col= adjustcolor("black", .05), border="gray", pch = 5)


###################################################
### code chunk number 15: zoomin-show (eval = FALSE)
###################################################
## a. <- (1:400)/256
## ra <- t(sapply(a., relE.l1e))
## ra2 <- ra[,-1]


###################################################
### code chunk number 16: boxcox
###################################################
da <- cbind(a = a., as.data.frame(ra2))
library(MASS)
bc1 <- boxcox(abs(expm1) ~ a, data = da, lambda = seq(0,1, by=.01), plotit=.plot.BC)
bc2 <- boxcox(abs(log1p) ~ a, data = da, lambda = seq(0,1, by=.01), plotit=.plot.BC)
c(with(bc1, x[which.max(y)]),
  with(bc2, x[which.max(y)]))## optimal powers
## ==> taking ^ (1/3) :
s1 <- with(da, smooth.spline(a, abs(expm1)^(1/3), df = 9))
s2 <- with(da, smooth.spline(a, abs(log1p)^(1/3), df = 9))


###################################################
### code chunk number 17: zoom-in-def-1 (eval = FALSE)
###################################################
## matplot(a., abs(ra2), type = "l", log = "y", # ylim = c(-1,1)*1e-12,
##         col=cc[-1], lwd=ll[-1], lty=lt[-1],
##         ylim = yl, xlab = "a", ylab = "", axes=FALSE)
## legend("topright", leg[-1], col=cc[-1], lwd=ll[-1], lty=lt[-1], bty="n")
## eaxis(1); eaxis(2); draw.machEps()
## lines(a., predict(s1)$y ^ 3, col=cc[2], lwd=2)
## lines(a., predict(s2)$y ^ 3, col=cc[3], lwd=2)


###################################################
### code chunk number 18: zoom-in-fig
###################################################
getOption("SweaveHooks")[["fig"]]()
cl2 <- "slateblue" # "gray35" # the color for "log(2)"
par(mar = c(4.1,4.1,0.6,1.6))

matplot(a., abs(ra2), type = "l", log = "y", # ylim = c(-1,1)*1e-12,
        col=cc[-1], lwd=ll[-1], lty=lt[-1],
        ylim = yl, xlab = "a", ylab = "", axes=FALSE)
legend("topright", leg[-1], col=cc[-1], lwd=ll[-1], lty=lt[-1], bty="n")
eaxis(1); eaxis(2); draw.machEps()
lines(a., predict(s1)$y ^ 3, col=cc[2], lwd=2)
lines(a., predict(s2)$y ^ 3, col=cc[3], lwd=2)

abline(v = log(2), col = cl2, lty=4)
axis(1, at=log(2), label=quote(a[0] == log~2), las=1,
     col.axis=cl2, col=cl2, lty=4)
## what system is it ?
sysInf <- Sys.info()[c("sysname", "release", "nodename", "machine")]
mtext(with(as.list(sysInf),
	   paste0(sysname," ",release,"(",substr(nodename,1,16),") -- ",
		  machine)),
      side=1, adj=1, line=2.25, cex = 3/4)


###################################################
### code chunk number 19: uniroot-x1
###################################################
## Find x0, such that  exp(x) =.= g(x) for x < x0 :
f0 <- function(x) { x <- exp(x) - log1p(exp(x))
                   x[x==0] <- -1 ; x }
u0 <- uniroot(f0, c(-100, 0), tol=1e-13)
str(u0, digits=10)
x0 <- u0[["root"]] ## -36.39022698 --- note that ~= \log(\eps_C)
all.equal(x0, -52.5 * log(2), tol=1e-13)

## Find x1, such that  x + exp(-x) =.= g(x) for x > x1 :
f1 <- function(x) { x <- (x + exp(-x)) - log1p(exp(x))
                   x[x==0] <- -1 ; x }
u1 <- uniroot(f1, c(1, 20), tol=1e-13)
str(u1, digits=10)
x1 <- u1[["root"]] ## 16.408226

## Find x2, such that  x =.= g(x)  for x > x2 :
f2 <- function(x) { x <- log1p(exp(x)) - x ; x[x==0] <- -1 ; x }
u2 <- uniroot(f2, c(5, 50), tol=1e-13)
str(u2, digits=10)
x2 <- u2[["root"]] ## 33.27835


###################################################
### code chunk number 20: log1pexp-plot
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfcol= 1:2, mar = c(4.1,4.1,0.6,1.6), mgp = c(1.6, 0.75, 0))
curve(x+exp(-x) - log1p(exp(x)), 15, 25,   n=2^11); abline(v=x1, lty=3)
curve(log1p(exp(x)) - x,       33.1, 33.5, n=2^10); abline(v=x2, lty=3)


###################################################
### code chunk number 21: def-test-pfuns
###################################################
t4p.l1e <- function(x)
{
    c(def   = log(1 + exp(x)),
      log1p = log1p(exp(x)),
      ## xlog1p = x + log1p(exp(-x)),
      xpexp = x + exp(-x),
      x = x)
}
leg <- local({ r <- body(t4p.l1e)[[2]]; r[[1]] <- `expression`; eval(r) })
##' The relative Error of log1pexp computations:
relE.pl1e <- function(x, precBits = 1024) {
    stopifnot(is.numeric(x), length(x) == 1, precBits > 50)
    dx <- t4p.l1e(x)  ## double precision
    x. <- mpfr(x, precBits=precBits)
    ## high precision *and* using the correct case:
    mMx <- if(x < 0) log1p(exp(x.)) else x. + log1p(exp(-x.))
    structure(as.numeric(1 - dx/mMx), names = names(dx))
}


###################################################
### code chunk number 22: comp-big
###################################################
x.s <- seq(-100, 750, by = 5)  # <- the big picture ==> problem for default
x.s <- seq( 5, 60, length=512) # <- the zoom in     ==> *no* problem for def.
rx.s <- t(sapply(x.s, relE.pl1e))
signif(cbind(x.s, rx.s),3)


###################################################
### code chunk number 23: bigpict-2-fig
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar = c(4.1,4.1,0.6,1.6), mgp = c(1.6, 0.75, 0))
cc <- adjustcolor(c(4,1,2,3),.8, red.f=.7, blue.f=.8)
lt <- c("solid","33","3262","dotdash")
ll <- c(.7, 1.5, 2, 2)
ym <- 1e-18
yM <- 1e-13
matplot(x.s, pmax(pmin(abs(rx.s),yM),ym), type = "l", log = "y", axes=FALSE,
        ylim = c(ym,yM), col=cc, lty=lt, lwd=ll, xlab = "x", ylab = "")
legend("topright", leg, col=cc, lty=lt, lwd=ll, bty="n")
eaxis(1, at=pretty(range(x.s), n =12)); eaxis(2)
draw.machEps(0.4)
x12 <- c(18, 33.3)
abline(v=x12, col=(ct <- adjustcolor("brown", 0.6)), lty=3)
axis(1, at=x12, labels=formatC(x12), padj = -3.2, hadj = -.1, tcl = +.8,
     col=ct, col.axis=ct, col.ticks=ct)


###################################################
### code chunk number 24: sessionInfo
###################################################
toLatex(sessionInfo(), locale=FALSE)


###################################################
### code chunk number 25: finalizing
###################################################
options(op.orig)


