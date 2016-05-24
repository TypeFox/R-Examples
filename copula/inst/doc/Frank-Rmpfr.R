### R code from vignette source 'Frank-Rmpfr.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
op.orig <-
options(width = 75,
        SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        useFancyQuotes = FALSE,
        ## for JSS, but otherwise MM does not like it:
        ## prompt="R> ",
        continue="  ")# 2 (or 3) blanks: use same length as 'prompt'
copDDir <- system.file('doc', package='copula')

Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")


###################################################
### code chunk number 2: nacopula-dDiagA (eval = FALSE)
###################################################
## copula:::dDiagA


###################################################
### code chunk number 3: nacopula-dDiagA-show
###################################################
writeLines(head(capture.output(print(
copula:::dDiagA
                     )), -1))


###################################################
### code chunk number 4: my-dDiagA
###################################################
dDiagA <- function(u, th, d, iPsi, absdPsi, absdiPsi, log = FALSE) {
    stopifnot(is.finite(th), d > 0, is.function(iPsi),
              is.function(absdPsi), is.function(absdiPsi))
    if(log) {
        log(d) + absdPsi(d*iPsi(u,th), th, log = TRUE) +
            absdiPsi(u, th, log = TRUE)
    } else {
        d* absdPsi(d*iPsi(u,th), th) * absdiPsi(u,th)
    }
}


###################################################
### code chunk number 5: def-orig-psi-func
###################################################
iPsi.0 <- function(u,theta) -log( (exp(-theta*u)-1) / (exp(-theta)-1) )
iPsi.1 <- function(u,theta) -log(expm1(-u*theta) / expm1(-theta))
absdiPsi.1 <- function(u, theta, log = FALSE)
  if(log) log(theta)-log(expm1(u*theta)) else theta/expm1(u*theta)


###################################################
### code chunk number 6: def-orig-absdPsi
###################################################
require("copula")# for polylog()
absdPsi.1 <- function(t, theta, log=FALSE) {
  p <- -expm1(-theta)
  Li. <- polylog(log(p) - t, s = 0, log=log,
                 method="negI-s-Eulerian", is.log.z=TRUE)
  if(log) Li. - log(theta) else Li. / theta
}


###################################################
### code chunk number 7: simpler-absdPsi
###################################################
absdPsi.2 <- function(t, theta, log=FALSE) {
  w <- log(-expm1(-theta)) - t
  Li. <- if(log) w - log(-expm1(w)) else -exp(w)/expm1(w)
  if(log) Li. - log(theta) else Li. / theta
}


###################################################
### code chunk number 8: dDiag-probl-big-theta
###################################################
curve(dDiagA(x, th = 38, d = 2, iPsi = iPsi.0,
             absdPsi=absdPsi.1, absdiPsi=absdiPsi.1, log = TRUE),
      ylab = "dDiagA(x, theta= 38, *, log=TRUE)",
      0, 1, col = 4, n = 1000)
## and using the slightly better   iPsi.1  does not help here:
curve(dDiagA(x, th = 38, d = 2, iPsi = iPsi.1,
             absdPsi=absdPsi.2, absdiPsi=absdiPsi.1, log = TRUE),
      add = TRUE, col = 2, n=1000)
legend("bottom", c("iPsi.0()","iPsi.1()"),col=c(4,2), lty=1, bty="n")


###################################################
### code chunk number 9: dDiag-ok-big-theta
###################################################
iPsi.2 <- function(u,theta) -log1p((exp(-u*theta)-exp(-theta)) / expm1(-theta))
curve(dDiagA(x, th = 38, d = 2, iPsi = iPsi.2,
             absdPsi=absdPsi.2, absdiPsi=absdiPsi.1, log = TRUE),
      ylab = "dDiagA(x, theta= 38, *, log=TRUE)",
      0, 1, col = 4, n = 1000)
## previously
curve(dDiagA(x, th = 38, d = 2, iPsi = iPsi.1,
             absdPsi=absdPsi.2, absdiPsi=absdiPsi.1, log = TRUE),
      add = TRUE, col = "darkgray", lwd=2, lty=3, n=1000)


###################################################
### code chunk number 10: ex-U-data
###################################################
d <- 5
(theta <- copFrank@iTau(tau = 0.75))
cop <- onacopulaL("Frank", list(theta, 1:d))
set.seed(1); for(l in 1:4) U <- rnacopula(n = 100, cop)
U. <- sort(apply(U, 1, max)) # build the max


###################################################
### code chunk number 11: negLogL
###################################################
mlogL <- function(theta)
    -sum(dDiagA(U., theta, d=d, iPsi = iPsi.2,
                absdPsi=absdPsi.2, absdiPsi=absdiPsi.1,
                log = TRUE))


###################################################
### code chunk number 12: llog-theta-plot
###################################################
p.mlogL <- function(th, mlogL, col= "red2", lwd = 1, lty = 1,
                    add= FALSE) {
  stopifnot(is.numeric(th), is.function(mlogL))
  nll <- vapply(th, mlogL, 0.)
  if(add) lines(nll ~ th, col=col, lwd=lwd, lty=lty)
  else plot(nll ~ th, xlab=expression(theta),
            ylab = expression(- logLik(theta, .)),
            type = "l", col=col, lwd=lwd, lty=lty)
  invisible(nll) # return invisibly
}
thet <- seq(11, 99, by = 1/4)
p.mlogL(thet, mlogL)

require("Rmpfr")## compute the same with *high* accuracy ...
## using three different precisions:
MPrecBits <- c(160, 128, 96)
mkNm <- function(bits) sprintf("%03d.bits", bits)
## As it takes a while, cache the result:
fnam <- sprintf("mlogL_mpfr_%s.rda", Sys.info()[["machine"]])
if (!file.exists(fn <- file.path(copDDir,fnam))) {
  print(system.time(
    nllMP <- lapply(MPrecBits, function(pBit) {
        nlM <- thM <- mpfr(thet, precBits = pBit)
        ## (vapply() does not work for "Rmpfr":)
        for(i in seq_along(thet)) nlM[i] <- mlogL(thM[i])
        nlM
    })
  )) ## 91.226 0.013 91.506 [nb-mm icore 5]
  names(nllMP) <- mkNm(MPrecBits)
  copSrcDDir <- if(Sys.getenv("USER") == "maechler")
      '~/R/Pkgs/copula/inst/doc' else ""
  if(file.exists(copSrcDDir))# <<- only for certain users; not on CRAN etc
      save(nllMP, file = file.path(copSrcDDir, fnam))
} else load(fn)

colB <- c("blue3","violetred4","tan3")
ltyB <- c(5:3)
lwdB <- c(2,2,2)

for(i in seq_along(nllMP)) {
    lines(thet, as.numeric(nllMP[[i]]),
          col=colB[i], lty = ltyB[i], lwd = lwdB[i])
}
leg <- c("double prec.", sprintf("mpfr(*, precBits = %d)", MPrecBits))
legend("top", leg,
       col= c("red3",colB), lty=c(1, ltyB), lwd=c(1,lwdB), bty="n")


###################################################
### code chunk number 13: empty-prompt
###################################################
op <- options(prompt = " ")


###################################################
### code chunk number 14: mlogL_2terms (eval = FALSE)
###################################################
##    absdPsi(d*iPsi(u,th), th, log = TRUE) +
##    absdiPsi(u,          th, log = TRUE)


###################################################
### code chunk number 15: 3f (eval = FALSE)
###################################################
## iPsi = iPsi.2; absdPsi = absdPsi.2; absdiPsi = absdiPsi.1


###################################################
### code chunk number 16: reset-prompt
###################################################
options(op)


###################################################
### code chunk number 17: check-iPsi
###################################################
stopifnot(all.equal(iPsi.2(U.,  50 ),
                    iPsi.2(U., mpfr(50, 96))),
          all.equal(iPsi.0(U., mpfr(50, 200)),
            pI.U <- iPsi.2(U., mpfr(50, 200)), tol=1e-50) )


###################################################
### code chunk number 18: plot-absdPsi-2
###################################################
psD.n <-            absdPsi.2(as.numeric(pI.U), 40)
psD.M <- as.numeric(absdPsi.2(pI.U,        mpfr(40, 200)))
matplot(U., cbind(psD.n, psD.M), type="l", log="y")
legend("top", c("double prec.", "mpfr(*, precBits = 200)"),
       col= 1:2, lty=1:2, bty="n")


###################################################
### code chunk number 19: plot-absdPsi-zero
###################################################
u0 <- 2^-(100:1)
psD.n <-            absdPsi.2(u0, 40)
psD.M <- as.numeric(absdPsi.2(u0, mpfr(40, 200)))
matplot(u0, cbind(psD.n, psD.M), type="l", log="xy")
legend("top", c("double prec.", "mpfr(*, precBits = 200)"),
       col= 1:2, lty=1:2, bty="n")


###################################################
### code chunk number 20: log1mexp
###################################################
log1mexp <- function(a) # accurately compute log(1-exp(-a))
{
    stopifnot(a >= 0)
    r <- a
    tst <- a <= log(2)
    r[ tst] <- log(-expm1(-a[ tst]))
    r[!tst] <- log1p(-exp(-a[!tst]))
    r
}


###################################################
### code chunk number 21: absdPsi.3
###################################################
absdPsi.3 <- function(t, theta, log=FALSE) {
  w <- log1mexp(theta) - t
  Li. <- if(log) w - log1mexp(-w) else -exp(w)/expm1(w)
  if(log) Li. - log(theta) else Li. / theta
}


###################################################
### code chunk number 22: nlogL-2-plot
###################################################
p.mlogL(th = seq(11, 99, by = 1/4),
        mlogL = (mlogL2 <- function(theta)
                 -sum(dDiagA(U., theta, d=d, iPsi = iPsi.2,
                            absdPsi = absdPsi.3, absdiPsi=absdiPsi.1,
                            log = TRUE))), lwd = 2)


###################################################
### code chunk number 23: llog-theta-plot-2
###################################################
thet <- 9:1000
nll <- p.mlogL(thet, mlogL = mlogL2, lwd=2)
(th0 <- thet[i0 <- max(which(is.finite(nll)))])
abline(v = th0, col="red2", lty="15", lwd=2)


###################################################
### code chunk number 24: dDiag-large-thet
###################################################
dDiagA(0.999, 715, d = d, iPsi = iPsi.2, absdPsi = absdPsi.3,
                          absdiPsi = absdiPsi.1, log = TRUE)


###################################################
### code chunk number 25: psiID1-large-thet
###################################################
absdiPsi.1(0.999, th = 715, log=TRUE)


###################################################
### code chunk number 26: def-absdiPsi-2
###################################################
absdiPsi.2 <- function(u, theta, log = FALSE)
  if(log) log(theta)- {y <- u*theta; y + log1mexp(y)} else theta/expm1(u*theta)


###################################################
### code chunk number 27: llog-theta-plot-3
###################################################
plot(nll ~ thet, xlab=expression(theta),
     ylab = expression(- logLik(theta, .)),
     type = "l", col="red2", lwd=2)
abline(v = th0, col="red2", lty="15", lwd=2)
nll3 <- p.mlogL(thet, mlogL = function(theta)
                 -sum(dDiagA(U., theta, d=d, iPsi= iPsi.2, absdPsi= absdPsi.3,
                             absdiPsi = absdiPsi.2, log = TRUE)),
                col = "blue2", lwd=3, lty=2, add = TRUE)
nll3[thet == 800]


###################################################
### code chunk number 28: true-iPsi-small
###################################################
pI <- iPsi.2(u=0.999, th= mpfr(800, 200))
cat(sapply(list(pI, .Machine$double.xmin),
           format, digits = 7), "\n")


###################################################
### code chunk number 29: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 30: finalizing
###################################################
options(op.orig)


