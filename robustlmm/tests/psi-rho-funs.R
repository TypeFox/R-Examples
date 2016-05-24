require(robustlmm)

## initialize device
pdf("psi-rho-funs.pdf")

## Test E... slots
t.fun <- function(obj) {
    rho <- obj@rho
    psi <- obj@psi
    Dpsi <- obj@Dpsi
    c(Erho=
      all.equal(integrate(function(x) rho(x)*dnorm(x), -Inf, Inf,
                          rel.tol = .Machine$double.eps^0.5)$value,
                obj@Erho()),
      Epsi2=
      all.equal(integrate(function(x) psi(x)^2*dnorm(x), -Inf, Inf,
                          rel.tol = .Machine$double.eps^0.5)$value,
                obj@Epsi2()),
      EDpsi=
      all.equal(integrate(function(x) Dpsi(x)*dnorm(x), -Inf, Inf,
                          rel.tol = .Machine$double.eps^0.5)$value,
                obj@EDpsi()))
}

stopifnot(t.fun(huberPsi))

t.plot.slot <- function(slot) {
    curve(slot(huberPsi, slot)(x), -5, 5,)
    curve(slot(smoothPsi, slot)(x), -5, 5, add = TRUE)
}

t.plot.slot("rho")
t.plot.slot("psi")
t.plot.slot("Dpsi")
t.plot.slot("wgt")
t.plot.slot("Dwgt")

stopifnot(t.fun(smoothPsi))

## TODO: add tests for derivatives like in robustbases lmrob-psifns.R

p.psiFun <- function(x, object,
                     col = c("black", "red3", "blue3", "dark green"),
                     leg.loc = "right", ...)
{
    ## Author: Martin Maechler, Date: 13 Aug 2010, 10:17
    m.psi <- cbind(rho    = object@rho(x),
                   psi    = object@psi(x),
                   "psi'" = object@Dpsi(x),
                   wgt    = object@wgt(x),
                   "wgt'" = object@Dwgt(x))
    fExprs <- quote(list(rho(x), psi(x), {psi*minute}(x), w(x) == psi(x)/x))
    matplot(x, m.psi, col=col, lty=1, type="l",
            main = substitute(FFF ~~ ~~ " with "~~ psi*"-type" == PSI(PPP),
                              list(FFF = fExprs, psi = object@name,
                                   PPP = paste(formals(object@rho)[-1], collapse=", "))),
            ylab = quote(f(x)), xlab = quote(x), ...)
    abline(h=0,v=0, lty=3, col="gray30")
    fE <- fExprs; fE[[1]] <- as.name("expression")
    legend(leg.loc, inset=.02, eval(fE), col=col, lty=1)
    invisible(cbind(x=x, m.psi))
}

mids <- function(x) (x[-1]+x[-length(x)])/2
chkPsiDeriv <- function(m.psi, tol = 1e-4) {
    stopifnot(length(tol) > 0, tol >= 0,
              is.numeric(psi <- m.psi[,"psi"]),
              is.numeric(dx  <- diff(x <- m.psi[,"x"])))
    if(length(tol) < 2) tol[2] <- 10*tol[1]
    xn0 <- abs(x) > 1e-5
    Dpsi <- diff(psi)/dx
    DnInf <- abs(Dpsi) < 50
    c(all.equal(diff(m.psi[,"rho"])/dx, mids(psi), tolerance=tol[1]), # rho'  == psi
      all.equal(Dpsi[DnInf], mids(m.psi[,"psi'"])[DnInf], tolerance=tol[2]),# psi'  == psip
      all.equal((psi/x)[xn0], m.psi[xn0,"wgt"], tolerance=tol[1]/10),# psi/x == wgt
      all.equal(m.psi[,"wgt'"],
                ifelse(x==0,0,m.psi[,"psi'"]/x - m.psi[,"psi"]/(x*x)),
                tolerance=tol[1])) # wgt' == psi'/x - psi/x^2
}

head(x. <- seq(-5, 10, length=1501))
## [separate lines, for interactive "play": ]
stopifnot(chkPsiDeriv(p.psiFun(x., huberPsi)))
head(x. <- seq(-5, 10, length=2501))
chkPsiDeriv(p.psiFun(x., smoothPsi))
stopifnot(chkPsiDeriv(p.psiFun(x., smoothPsi)))
head(x. <- seq(-5, 10, length=1501))

## intPsi <- Vectorize(function(x) integrate(function(y) smoothPsi@psi(y), 0, x)$value)
## curve(intPsi, 0, 10)
## curve(smoothPsi@rho(x), 0, 10, add=TRUE)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
