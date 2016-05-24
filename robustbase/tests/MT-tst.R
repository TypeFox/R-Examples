require("robustbase")

##---> ./poisson-ex.R
##     ~~~~~~~~~~~~~~  for more glmrobMT() tests

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...
source(system.file("xtraR/ex-funs.R", package = "robustbase"))
## -> newer assert.EQ()  {TODO: no longer needed in 2015}

if(!require("sfsmisc")) {
    eaxis <- axis  # so we can use  eaxis() below
}


(doExtras <- robustbase:::doExtras())

## Explore the espRho() function: ---------------------------------------------
if(!dev.interactive(orNone=TRUE)) pdf("MT-E_rho.pdf")
E.rho <- robustbase:::espRho
lambdas <- ((1:10)/2)^2
cws <- c(1, 1.5, 1.75, 2, 2.25, 3)
(gr <- expand.grid(lam = lambdas, cw = cws))

Egr <- apply(gr, 1, function(r) {
    lam <- r[["lam"]]; cw <- r[["cw"]]; sL <- sqrt(lam)
    xx <- seq(lam - 2*sL, lam + 2*sL, length=17)
    vapply(xx, function(X) E.rho(lam, xx=X, cw=cw), NA_real_)
})
str(Egr)# 17 x 60
mLeg <- function(pos, type="o")
    legend(pos, legend=paste("lambda = ", format(lambdas, digits=2)),
           lty=1:5, col=1:6, pch= c(1:9, 0, letters, LETTERS), bty="n")
matplot(Egr[, gr[,"cw"]== 1.0 ], type="o",main="c_w = 1.0" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"]== 1.5 ], type="o",main="c_w = 1.5" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"]== 1.75], type="o",main="c_w = 1.75"); mLeg("bottomright")
matplot(Egr[, gr[,"cw"]== 2.0 ], type="o",main="c_w = 2.0" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"]== 2.25], type="o",main="c_w = 2.25"); mLeg("bottomright")
matplot(Egr[, gr[,"cw"]== 3.0 ], type="o",main="c_w = 3.0" ); mLeg("bottomright")

dev.off()


## Explore the m() function: ---------------------------------------------
if(!dev.interactive(orNone=TRUE)) pdf("MT-m_rho.pdf")

mkM <- robustbase:::mk.m_rho # itself calling splinefun(*, "monoH.FC")
getSpline.xy <- function(splfun) {
    ## Depending on the version of R, the
    ## environment of splinefun() slightly changes:
    stopifnot(is.function(splfun), length(e <- environment(splfun)) > 0)
    if("x0" %in% ls(e))
	list(x = e$x0, y = e$y0)
    else list(x = e$x, y = e$y)
}

m21 <- mkM(2.1, recompute=TRUE)# the default 'cw = 2.1'
m16 <- mkM(1.6, recompute=TRUE)
p.m2 <- function(mrho, from = 0, to, col=2, addKnots=TRUE, pchK=4, cexK=1.5, ...) {
    stopifnot(is.function(mrho))
    curve(mrho, from, to, col=col, ...)
    curve(sqrt(x), add=TRUE, col=adjustcolor("gray",.5), lwd=2)
    if(addKnots) points(getSpline.xy(mrho), pch=pchK, cex=cexK)
}
p.m.diff <- function(mrho, from = 0, to, col=2, addKnots=TRUE, pchK=4, cexK=1.5, ...) {
    stopifnot(is.function(mrho))
    curve(mrho(x) - sqrt(x), from=from, to=to, n=512, col=col, ...)
    abline(h=0,lty=3)
    if(addKnots) {
	xy <- getSpline.xy(mrho)
	if(is.numeric(x <- xy$x))
	    points(x, xy$y - sqrt(x), pch=pchK, cex=cexK)
        else warning("'addKnots' not available: No knots in function's environment")
    }
}

p.m2(m21, to=10)
p.m2(m16, to=10)
p.m2(m21, to=50)
p.m2(m21, to=120, cexK=.8)
p.m.diff(m21, to=120, cex=.5)# pchK="."
p.m.diff(m16, to=120, cex=.5)# pchK="."

mm21 <- function(.) robustbase:::mm(., m21)
environment(mm21) <- environment(m21)# <- for p.m()
p.m2(mm21, to=120, cexK=.8)
p.m.diff(mm21, to=120, cexK=.8)#-- discontinuity at 100 !!
## TODO: ways to improve!

## Here: look at "larger lambda" (and more cw)

la2 <- 5*2^seq(0, 10, by = 0.25)
c.s <- .25*c(1:10, 15, 50)
mL <- lapply(c.s, function(cc) mkM(cc, lambda = la2, recompute=TRUE))
str(mL, max=1) # a list of functions..
assert.EQ(la2, getSpline.xy(mL[[1]])$x)
mmL <- sapply(mL, function(F) getSpline.xy(F)$y)
matplot(la2, mmL, type ="l") # "all the same" from very far ...
mm.d. <- mmL - sqrt(la2)
matplot(la2, mm.d., type ="l", xlab=quote(lambda)); abline(h=0, lty=3)
legend("bottom", legend= paste("cw=",c.s), col=1:6, lty=1:5, ncol = 3, bty="n")

matplot(la2, -mm.d., type ="l", xlab=quote(lambda), log = "xy", axes=FALSE)
eaxis(1); eaxis(2)
legend("bottom", legend= paste("cw=",c.s), col=1:6, lty=1:5, ncol = 3, bty="n")
## ok, that's the correct scale
c.s2  <- c.s  [c.s >= .75]
mm.d2 <- mm.d.[, c.s >= .75]

matplot(la2, -mm.d2, type ="l", xlab=quote(lambda), log = "xy", axes=FALSE)
eaxis(1); eaxis(2)
legend("bottomleft", legend= paste("cw=",c.s2), col=1:6, lty=1:5, ncol = 3, bty="n")

##->   log (sqrt(lam) - m(lam)) = a[c] - beta * log(lam) :
dd2 <- data.frame(m.d = c(mm.d2),
                  cw = rep(c.s2, each = length(la2)),
                  lambda = rep(la2, length(c.s2)))

## gives a pretty nice picture:
summary(fm <- lm(log(-m.d) ~ 0+factor(cw) + log(lambda),
                 data = dd2, subset = lambda >= 50))
##=> slope of log(lambda) = -1/2
dd3 <- within(dd2, { ld2 <- log(-m.d) + 1/2 * log(lambda) })[dd2[,"lambda"] >= 50,]
plot(ld2 ~ cw, data = dd3, type = "b")
plot(ld2 ~ cw, data = dd3, type = "b", log="x")
coplot(ld2 ~ cw|lambda, data = dd3)
coplot(ld2 ~ cw|log(lambda), data = dd3)
coplot(ld2 ~ log10(cw) | log10(lambda), data = dd3)

dev.off()
##-------------------------------------------------------- end m(.) -------------


## The simple intercept example from  ./glmrob-1.R
set.seed(113)
y <- rpois(17, lambda = 4)
y[1:2] <- 99:100 # outliers
y.1 <- y
x.1 <- cbind(rep(1, length(y)))

options("robustbase:m_rho_recompute" = TRUE)#-> recompute in any case:
showSys.time( r <- glmrob(y ~ 1, family = poisson, method = "MT", nsubm=100) )# some output
str(r)

## was   c(ini = 1.30833281965018, est = 1.29369680430613)
## then  c(ini = 1.30833281965018, est = 1.29369680422799)
##       c(ini = 1.30833281965018, est = 1.29369680430627)
r.64b <- c(ini = 1.30833281965018, est = 1.29369680452016)
stopifnot(r$converged)
assert.EQ(r$initial,      r.64b[["ini"]], check.attributes=FALSE, tol = 1e-13)# rel.diff: 3.394.e-16
assert.EQ(r$coefficients, r.64b[["est"]], check.attributes=FALSE, tol = 1e-09)# as long we use different optim())


## now, as the algorithm has a random start:
set.seed(7)
nSim <- if(doExtras) 20 else 2
showSys.time(LL <- replicate(nSim,
     glmrob(y ~ 1, family = poisson, method = "MT"),
                             simplify=FALSE))
ini <- sapply(LL, `[[`, "initial")
est <- sapply(LL, `[[`, "coefficients")
## surprise:  all the 20 initial estimators are identical:
stopifnot(diff(range(ini)) == 0,
          diff(range(est)) == 0)
## probably too accurate ... but ok, for now
assert.EQ(est[1], r.64b[["est"]], check.attributes=FALSE, tol = 1e-10)# Winbuilder needed ~ 2e-11
assert.EQ(ini[1], r.64b[["ini"]], check.attributes=FALSE, tol = 1e-10)

ccvv <- sapply(LL, `[[`, "cov")
stopifnot(ccvv[1] == ccvv)
assert.EQ(print(ccvv[1]), 0.0145309081924157, tol = 1e-7, giveRE=TRUE)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
## "Platform" info
(SysI <- Sys.info()[c("sysname", "release", "nodename", "machine")])
if(require("sfsmisc") && SysI[["sysname"]] == "Linux")
    ## not on the Mac (yet)
    c(SysI, MIPS=Sys.MIPS(), Sys.sizes()) else SysI
