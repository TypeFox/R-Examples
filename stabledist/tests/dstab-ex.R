require("stabledist")
dPareto <- stabledist:::dPareto

source(system.file("test-tools-1.R", package = "Matrix"), keep.source=interactive())
					#-> identical3(), showProc.time(),...
(doExtras <- stabledist:::doExtras())
if(!require("sfsmisc")) eaxis <- axis # use sfsmisc::eaxis if available

stopifnot(0 < print(dstable(4000., alpha=1.00001, beta=0.6)))
## gave error in fBasics::dstable()
## now 18 warnings from uniroot()

pdf("dstab-ex.pdf")

x <- 2^seq(0, 20, length= if(doExtras) 200 else 64)
## Regression check for alpha=2: <==> *norm() :
x. <- x/1024
fx <- dstable(x., alpha = 2, beta = 0, gamma = 1.1, delta=2, pm=2)
lf <- dstable(x., alpha = 2, beta = 0, gamma = 1.1, delta=2, pm=2, log=TRUE)
stopifnot(
    local({i <- is.finite(log(fx)); all.equal(log(fx[i]), lf[i])}),
    all.equal(fx, dnorm(x., m=2, s=1.1)),
    all.equal(lf, dnorm(x., m=2, s=1.1, log=TRUE)))

fx <- dstable(x, alpha = 1.0001, beta = 0.6)

plot(x,fx, log="x", type="l")# looks good
plot(x,fx, log="xy", type="l")# --- perfect now
stopifnot((dlx <- diff(log(fx))) < 0, # decreasing
      { dl <- dlx[x[-1] > 1000] # the negative slope:
	relErr(if(doExtras) -0.13934 else -0.44016, dl) < 4e-4},
	  diff(dl) > -1e-6)#  > 0 "in theory"; > -6e-7, ok on 64-bit Linux

nc <- if(doExtras) 512 else 101 # number of points for curve()

zeta <- function(alpha,beta) if(alpha==1) 0 else -beta*tan(pi/2*alpha)

## negative beta:
cx <- curve(dstable(x, 0.75, -.5), -.5, 1.5, n=nc)# ok, now
m <- stableMode(0.75, -.5, tol=1e-14)
stopifnot(all.equal(m, 0.35810298366, tol = 1e-7),
	  all.equal(dstable(m, 0.75, -.5), 0.3554664043, tol=1e-6))

showProc.time()

###-------- "small" alpha -----------------
## alpha --> 0 --- very heavy tailed -- and numerically challenging.

## symmetric (beta = 0)
(x0 <- (-16:16)/256)
fx0 <- dstable(x0, alpha = 0.1, beta=0, gamma = 1e6)
plot(x0, fx0, type = "o",
     main = expression(f(x, alpha== 0.1, beta == 0, gamma == 10^6)))
stopifnot(all.equal(fx0[17],1.15508291498374),
	  all.equal(fx0[ 1],0.02910420736536),
	  all.equal(range(diff(fx0[1:8])),
		    c(0.0011871409, 0.0025179435), tol=1e-6)
	  )

## beta > 0
r3 <- curve(dstable(x, alpha = 0.3, beta = 0.5, tol=1e-7),
	    -1, 1)
m3 <- stableMode(0.3, 0.5, tol=1e-14)# still with 3 warnings
stopifnot(all.equal(m3, -0.2505743952946, tol = 1e-10))
## zoom into the above
r3. <- curve(dstable(x, alpha = 0.3, beta = 0.5, tol=1e-7),
	    -.27, -.22)
abline(v = m3, col="gray40", lty=2)

r1 <- curve(dstable(x, alpha = 0.1, beta = 0.5, tol=1e-7),
	    -.4, .2, ylim = c(0, 20), n = nc)
m1 <- stableMode(0.1, 0.5, tol=1e-15)
abline(v=m1, h=0, col="gray40", lty=2)
stopifnot(all.equal(m1, -0.079193, tol=1e-5)) # -0.079193188, was -0.079192219, and -0.0791921758
title(main = expression(f(x, alpha== 0.1, beta == 0.5)))
## check mode *and* unimodality
i. <- r1$x > m1
stopifnot(## decreasing to the right:
	  diff(r1$y[ i.]) < 0,
	  ## increasing on the left:
	  diff(r1$y[!i.]) > 0)

## Now *really* small alpha --> 0:
##     --------------------------
## Note that
if(require("FMStable")) {
    try( FMStable::setParam(alpha = 1e-18, location=0, logscale=0, pm = 0) )
## gives
## Error in .... setParam: S0=M parametrization not suitable for small alpha
}
## now if this is true (and I think we can trust Geoff Robinson on this),
## we currently have a  "design bug - problem":
## as internally, we always translate to  'pm=0' parametrization __FIXME_(how ??)__
## --> solution: see below:  there's a simple  (alpha -> 0) asymptotic

## These now "work": .... well with integration warnings !!

pdstab.alpha <- function(x, beta, alphas = 2^-(40:2),
                         type = "o", add=FALSE, log = "xy", ...)
{
    stopifnot(is.numeric(x), length(x) == 1, length(beta) == 1,
              is.numeric(beta), -1 <= beta, beta <= 1,
              is.numeric(alphas), 0 <= alphas, alphas <= 2,
              is.logical(add))
    if(is.unsorted(alphas)) alphas <- sort(alphas)
    dst.alp <- vapply(alphas, function(aaa)
                    dstable(x, aaa, beta = beta, pm = 0), 1.) ## with warnings
    if(add)
        lines(alphas, dst.alp, type=type, ...)
    else {
        plot(alphas, dst.alp, type=type, log=log, axes=FALSE,
             xlab = quote(alpha),
             ylab = bquote(f(x == .(x), alpha)),
             main = bquote(dstable(x == .(x), alpha, beta == .(beta), pm == 0)
             ~~~"for (very) small" ~ alpha))
        if(!require("sfsmisc")) eaxis <- axis # use sfsmisc::eaxis if available
        eaxis(1); eaxis(2)
    }
    invisible(cbind(alphas, dstable = dst.alp))
}

## vary beta,  keep x :
pdstab.alpha(x = -1, beta = 0.1)
if(doExtras) {
pdstab.alpha(x = -1, beta = 0.3, add=TRUE, col=2, type="l")
pdstab.alpha(x = -1, beta = 0.5, add=TRUE, col=3, type="l")
pdstab.alpha(x = -1, beta = 0.7, add=TRUE, col=4, type="l")
pdstab.alpha(x = -1, beta = 0.9, add=TRUE, col=5, type="l")

## vary x,  keep beta
pdstab.alpha(x =  -1, beta = 0.3)
pdstab.alpha(x =  +1, beta = 0.3, add=TRUE, col=2, type="l")
pdstab.alpha(x =  +5, beta = 0.3, add=TRUE, col=3, type="l")
pdstab.alpha(x = +50, beta = 0.3, add=TRUE, col=4, type="l")
pdstab.alpha(x = -10, beta = 0.3, add=TRUE, col=5, type="l")
}
## Plots suggest a simple formula
##  log(f(x, alpha, beta))=  c(x,beta) + log(alpha)
##      f(x, alpha, beta) =  C(x,beta) * (alpha)   -- ???

## for normal alpha, it looks different {which is reassuring!}
pdstab.alpha(x = -1, beta = 0.3, alphas = seq(1/128, 2, length=100))

## Show the formula, we derived:
## f(x, \alpha, \beta) ~  \alpha * \frac{1 + \beta}{2e \abs{x + \pi/2 \alpha\beta}}
## ONLY ok, when  x > zeta := - beta * tan(pi/2 *alpha)
## otherwise, "swap sign" of (x, beta)
dst.smlA <- function(x, alpha, beta, log = FALSE) {
    pa <- pi/2*alpha
    i. <- x < -pa*beta
    if(any(i.)) {
        beta <- rep(beta, length.out=length(x))
        beta[i.] <- -beta[i.]
        x   [i.] <- -x   [i.]
    }
    ## alpha*(1+beta) / (2*exp(1)*(x+ pa*beta))
    r <- log(alpha) + log1p(beta) - (1 + log(2*(x+ pa*beta)))
    if(log) r else exp(r)
}

set.seed(17)

alpha <- 1e-15 ## must be larger than cutoff in .fct1() in ../R/dpqr-stable.R :
stopifnot(alpha > stabledist:::.alpha.small.dstable)
for(beta in runif(if(doExtras) 20 else 2, -1, 1)) {
 cat(sprintf("beta =%10g: ", beta))
 for(N in 1:(if(doExtras) 10 else 1)) {
     x <- 10*rnorm(100); cat(".")
     stopifnot(all.equal(dstable (x, alpha, beta),
                         dst.smlA(x, alpha, beta), tol = 1e-13))
 }; cat("\n")
}


curve( dstable (x, alpha=1e-10, beta=.8, log=TRUE) , -10, 10)
curve( dst.smlA(x, alpha=1e-10, beta=.8, log=TRUE), add=TRUE, col=2)

## "same" picture (self-similar !)
curve( dstable (x, alpha=1e-10, beta=.8, log=TRUE), -1, 1)
curve( dst.smlA(x, alpha=1e-10, beta=.8, log=TRUE), add=TRUE, col=2)



## Testing stableMode here:


### beta = 1 (extremely skewed)  and small alpha: ---------
##  --------
## Problem at *left* ("less problematic") tail, namely very close to where the
## support of the density becomes mathematically exactly zero :
##
## clear inaccuracy / bug --- maybe *seems* curable
##
curve(dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE), -40, 10)
##            ------
## --> warnings both from uniroot ("-Inf") and .integrate2()
## about equivalent to
curve(dstable(x, alpha= 0.1, beta=1, pm=1, log=TRUE), 1e-15, 4e4,
      log="x", xaxt="n"); eaxis(1)
## If we decrease  zeta.tol "to 0", we get better here:
curve(dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE, zeta.tol=1e-100), -40, 20)
## or here, ... but still not good enough
curve(dstable(exp(x), alpha= 0.1, beta=1, pm=1, log=TRUE, zeta.tol=1e-200), -45, 30)



showProc.time()

##------ NB: Pareto tail behavior -- see more in ./tails.R
##						   =======

## alpha ~= 1  ---- and x ~ zeta(a,b) -----------
## ==========
f1 <- dstable(6366.197,	 alpha= 1.00001, beta= .1)
f2 <- dstable(-50929.58, alpha= 1.00001, beta= -.8)
stopifnot(f1 > 0, f2 > 0)

## these all work (luck):
zet <- zeta(alpha= 1.00001, beta= -.8)# -50929.58
## here, we must have larger zeta.tol: = 5e-5 is fine; now using "automatic" default
r4 <- curve(dstable(zet+x, alpha= 1.00001, beta= -.8), -3, 7,
	    xlab = expression(zeta(alpha,beta) - x),
	    ylim=c(2.1, 2.3)*1e-10, n = nc)
cc <- "pink3"
abline(v=0, col=cc)
z.txt <- bquote(paste(x == zeta(.), phantom() == .(signif(zet,6))))
mtext(z.txt, at=0, line = -1, adj = -.1, padj=.2, col=cc)
## no longer much noise (thanks to zeta.tol = 1e-5):
curve(dPareto(zet+x, alpha= 1.00001, beta= -.8), add=TRUE, col=2)
stopifnot({ rr <- range(r4$y)
	    2.1e-10 < rr & rr < 2.3e-10 })

showProc.time()

### ---- alpha == 1 ---------
curve(dstable(x, alpha = 1, beta = 0.3), -20, 20, log="y", n=nc)
curve(dstable(x, alpha = 1, beta = 0.3, log=TRUE), -200, 160, n=nc)

curve(dPareto(x, alpha = 1, beta = 0.3, log=TRUE), add=TRUE, col=4)
## "works", but discontinuous --- FIXME
## ditto:
curve(dstable(x, alpha=1, beta= 0.1, log=TRUE), -70,80, col=2)
curve(dPareto(x, alpha=1, beta= 0.1, log=TRUE), add=TRUE)

showProc.time()

dstable(-44, alpha=1, beta= .1)# failed
## large x gave problems at times:
dstable(-1e20, alpha = 0.9,  beta = 0.8)

chkUnimodal <- function(x) {
    ## x = c(x1, x2)  and  x1 is *increasing*  and x2 is *decreasing*
    stopifnot((n <- length(x)) %% 2 == 0,
	      (n2 <- n %/% 2) >= 2)
    if(is.unsorted(x[seq_len(n2)])) stop("first part is *not* increasing")
    if(is.unsorted(x[n:(n2+1)]))    stop("seconds part is *not* decreasing")
    invisible(x)
}

showProc.time()

xLrg <- c(10^c(10:100,120, 150, 200, 300), Inf)
xLrg <- sort(c(-xLrg, xLrg))
d <- dstable(xLrg, alpha = 1.8,	  beta = 0.3 ); chkUnimodal(d)
d <- dstable(xLrg, alpha = 1.01,  beta = 0.3 ); chkUnimodal(d) # (was slow!)
## look at the problem:
dstCurve <- function(alpha, beta, log=TRUE, NEG=FALSE,
		     from, to, n=nc, cLog=NULL, ...)
{
    if(NEG) {
	r <- curve(dstable(-x, alpha=alpha, beta=beta, log=log),
		   from=from, to=to, n=n, log = cLog, ...)
	curve(dPareto(-x, alpha=alpha, beta=beta, log=log), add=TRUE,
	      col=2, lwd=2, lty=2)
    } else {
	r <- curve(dstable(x, alpha=alpha, beta=beta, log=log),
		   from=from, to=to, n=n, log = cLog, ...)
	curve(dPareto(x, alpha=alpha, beta=beta, log=log), add=TRUE,
	      col=2, lwd=2, lty=2)
    }
    leg.ab <- paste0("(", if(NEG) "-x" else "x",
		     ", a=",formatC(alpha, digits=3),
		     ", b=",formatC(beta, digits=3),")")
    legend("topright", paste0(c("dstable ", "dPareto"), leg.ab),
	   col=1:2, lty=1:2, lwd=1:2, bty="n")
    invisible(r)
}

## (was *S.L.O.W.* on [2010-03-28] !)
r <- dstCurve(alpha = 1.01, beta = 0.3, NEG=TRUE,
	      from=1e10, to=1e20, cLog="x", ylim = c(-100, -45))
## zoom in:
r <- dstCurve(alpha = 1.01, beta = 0.3, , , .1e13, 9e13, ylim = c(-80, -55))
showProc.time()

d <- dstable(xLrg, alpha = 1.001, beta = -0.9) # >= 50 warnings
try( chkUnimodal(d) ) # FIXME
## look at the problem:
dstCurve(alpha = 1.001, beta = -0.9, log=FALSE, NEG=TRUE,  1e10, 1e20, cLog="xy")
## and at the right tail, too:
dstCurve(alpha = 1.001, beta = -0.9, log=FALSE, NEG=FALSE, 1000, 1e17, cLog="xy")

d <- dstable(xLrg, alpha = 1. ,	  beta = 0.3 ); chkUnimodal(d) # "ok" now
d <- dstable(xLrg, alpha = 0.9,	  beta = 0.3 ) # 10 warnings (had 11)
try( chkUnimodal(d) ) # FIXME
d <- dstable(xLrg, alpha = 0.5,	  beta = 0.3 ) # 19 warnings (had 22)
chkUnimodal(d)
d <- dstable(xLrg, alpha = 0.1,	  beta = 0.3 ) # 26 warnings (had 21)
chkUnimodal(d)

showProc.time()

##-------------	 beta = 1  ---------------------
options(dstable.debug = TRUE)
dstable(1, alpha=1.2,	beta= 1 - 1e-7)#ok
dstable(1, alpha=1.2,	beta= 1)# gave error, because	g(pi/2) < 0
dstable(0, alpha=13/16, beta= 1 -2^-52)# had error as	g(-theta0)  |->	 NaN
dstable(0, alpha=19/16, beta= 1)       # had error as	g(pi/2)	    |->	 NaN
options(dstable.debug = FALSE)

## NB: (beta=1, alpha = 1/2) is 'Levy' ---> dLevy() and some checks
## -- in ./pstab-ex.R
##	   ~~~~~~~~~~

if(doExtras) { ## actually "very-Extras" (checkLevel == "FULL")
 ## This needs 65 seconds (nb-mm3: 54*32*11 dstable() values)

 ep <- 2^-(1:54)## beta := 1 - ep ---> 1  {NB: 1 - 2^-54 == 1  numerically}
 alph.s <- (1:32)/16   # in (0, 2]
 f.b1 <- sapply(alph.s, function(alf)
		sapply(1 - ep, function(bet)
		       dstable(0:10, alpha = alf, beta = bet)),
		simplify = if(getRversion() >= "2.13") "array" else TRUE)
 print(summary(f.b1))
 r.b1 <- range(f.b1)
 stopifnot(0 < r.b1[1], r.b1[2] < 0.35)
 ## "FIXME" test more: monotonicity in x {mode is < 0 = min{x_i}}, beta, alpha, ...
 showProc.time()

} else message("expensive dstable() checks  have been skipped")

cat('Time elapsed: ', proc.time(),'\n') # "stats"
