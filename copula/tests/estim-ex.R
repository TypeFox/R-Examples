## Copyright (C) 2012, 2015 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

###__ Estimation methods for Archimedean Copulas __

require(copula)

(doExtras <- copula:::doExtras())

## From source(system.file("test-tools.R", package = "Matrix")) :
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

### We have other implicit estimation "tests" in  demo(estimation.gof)
## i.e., ../demo/estimation.gof.R
##       ~~~~~~~~~~~~~~~~~~~~~~~~

## Frank & mle ---- log-density is *REALLY* not good enough:
## -----------
## at least not for "GoF" where we may have "tail-dependent"
## observations which are highly improbable under Frank and would
## "need" a very large theta :

p.log.f <- function(t0, t1, n.th = 256,
                    thet.vec = seq(t0, t1, length.out= n.th),
                    d, u0, uu = rbind(rep(u0, length.out = d)),
                    cop = "Frank", legend.x = "bottomright",
                    col = mapply(adjustcolor, col=c(1,3,2,4,1,2), alpha=c(2, 5:7, 4,2)/10),
                    lwd = c(5,4,3,1,5,5), lty = 1:6)
{
    stopifnot(d == as.integer(d), length(d) == 1, d >= 2,
              NCOL(uu) == d, is.numeric(thet.vec), length(thet.vec) > 1)
    cop <- onacopulaL(getAcop(cop), list(NA, 1:d))

    stopifnot(is.function(.f. <- cop@copula@dacopula))
    ## meths <- c("negI-s-Stirling", "negI-s-Eulerian")
    meths <- local({m <- eval(formals(polylog)$method);m[grep("^negI-s", m)]})
    l.arr <- sapply(meths, function(meth)
                {
                    sapply(c(FALSE,TRUE), function(Li.log)
                           sapply(thet.vec, .f., u = uu, log=TRUE,
                                  method = meth, Li.log.arg=Li.log))
                }, simplify = "array")
    dim(l.arr) <- local({D <- dim(l.arr); c(D[1], prod(D[-1]))})
    ## n x 4
    nms <- c(t(outer(substring(meths, 8),
                     paste("Li_n(", c("<reg.>", "log(arg)"), ")", sep=""),
                     ##paste("Li.log.A=",c(FALSE,TRUE),sep=""),
                     paste, sep="_")))
    colnames(l.arr) <- nms

    tit <- bquote(.(cop@copula@name)*" copula -  log density " ~~
                  log~f(theta, bold(u)[d==.(d)]))
    matplot(thet.vec, l.arr, type = "l", lwd=lwd, col=col, lty=lty,
            main = tit, xlab=expression(theta), ylab = expression(log~f(theta, .)))
    legend(legend.x, nms, col=col, lwd=lwd, lty=lty,
           title = "polylog(.., method = \"*\")", inset = .01)
    mkU <- function(u, n1 = 4, n2 = 2) {
        n <- length(u <- c(u))
        cu <- if(n > (n1+n2)) ## use "..."
            c(format(head(u, n1)), "...", format(tail(u, n2))) else format(u)
        paste("(", paste(cu, collapse = ", "), ")", sep="")
    }
    mtext(bquote(bold(u)[d==.(d)] == .(mkU(uu))), line=1/8)
    invisible(structure(cbind(theta = thet.vec, l.arr),
			u = uu, dacopula = .f.))
}

unattr <- function(obj) { mostattributes(obj) <- NULL ; obj }
show. <- function(pmat) unattr(pmat)[, c(1,3,5,7)]

if(!dev.interactive(orNone=TRUE)) pdf("estim-ex.pdf")

th <- seq(1,60, by=1/4)
show.(m1 <- p.log.f(thet.vec = th, d = 5, u0 = .987))
## wow! asymptotic is quite good even here!

## using  MC  (instead of polylog(.))  also ``works'' {the function is noisy ..}:
l.th.MC <- sapply(th, attr(m1, "dacopula"), u = attr(m1, "u"),
                  n.MC = 2000, log=TRUE) ## 2000: be speedy in test
lines(th, l.th.MC, lwd=3, col=adjustcolor("sandybrown", .5))
legend("right", "n.MC = 2000", lwd=3, col=adjustcolor("sandybrown", .5),
       inset=.01)
rm(th)

showProc.time()

## Extend the range:
show.(m2 <- p.log.f(1, 200, n.th=401, d = 5, u0 = .987))

if(doExtras) {
## Extend the range even more -- change u0
show.(m3 <- p.log.f(10, 500, d = 5, u0 = .96))
## and more
show.(m3.2 <- p.log.f(10, 800, d = 5, u0 = .96))#-> breakdown at ~ 775..
}
showProc.time()


## higher d:
if(doExtras) {
show.(m4 <- p.log.f(10, 500, d = 12, u0 = 0.95))
m5 <- p.log.f(10, 500, d = 12, u0 = 0.92)
m6 <- p.log.f(1,  200, d = 12, u0 = c(0.8, 0.9), legend.x="topright")
m7 <- p.log.f(1,  400, d = 12, u0 = c(0.88, 0.9))
}

## whereas this now also overflows for Eulerian *AND* asymp. + log.arg:
show.(mm <- p.log.f(10, 1000, d = 12, u0 = 0.9))

##--> investigation shows that  w is *also* underflowing to 0  (!!!) :

myTracer <- quote({
    cat(sprintf("Tracing %s: variables at entry are\n",
                ## the  -4   is purely by trial and error -- better?
                deparse(sys.call(-4)[1])))
    print(ls.str())
})
trace(polylog, tracer = myTracer, print=FALSE)
## Tracing function "polylog" in package "copula"

f <- attr(mm,"dacopula")

f(attr(mm,"u"), 825, log=TRUE, method="negI-s-asymp")
## ....
## z : num -4.15e-322
## [1] 61.48259
f(attr(mm,"u"), 800, log=TRUE, method="negI-s-asymp")
## z : num -2.44e-312

if(FALSE)
    debug(f)
## and then realize that in f(), i.e.,  Frank's dacopula() final line,
##     (d-1)*log(theta) + Li. - theta*u.sum - lu
## the terms     Li. - (theta * u.sum)
## 'more or less cancel' -- so if we could compute a
##  *different* rescaled polylog() .. maybe we can get better

untrace(polylog)

if(doExtras)
## and similarly here:
ll <- p.log.f(1, 12000, d = 12, u0 = 0.08)

showProc.time()

###------- Diagonal MLE --- dDiag() , etc ----------------------------------

## Some basic  dDiag() checks  {that failed earlier}:
stopifnot(identical(0, dDiag(0, onacopulaL("AMH", list(0.5, 1:4)))),
          all.equal(dDiag(.9, onacopulaL("Frank", list(44, 1:3))),
                    1.008252, tolerance= 1e-5),
          all.equal(dDiag(.9, onacopulaL("Joe",   list(44, 1:3))),
                    1.025283, tolerance= 1e-5),
          TRUE)

demo("dDiag-plots", package = "copula")
##    -----------> ../demo/dDiag-plots.R

showProc.time()

## A problem to solve:
curve(copAMH@dDiag(exp(-x), theta=.9, d=4, log=TRUE), 0, 500, n=2001,
      col = 2, lwd= 2,
      main = "FIXME (*not* urgent): dDiag(u, <AMH>,  log=TRUE)   for small u")

showProc.time()

##----------- Test all estimation methods ----------------

## all methods
(estMeth <- eval(formals(enacopula)$method))

## Use selected (the best 6) estimation methods:
## (estMeth <- c("mle", "smle", "dmle", "tau.tau.mean", "mde.chisq.CvM", "mde.chisq.KS"))

n <- if(doExtras) 128 else 8
d <- 6 ; theta <- 2
(cop <- onacopula("Clayton", C(theta, 1:d)))
set.seed(1)
U <- rnacopula(n, cop)
## Show that 'racopula()' is indeed identical {and faster for the case of non-nested}
stopifnot(identical(U,
                {set.seed(1); copula:::racopula(n, copClayton, theta=theta, d=d)}))

rr <- sapply(estMeth, function(e) {
    enacopula(U, cop, e, n.MC = if(e == "smle") 1000 else 0)
})

round(cbind(rr, bias = rr - theta), 3)

showProc.time()

### 2-level nested Archimedean copulas :
if(doExtras) {
    ## was simply (!) :   demo("dnac-demo")
    v.dnac <- vignette("dNAC", package="copula")
    R.dnac <-
	if(nzchar(v.dnac$R))
	    with(v.dnac, file.path(Dir, "doc", R))
	else { ## default: knitr::rmarkdown() unfortunately does *not* provide R
	    stopifnot(file.exists(ff <- with(v.dnac, file.path(Dir, "doc", File))))
	    oR <- file.path(tempdir(), sub("Rmd$", "R", basename(ff)))
	    knitr::purl(ff, output = oR)
	    oR
	}
    cat("R file: ", R.dnac, "\n -- now source()ing it :\n")
    source(R.dnac, echo = TRUE, max.deparse.length = Inf,
           keep.source = TRUE, encoding = getOption("encoding"))
}


showProc.time()
