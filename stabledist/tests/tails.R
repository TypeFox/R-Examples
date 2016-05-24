require("stabledist")

###--- Tail approximations etc  -- both for pstable() and dstable()
dPareto <- stabledist:::dPareto

source(system.file("test-tools-1.R", package = "Matrix"), keep.source=interactive())
					#-> identical3(), showProc.time(),...
(doExtras <- stabledist:::doExtras())
nc <- if(doExtras) 512 else 64 # number of points for curve()

pdf("stable-tails.pdf")

pstab.tailratio <- function(alpha, beta, n = nc, prob = 1/4096,
                            xmin = qstable(0.95, alpha,beta, tol = 0.01),
                            xmax = qstable(1 - prob, alpha,beta))
{
  ## Purpose: Explore eps(x) where  (1 - pstable(x))/(1- pPareto()) = 1 + eps(x)
  ##  <==>
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 21 Mar 2011, 10:09
    cl <- match.call()
    stopifnot(0 < prob, prob < .5,
              0 < xmin, xmin < xmax)
    x <- exp(seq(log(xmin), log(xmax), length.out = n))
    iF <-               pstable(x, alpha,beta, lower.tail=FALSE)
    ok <- iF > 0
    iF <- iF[ok]
    x <- x[ok]
    iFp <- stabledist:::pPareto(x, alpha,beta, lower.tail=FALSE)
    eps <- (iF - iFp)/iFp
    structure(list(x=x, eps=eps, call = cl, alpha=alpha, beta=beta),
              class = "pstableTailratio")
}

plot.pstableTailratio <- function(x, type="l", col="blue3",
                                  lin.col = adjustcolor("red2",.6), ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
    stopifnot(is.list(x), is.numeric(x$eps))
    tit <- substitute("tail ratio approx. for"~~ pstable(alpha==A, beta==B),
                      list(A=x[["alpha"]], B=x[["beta"]]))
    dat <- as.data.frame(x[c("x","eps")])
    ## NOTA BENE: Empirically,  I see that  eps > 0 <==> alpha > 1
    ##                                      eps < 0 <==> alpha < 1
    dat <- dat[dat[,"eps"] > 0, ] ## drop invalid  eps[]
    plot(eps ~ x, log = "xy", data = dat, col=col,
         ylab = expression(epsilon ~~~ "('eps')"),
         main = tit, type=type, ...)
    mtext( expression(epsilon(x) == (bar(F)(x,.) - bar(F)[P](x,.)) / bar(F)[P](x,.)) )
    fm <- lm(log(eps) ~ log(x), weights = x^2, data = dat)
    lines(dat[["x"]], exp(predict(fm)), col=lin.col)
    Form <- function(x) formatC(x, digits=4, wid=1)
    leg.line <-
        substitute(log(epsilon) == A + B * log(x),
                   list(A = Form(coef(fm)[["(Intercept)"]]),
                        B = Form(coef(fm)[["log(x)"]])))
    legend("topright", legend=leg.line, bty = "n", lty=1, col=lin.col)
}


plot(tr0  <- pstab.tailratio(1, 0.5))
plot(tr1  <- pstab.tailratio(1.1, 0.25))
plot(tr2 <- pstab.tailratio(0.99, +0.992))

showProc.time()

plot(tr   <- pstab.tailratio(1.2, 0.5))

plot(tr3 <- pstab.tailratio(0.7, +0.9))

plot(tr4 <- pstab.tailratio(1.7, +0.6))# not really useful: pstable(.) = 1 too early

showProc.time()

##---------------- Now the density

##' @title Explore eps(x) where   dstable(x)/dPareto(x) = 1 + eps(x)
##' @param alpha
##' @param beta
##' @param n
##' @param prob
##' @param xmin
##' @param xmax
##' @return an object of \code{\link{class} "dstableTailratio"}, ...
##' @author Martin Maechler, 21 Mar 2011
dstab.tailratio <- function(alpha, beta, n = nc, prob = 1/4096,
                            xmin = qstable(0.95, alpha,beta, tol = 0.01),
                            xmax = qstable(1 - prob, alpha,beta))
{
    cl <- match.call()
    stopifnot(0 < prob, prob < .5,
              0 < xmin, xmin < xmax)
    x <- exp(seq(log(xmin), log(xmax), length.out = n))
    f <-               dstable(x, alpha,beta)
    ok <- f > 0
    f <- f[ok]
    x <- x[ok]
    fp <- stabledist:::dPareto(x, alpha,beta)
    eps <- (f - fp)/fp
    structure(list(x=x, eps=eps, call = cl, alpha=alpha, beta=beta),
              class = "dstableTailratio")
}

##' @title plot() method for  dstableTailratio() results
##' @param x object of \code{\link{class} "dstableTailratio"}.
##' @param type plot type; default simple lines
##' @param col
##' @param lin.col
##' @param ... optional further arguments passed to \code{\link{plot.formula}()}.
##' @return
##' @author Martin Maechler
plot.dstableTailratio <- function(x, type="l", col="blue3",
                                  lin.col = adjustcolor("red2",.6), ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
    stopifnot(is.list(x), is.numeric(x$eps))
    tit <- substitute("tail ratio approx. for"~~ dstable(alpha==A, beta==B),
                      list(A=x[["alpha"]], B=x[["beta"]]))
    dat <- as.data.frame(x[c("x","eps")])
    ## NOTA BENE: Empirically,  I see that  eps > 0 <==> alpha > 1
    ##                                      eps < 0 <==> alpha < 1
    dat <- dat[dat[,"eps"] > 0, ] ## drop invalid  eps[]
    plot(eps ~ x, log = "xy", data = dat, col=col,
         ylab = expression(epsilon ~~~ "('eps')"),
         main = tit, type=type, ...)
    mtext( expression(epsilon(x) == (f(x,.) - f[P](x,.)) / f[P](x,.)) )
    fm <- lm(log(eps) ~ log(x), weights = x^2, data = dat)
    lines(dat[["x"]], exp(predict(fm)), col=lin.col)
    Form <- function(x) formatC(x, digits=4, wid=1)
    leg.line <-
        substitute(log(epsilon) == A + B * log(x),
                   list(A = Form(coef(fm)[["(Intercept)"]]),
                        B = Form(coef(fm)[["log(x)"]])))
    legend("topright", legend=leg.line, bty = "n", lty=1, col=lin.col)
}

plot(fr   <- dstab.tailratio(1.01, 0.8))
plot(fr   <- dstab.tailratio(1.05, 0.4))
plot(fr   <- dstab.tailratio(1.1,  0.4))
plot(fr   <- dstab.tailratio(1.2,  0.5))
plot(fr   <- dstab.tailratio(1.3,  0.6))

showProc.time()

plot(fr   <- dstab.tailratio(1.4, 0.7))
plot(fr   <- dstab.tailratio(1.5, 0.8))

plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1000))
plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1e4));abline(v=1000, lty=2)
plot(fr   <- dstab.tailratio(1.5, 0.8, xmax= 1e5));abline(v=1e4, lty=2)

showProc.time()

plot(fr   <- dstab.tailratio(1.6, 0.9))
plot(fr   <- dstab.tailratio(1.7, 0.1))
plot(fr   <- dstab.tailratio(1.8, 0.2))

showProc.time()

##------ Some explicit tail problems visualized:

I <- integrate(dstable, 0,Inf, alpha=0.998, beta=0, subdivisions=1000)
str(I)
stopifnot(abs(I$value - 0.5) < 1e-4)

curve(dstable(x, alpha=.998, beta=0,    log=TRUE), 10, 1e17, log="x")
curve(dstable(x, alpha=.999, beta=0.1,  log=TRUE), 10, 1e17, log="x")
curve(dstable(x, alpha=.999, beta=0.9,  log=TRUE), 10, 1e17, log="x")
curve(dstable(x, alpha=.999, beta=0.99, log=TRUE), 10, 1e17, log="x")
curve(dstable(x, alpha=.999, beta=0.99, log=TRUE), 10, 1e170, log="x")
showProc.time()

## less problems when alpha > ~= 1 (but it's  __S..L..O..W__ !)
curve(dstable(x, alpha=1.001,beta=0.99, log=TRUE), 10,  1e7, log="x")
curve(dstable(x, alpha=1.001,beta=0.99, log=TRUE), 10, 1e17, log="x")
## -> problem --> zoom in:
curve(dstable(x, alpha=1.001,beta=0.99, log=TRUE), 1e12, 160e12)
curve(dPareto(x, alpha=1.001,beta=0.99, log=TRUE), add=TRUE, lty=3, col=4)

curve(dstable(x, alpha=1.001,beta=0.99, log=TRUE), 10, 1e40, log="x")
showProc.time()

## NB: alpha == 1   also has problems in tail  --- only as long as "old R"s wrong uniroot is used:
curve(dstable(x, alpha=1.   ,beta=0.99, log=TRUE), 1, 20)
curve(dstable(x, alpha=1.   ,beta=0.99, log=TRUE), 1,100)


showProc.time()
