## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

###  Generating Exponentially Tilted Stable Random Vars (r ET stable)
###  ================================================================
###  --> Experiments with retstable*() versions

require(copula)
## source(system.file("Rsource", "utils.R", package="copula"))
##--> tryCatch.W.E(), canGet()
## NB:  Since  Feb.2012 -- the package ../inst/rData/ has the saveFile s
##        Size  Name
## ----------------------------------
##       16205  retstable_CPU2.rda
##      304463  retstable_Nstat.rda
##       16415  retstable_st2.rda

### --- ======      --------------------------#--
### --- Part I ---  Experiments with retstableR()
### --- ======      --------------------------#--

##' Zolotarev's function A(.) to the power 1-alpha,
##' i.e. sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x)
##' and other things
##'
##' @title Zolotarev's function A
##' @param ialpha parameter 1-alpha
##' @param x sequence of points in (0,pi)
##' @param col color
##' @return a plot
##' @author Martin Maechler
p.A <- function(ialpha, x = seq(0, 3.1, length=100), col = "tomato") {
    stopifnot(0 <= ialpha, ialpha <= 1, 0 <= x, x <= pi)
    if(is.unsorted(x)) x <- sort(x)
    tit <- substitute(list("Zolotarev's  "* {A(x, alpha) - 1} ,
                           "  "* 1-alpha == IA), list(IA = ialpha))
    alpha <- 1 - ialpha
    A1 <-  A..Z(x, alpha, ialpha) - 1
    A1d <- A..Z(x, alpha, 1-alpha) - 1
    plot(x, A1, ylim = range(A1, A1d), ylab = expression(A(x,alpha) - 1),
         col = col, main = tit, type ="o")
    abline(h=0, lty=2)
    ## add the "dumb" version that just works with alpha
    gray <- rgb(.4, .5, .6, alpha = .3) # opaque color
    lines(x, A1d, col = gray, lwd = 3)
}

if(doPDF <- !dev.interactive(orNone=TRUE))
    pdf("retstable-ex.pdf")

## A(.) --> 1
par(mfrow = c(2,2), mar = .1+c(4,4,4,1), mgp = c(1.5, .6, 0))
p.A(1e-7)
p.A(1e-9)
p.A(1e-12)
p.A(1e-14) # still fine, visually
p.A( 1e-15) # wiggliness!
p.A(.8e-15) # more wiggliness; dumb version very slightly differs
p.A(.5e-15) # more -- dumb version differs clearly
p.A(.2e-15) # ditto
p.A(.1e-15) # even more -- but this is already < .Machine$double.eps:
.Machine$double.eps / .1e-15
p.A(.08e-15) # = 8e-17
p.A(.06e-15) # = 6e-17 still just wiggliness
p.A(.05e-15) # complete breakdown

## Conclusion: for the small range of  1-alpha in [0.6, 8] * 10^-17
## ----------  using 'accurate 1-alpha' instead of just compute "1 - alpha"
##             *is* somewhat more accurate

## investigate the  m(V0) = "optimal m" : ----------------------

V0 <- seq(0, 40, by = 1/256)
nm <- length(m <- sapply(V0, copula:::m.opt.retst))
iJmp <- which(diff(m) != 0)
stopifnot( all.equal(m[c(iJmp,nm)], 1:40) )
iKept <- c(1, rbind(iJmp, iJmp+1), nm)

V0. <- V0[iKept]
m. <-  m[iKept]

V0[iJmp] - m[iJmp] ## 0.3828125 0.4296875 0.4492188 0.4609375 => conv. to 0.5
V0[iJmp+1] - m[iJmp+1] ## .... converging to -0.5


## ------------ How many times is rstable1(1, ...) called ? ------------------
## --> we *really* should call  rstable1(N, ....)
## and then use these values where N = N(V0):

saveFile <- system.file("rData", "retstable_Nstat.rda",
                        package = "copula")

## This needs roughly two and a half hours -- elapsed time !!
if(nzchar(saveFile) && file.exists(saveFile)) {
    load (saveFile)
} else {
    .e <- new.env()
    ## our tracer function is just a "counter":
    rstTR <- function() .e$N <- if(is.null(.e$N)) 0L else .e$N + 1L
    trace(rstable1, rstTR, print=FALSE)

    ## For now, we use the R (reliable / correct) version:

    nV <- length(V0s <- c(1:4, 5*(1:10), 10*(6:20)))
    nAlph <- length(alphas <- (1:9)/10)
    nSim <- 1024

    Nstat <- array(dim = c(nSim, nAlph, nV))
    for(iV0 in 1:nV) {
        V0 <- V0s[iV0]
        cat(sprintf("V0 = %5d :", V0))
        for(ia in 1:nAlph) {
            alpha <- alphas[ia]
            Nstat[, ia, iV0] <- replicate(nSim,
                                      { .e$N <- 0L # start counting ...
                                        RR <- retstableR(alphas[ia], V0 = V0)
                                        .e$N }
                                          )
            cat(".")
        }
        cat("\n")
    }

    save(V0s,alphas, nSim, Nstat, file = saveFile)
}

## Visualize  'Nstat' -- can do this after loading the file

boxplot(Nstat[,1,], range=2)

## Log log scale
op <- if(require("sfsmisc"))
    sfsmisc::mult.fig(length(alphas))$old.par else { ## less nicely looking
        m <- length(alphas); par(mfrow= c(3, ceiling(m/3))) }
for(i in seq_along(alphas)) {
    N.st <- Nstat[,i,]
    plot(NA,NA, xlim = range(V0s), ylim=range(N.st),
         main = substitute(N == no({"calls to_ " * rstable1()}) * " for  "
         * { alpha == A },
         list(A = alphas[i])),
         xlab=expression(V[0]), ylab = "Nstat", type = "n", log="xy")
    boxplot(N.st, at = V0s, range=2, add = TRUE)
    lines(apply(N.st, 2, mean, trim=.05) ~ V0s, col = 2, lwd=2)
}
par(op)

str(Nst.mns <- apply(Nstat, 3:2, mean, trim=.05))
str(Nst.q95 <- apply(Nstat, 3:2, quantile, prob=.95))
str(Nst.max <- apply(Nstat, 3:2, max))
str(Nst.min <- apply(Nstat, 3:2, min))

## need  k * V0 :  what k
## in the mean:
rowMeans(Nst.mns / V0s) ## ?== Euler's  e == exp(1)
## maximally
rowMeans(Nst.max / V0s)

## is it  'e'  ??
plot(V0s, rowMeans(Nst.mns / V0s), log="x"); rug(V0s); abline(h=exp(1), col=2)
L <- V0s >= 20
plot(V0s[L], rowMeans(Nst.mns[L,] / V0s[L]), log="x"); rug(V0s[L])
abline(h=exp(1), col=2)
## well, as a limit probably

mtit <- substitute(N == `# `*group("{", rstable1() * " calls", "}") *
                   "   for  " * alpha %in% group("{", list(A,B,...,Z), "}"),
                   list(A=alphas[1], B=alphas[2], Z=tail(alphas,1)))

## normal scale
matplot (V0s, Nst.max, ylim=range(Nstat), col = "gray20",
         ylab="", main=mtit, type="l", lty=2); rug(V0s)
matlines(V0s, Nst.min, col = "gray20", lty=2)
abline(a=0,b=1, col = "tomato"); text(20,20, expression(y == x), col="tomato",
                adj=-.2)
matlines(V0s, Nst.mns, type = "b", lty=1)
matlines(V0s, Nst.q95, col = "gray80", lty=3)

## log log scale
matplot (V0s, Nst.max, ylim=range(Nstat), col = "gray20",
         ylab="", main=mtit, type="l", lty=2, log="xy"); rug(V0s)
matlines(V0s, Nst.min, col = "gray20", lty=2)
abline(a=0,b=1, col = "tomato"); text(20,20, expression(y == x), col="tomato",
                adj=-.2)
matlines(V0s, Nst.mns, type = "b", lty=1)
matlines(V0s, Nst.q95, col = "gray80", lty=3)


## or just the "residuals"
matplot(V0s, Nst.mns / V0s, type = "b", lty=1, log = "x")
rug(V0s)

if(doPDF) { dev.off(); pdf("retstable-ex-2.pdf") } else par(mfrow=c(1,1))

### --- =======      --------------------------#-----------
### --- Part II ---  Experiments with retstableC() methods
### --- =======      --------------------------#-----------

##' mean function
meanretstable <- function(alpha,h,V0) alpha*V0*h^(alpha-1)

nalpha <- length(alpha <- c(0.05, (1:9)/10, 0.95))
nV0    <- length(  V0 <- c(0.2,0.5,1,2:5,10))
nh     <- length(   h <- c(0.5, 0.75, 1, 1.5, 2, 5, 10))
meth <- c("MH", "LD")
nsim <- 100000 # a lot
hist.breaks <- 200 # == function(nsim, <parameters>)

saveFile2 <- system.file("rData", "retstable_st2.rda", package = "copula")
saveFile3 <- system.file("rData", "retstable_CPU2.rda",package = "copula")

if(file.exists(saveFile2) && file.exists(saveFile3)) {
    ## we have precomputed it ...
    load (saveFile2)
    load (saveFile3)
} else {

    set.seed(47)
    dnS <- list(alpha = paste("alpha", alpha, sep="="),
		V0 = paste("V0", formatC(V0), sep="="),
		h = paste("h", formatC(h), sep="="),
                meth = meth)
    dim.S <- c(nalpha, nV0, nh, length(meth))
    St.c <- vector("list", prod(dim.S))
    dim(St.c) <- dim.S
    dimnames(St.c) <- dnS
    St.. <- array(dim = c(nsim, length(meth)),
                  dimnames = list(NULL, meth))
    CPU.c   <- array(dim = dim(St.c), dimnames = dnS)
    tt.pval <- array(dim = dim(St.c), dimnames = dnS)
    KS.pval <- array(dim = dim(St.c)[1:3], dimnames = dnS[1:3])

    count <- 1
    for(ia in 1:nalpha) {
	alph <- alpha[ia]
	for(iV0 in 1:nV0) {
	    V0. <- V0[rep.int(iV0,nsim)]
	    cat("alpha=",alph,", V0=",V0[iV0],", h= ")
	    for(ih in 1:nh) {
		cat(h[ih],"")
		mu. <- meanretstable(alph, h = h[ih], V0 = V0.[1])
		for(met in meth) {
		    ## generate data, measure run time
		    CPU.c[ia, iV0, ih, met] <-
			system.time(
				    St..[,met] <-
				    retstable(alph, V0= V0., h= h[ih], method = met)
				    )[1]
		    ## (asymptotical) two-sided t-test
		    tt.pval[ia, iV0, ih, met] <- t.test(St..[,met], mu = mu.)$p.value
		    St.c[ia, iV0, ih, met] <-
			list(hist(log(St..[,met]),
				  breaks = hist.breaks, plot = FALSE))
		}
		KS.pval[ia, iV0, ih] <- ks.test(St..[,"MH"],
						St..[,"LD"])$p.value
	    }
	    cat("\nProc.time(): ",format(proc.time()[1]),"; ",
		format(100*count/(nalpha*nV0), digits=3), "% done\n\n",sep="")
	    count <- count + 1
	}
    }
    save(St.c,  alpha, V0, h, meth, nsim, file = saveFile2)
    save(CPU.c, KS.pval, tt.pval,
         alpha, V0, h, meth, nsim, file = saveFile3)
}


## ---> Compare  CPU.c for the two methods,
## --- statistically compare the two methods for St.c

if(getOption("width") < 100) options(width=100)


##' check the random variates via histogram plots
##'
##' @title Graphical tests via histogram plots
##' @param alphalab alpha label
##' @param V0lab V0 label
##' @param hlab h label
##' @param main title
##' @param nBreaks number of breaks
##' @param log log-scale
##' @return histogram
##' @author Martin Maechler
histSt <- function(alphalab, V0lab, hlab, main, nBreaks = 100, log = TRUE) {
    stopifnot(## the "globals":
              is.list(St.c), is.array(St.c), length(dnS <- dimnames(St.c)) == 4,
              is.character(meth), is.character(alphalab), is.character(V0lab),
              is.character(hlab), any(alphalab == dnS[[1]]),
              any(V0lab == dnS[[2]]), any(hlab == dnS[[3]]))
    if(missing(main) || is.null(main))
	main <- paste(format(nsim, sci=FALSE),
		      "	 expo.tilted stable vars w/ ",
		      alphalab,", ",V0lab,", ",hlab,
                      if(log) "  LOG scale", sep = "")
    nmeth <- length(meth)
    op <- if(require("sfsmisc"))
        mult.fig(mfrow= c(nmeth,1), main = main)$old.par else { ## less nicely looking
            cex.m <- par("cex.main")
            pp <- par(mfrow= c(nmeth,1), oma = c(0,0, 1+ 1.5*cex.m, 0),
                      mgp = c(1.5, 0.6, 0))
            plot.new()
            mtext(main, side = 3, outer = TRUE, line = cex.m - .5,
                  cex = cex.m, font = par("font.main"), col = par("col.main"))
            par(new = TRUE)
            pp
        }
    on.exit(par(op))

    S. <- St.c[alphalab, V0lab, hlab, ]
    ##f if(log) {
    ##f    S. <- log(S.)
    xl <- expression( log( tilde(S) ) )
    ##f } else
    ##f xl <- expression( tilde(S) )
    ##f breaks <- pretty(range(S.), n = nBreaks)
    for(me in meth) {
        ## hist(S.[me,], breaks = breaks, main = me, xlab = xl, freq=FALSE)
        plot(S.[[me]], main = me, xlab = xl, freq=FALSE)
        ##f lines(density(S.[me]), col=2, lwd=2)
    }
}

histSt("alpha=0.3", "V0=5", "h=1")

histSt("alpha=0.5", "V0=1", "h=1")
histSt("alpha=0.1", "V0=0.5", "h=1")

## For  h != 1
## Hmm:  MH and LD  look different here : this is suspicious:
histSt("alpha=0.1", "V0=10", "h=2")
histSt("alpha=0.1", "V0=5", "h=5")

histSt("alpha=0.8", "V0=10", "h=2")

histSt("alpha=0.8", "V0=5",  "h=10")
histSt("alpha=0.9", "V0=5",  "h=10")

## no problem for  h == 1 :
histSt("alpha=0.9", "V0=5",  "h=1")

## Instead of "just" the plots, now also look at the P-values:
## For both tests, the P-values look like ~ U[0,1]  as they should:
hist(c(KS.pval))
hist(c(tt.pval))
## More formally, they do look like uniform ([0,1]):
ks.test(tt.pval, "punif")
ks.test(KS.pval, "punif")

## Times for nsim replicates for given alphas, V0s, hs, and methods
## ratio of the two methods : decision at  r == 1   <==>   log(r) == 0
CPUr <- CPU.c[,,,"MH"] / CPU.c[,,,"LD"]

plot(density(log10(CPUr))); rug(log10(CPUr))
if(FALSE) { # hmm, not sensical yet
    plot(density(log10(CPUr)), log="x", xaxt = "n")
    rug(log10(CPUr))
}
## x-range  in log-scale and back-transformed
x.r <- 10^(xLr <- par("usr")[1:2])
(x. <- floor(xLr[1]):ceiling(xLr[2]))
x. <- sort(10^x. %*% t(c(1,2,5))); x. <- x.[(10^xLr[1] <= x.) &
                                             (x. <= 10^xLr[2])]
x.
axis(1, at = log10(x.), label = x.)
abline(v=0, col = "skyblue", lty=2)

signif(CPUr[,, "h=1"], 3)
## --> the boundary is simply between V0=2 and V0 =5 (!)
signif(CPUr[,, "h=5"], 2) # slightly different picture


## Workaround for now:
(V0 <- as.numeric(sub("^V0=","", dimnames(CPUr)[[2]])))
(h  <- as.numeric(sub("^h=","",  dimnames(CPUr)[[3]])))

## or use the lattice equivalent
filled.contour(alpha, V0, log10(CPUr[,, "h=1"]))
filled.contour(alpha, V0, log10(CPUr[,, "h=5"]))


## Theory:  expect boundary   h^alpha * V0  <= c
##                            ------------------
str(dCPU <- cbind(expand.grid(alpha=alpha, V0=V0, h=h), CPUr = c(CPUr)))
rm(CPUr)# (so we know that dCPU is used below)
## original scale
plot  (CPUr ~ I(h^alpha * V0), data = dCPU); abline(h=1, col="tomato")

##' sophisticated plot of CPU times
##'
##' @title CPU times plot
##' @param log scale
##' @param do.h.eq.1 h == 1
##' @param main title
##' @param data a dataframe "like 'dCPU'" (typically subset of it)
##' @param col vector of colors,
##' @param pch vector of point characters
##' @param ... further arguments to plot()
##' @return plot of CPU times
##' @author Martin Maechler
p.cpu <- function(log = "", do.h.eq.1 = TRUE,
		  main = expression(t[CPU](MH) / t[CPU](LD) *
		      " as " * f(h, alpha, V[0])),
		  data = dCPU, colT = NULL, pchT = NULL,
		  pal = { if(requireNamespace("RColorBrewer"))
			      RColorBrewer::brewer.pal(nalpha, "RdYlGn") # "Spectral"
			  else c("#A50026", "#D73027", "#F46D43", "#FDAE61",
				 "#FEE08B", "#FFFFBF", "#D9EF8B", "#A6D96A",
				 "#66BD63", "#1A9850", "#006837")
		      },
		  ... )
{
    doCol <- !is.null(colT) ; if(!doCol) col <- par("col")
    doPch <- !is.null(pchT) ; if(!doPch) pch <- par("pch")
    if(doCol || doPch) {
	stopifnot(colT %in% c("V0","alpha","h"))
	stopifnot(pchT %in% c("V0","alpha","h"))
	a.c <- as.numeric(a.f <- as.factor(data[,"alpha"]))
	V.c <- as.numeric(V.f <- as.factor(data[,"V0"]))
	h.c <- as.numeric(h.f <- as.factor(data[,"h"]))
        a.f <- paste("alpha", levels(a.f), sep="=")
        V.f <- paste("V0", levels(V.f), sep="=")
        h.f <- paste("h", levels(h.f), sep="=")
        nalpha <-
	getC <- function(typ)
	    switch(typ, "V0" = V.c, "alpha" = a.c, "h" = h.c)
	getL <- function(typ)
	    switch(typ, "V0" = V.f, "alpha" = a.f, "h" = h.f)
	if(doPch) {
	    pch <- getC(pchT)
	    pLabs <- getL(pchT)
	}
	if(doCol) {
	    col <- 1 + getC(colT)
	    cLabs <- getL(colT)
	    ## also set palette()!
	    oPal <- palette(pal)
	    ## on exit, revert to previous color palette:
	    on.exit(palette(oPal))
	}
    }
    plot(CPUr ~ I(h^alpha * V0), main=main, xlab =
	 expression(h ^ alpha %*% V[0]), data=data, log=log,
         col=col, pch=pch, ...)
    abline(h = 1, lty=2, col="tomato") # ratio == 1  <==> "MH" as fast as "LD"
    if(do.h.eq.1) {
	points(CPUr ~ I(h^alpha * V0), data=data, pch = 2, col = "red",
	       subset= h == 1)
	legend("bottom", "h = 1", pch = 2, col = "red", inset=.01)
    }
    if(doPch)
	legend("bottomright", legend = pLabs,
	       col = par("col"), pch = unique(pch), inset=.01)
    if(doCol)
	legend("topleft", legend = cLabs,
	       pch = par("pch"), col = unique(col), inset=.01)
}

p.cpu()
p.cpu(log="xy")
## cool! -- the formula seems correct!

p.cpu(log="xy", colT = "alpha", pchT = "V0")
## Zoom in:
p.cpu(log="xy", data = subset(dCPU, subset = 3 <= V0 & V0 <= 5),
      do.h.eq.1=FALSE, colT = "alpha", pchT = "V0")# still a bit wide x-range
mtext("3 <= V0 <= 5")

## zoom in more, no longer log-scaled:
p.cpu(data = subset(dCPU, subset = 3 <= V0 & V0 <= 5),
      xlim = c(2,5), ylim=c(0.5,2), # restricting range
      do.h.eq.1=FALSE, colT = "alpha", pchT = "V0")
mtext("3 <= V0 <= 5")




summary(mod1 <- lm(log(CPUr) ~ I(alpha*log(h)) + log(V0), data = dCPU))
plot(mod1)
summary(mod2 <- lm(log(CPUr) ~ alpha*log(h) + log(V0), data = dCPU))
plot(mod2)
anova(mod1,mod2) # -> second model is not better
plot(residuals(mod1) ~ alpha, data=dCPU)
plot(residuals(mod1) ~ h, data=dCPU)
plot(residuals(mod1) ~ I(alpha*log(h)), data=dCPU)

## Only look "around"  CPUr == 1  as we want to choose for that
dim(s.CPU <- with(dCPU, dCPU[1 <= V0 & V0 <= 5,])) # - only 165 obs
summary(mod.1 <- lm(log(CPUr) ~ I(alpha*log(h)) + I((alpha*log(h))^2) + log(V0),
                    data = s.CPU))
## coefficient of log(V0) is ~= 1 ---> set it == 1  :
summary(mod.2 <- lm(log(CPUr / V0) ~ I(alpha*log(h)) + I((alpha*log(h))^2),
                    data = s.CPU))
## or even better
summary(mod.3 <- lm(log(CPUr / V0) ~ I(alpha*log(h)^2), data = s.CPU))

plot(mod.1)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
