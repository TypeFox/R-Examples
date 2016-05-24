#### Psi(), Rho(), weight() etc functions for  M-Estimation and extensions

## Use an  S4 class for such function classes
## Follow a similar idea as  nlsModel() {in "stats"} which returns
##  a list of functions sharing a common {non-small!} environment

## NOTA BENE:  Experiments etc are currently in ../misc/experi-psi-rho-funs.R
## ---------   (FIXME: move those to ../tests/psi-rho-etc.R and the vignette
## ../vignettes//psi_functions.Rnw  (and see ../inst/xtraR/plot-psiFun.R)

## ---> look for 'FIXME' below !!!
##               -------

### A.  (Symmetric) Location / Regression

## A single   function(x, tuningPars)
## a. 1st argument 'x', numeric; must work vectorized on x
## b. further arguments: tuning parameters *with a default*
setClass("functionX", contains = "function",
         validity = function(object) {
             ## "function" is already because of 'contains'
             if(names(ff <- formals(object))[1] != "x")
                 return("first argument must be 'x'")
             f0 <- object(0)
             fI <- object(Inf)
             if(!identical(c(f0,fI), object(c(0,Inf))))
                 return("F(x, *) does not vectorize in 'x'")
             ## Otherwise : valid
             TRUE
         })

## A functional --- i.e. function of "tuning pars only",  such as
##  Ep(hc) = Int_{-Inf}^{+Inf} psi(x; hc)^2 dnorm(x) dx

##' This one is *not* checked for vectorization: needed when length(k) > 1
setClass("functionXal", contains = "function")

##' Here F(k) must vectorize in k
setClass("functionXal1", contains = "functionXal",
         validity = function(object) {
             f0 <- object(0)
             fI <- object(Inf)
             if(!identical(c(f0,fI), object(c(0,Inf))))
		 return("F(k) = I_k[f(.)] does not vectorize in 'k'")
             ## Otherwise : valid
             TRUE
         })

setClass("psi_func",
         representation(rho = "functionX",
                        psi = "functionX", ## psi(x) == d/dx rho(x) = x * wgt(x)
                        wgt = "functionX", ## wgt(x) == psi(x) / x
                        Dpsi = "functionX",## psi'(x) == d/dx psi(x) = rho''(x)
                        Dwgt = "functionX", ## wgt'(x) == d/dx wgt(x)
                        ## tuning parameters, i.e., formals(rho)[-1]
                        tDefs = "numeric",## *named* values of tuning parameters
                        ## FIXME !! {see 4 lines below}
                        Erho =  "functionXal", # = E_X[rho(X)];   X~N(0,1);
                        Epsi2 = "functionXal", # = E_X[psi(X)^2]; X~N(0,1); 'A'
                        EDpsi = "functionXal", # = E_X[psi'(X)];  X~N(0,1); 'B'
                        ##
                        name = "character",
                        xtras = "list" ## for flexible extensions..
                        ))
## FIXME: need other E[] than just wrt N(0,1)
## -----  e.g. for robglm(), need  E[] wrt Gamma(.)

### Constructors / "Examples" [the examples are the objects, we'll really use!]

psiFunc <- function(rho,psi,wgt, Dpsi,Dwgt, Erho=NULL, Epsi2=NULL, EDpsi=NULL, name, ...)
{
    lent <- length(dotsargs <- list(...))
    ## '...'  must contain all tuning parameters and their defaults:
    ## NOTA BENE: Now want at least one tuning parameter.. "worst case": a dummy
    stopifnot(lent >= 1, length(nt <- names(dotsargs)) == lent,
              all(nchar(nt)) >= 1)

    ## Definition of Dwgt is optional
    if (missing(Dwgt)) Dwgt <- .defDwgt(psi, Dpsi)

    ## rho, psi,... checking: must have argument names
    argn <- c("x", nt)
    for(fnam in list("rho", "psi", "wgt", "Dpsi", "Dwgt",
                     "Erho", "Epsi2", "EDpsi")) {
        f <- get(fnam, inherits = FALSE)
        ef <- environment(f)
        nf <- names(formals(f))   # "x" and "k" for Huber's
        if (fnam %in% c("Erho", "Epsi2", "EDpsi")) {
            if(!identical(nf, argn[-1]))
                stop("arguments of function '",fnam,"' are (",
                     paste(nf,  collapse=","),") but should be (",
                     paste(argn[-1],collapse=","),").")

            formals(f) <- dotsargs
        } else {
            if(!identical(nf, argn))
                stop("arguments of function '",fnam,"' are (",
                     paste(nf,  collapse=","),") but should be (",
                     paste(argn,collapse=","),").")

            formals(f)[-1] <- dotsargs
        }
        environment(f) <- ef
        assign(fnam, f, inherits = FALSE)
    }
    fnctl.typ <- if(lent == 1 && length(dotsargs[[1]]) == 1)
        "functionXal1" else "functionXal"
    new("psi_func",
	rho = new("functionX", rho),
	psi = new("functionX", psi),
	wgt = new("functionX", wgt),
	Dpsi= new("functionX", Dpsi),
        Dwgt= new("functionX", Dwgt),
	## tNams = if(lent) nt else character(0),
	tDefs = unlist(dotsargs),
	Erho = new(fnctl.typ, Erho),
	Epsi2= new(fnctl.typ, Epsi2),
	EDpsi= new(fnctl.typ, EDpsi),
        name = if (missing(name)) character(0) else name,
	xtras= list(tuningP = dotsargs))
}

## Generate default Dwgt function

## Unfortunately, MM can't see how to make this works nicely;
## ._.. = args should really be something like  'x, k' {no parens}:
.defDwgt <- function(psi, Dpsi) {
    args <- formals(Dw <- psi)# -> same formals
    body(Dw) <- substitute({
        y <- .X.
        .X. <- .X.[not0 <- .X. != 0]
        y[not0] <- ( Dpsi(._..) - psi(._..)/.X. ) / .X.
        y
    }, list(.X. = as.name(names(args[1])), ._.. = args))
    environment(Dw) <- environment()
    Dw
}
## so we use this "less nice" variant:
.defDwgt <- function(psi, Dpsi) {
    nf <- names(formals(psi))
    eval(parse(text =
	       gsub("_,_", paste(nf, collapse=","),
		    gsub("x", nf[1], "function(_,_) {
        y <- x
        x <- x[not0 <- x != 0]
        y[not0] <- ( Dpsi(_,_) - psi(_,_)/x ) / x
        y
    }"))))
}

chgDefaults <- function(object, ...)
    standardGeneric("chgDefaults")

setMethod("chgDefaults", signature("psi_func"),
          function(object, ...)
      {
          lent <- length(dotsargs <- list(...))
          ## '...'  must contain all tuning parameters and their defaults:
          stopifnot(lent >= 1, length(nt <- names(dotsargs)) == lent,
                    all(nchar(nt)) >= 1)
          ## rho "..." must conform to rho, etc:
          nf <- names(ff <- formals(object@rho))
          if(!identical(nf[-1], nt))
              stop("invalid tuning parameter names: ",
                   paste(nt,    collapse=",")," instead of ",
                   paste(nf[-1],collapse=","),".")

          for(fnam in list("rho", "psi", "wgt", "Dpsi", "Dwgt",
                           "Erho", "Epsi2", "EDpsi")) {
              f <- slot(object, fnam)
              ef <- environment(f)
              if (is(f, "functionXal"))
                  formals(f) <- dotsargs else formals(f)[-1] <- dotsargs
              environment(f) <- ef
              ## lowlevel {faster than}: slot(..) <- new("functionX", f)
              slot(object, fnam)@.Data <- f
          }
          object@tDefs <- unlist(dotsargs)
          if(identical(nt, names(object@xtras$tuningP)))# TODO: should update even if there are others
              object@xtras$tuningP <- setNames(eval(dotsargs), nm=nt)
          object
      })

.sprintPsiFunc <- function(x, short=FALSE, round=3) {
    v <- x@tDefs
    n <- names(v)
    ## do not print a single dummy parameter "."
    if (length(n) == 1 && n == ".") v <- numeric(0)
    if (!length(name <- x@name)) name <- "<unnamed>"
    if (!short) name <- sprintf("%s psi function", name)
    if (length(v) >= 1) {
        if (short)
            paste(name, paste(n, round(v, round), sep = "=", collapse = "\n"),
                  sep = "\n")
        else
	    paste0(name, " (", pasteK(n, round(v, round), sep = " = "), ")")
    } else name
}

setMethod("show", signature("psi_func"),
          function(object) cat(.sprintPsiFunc(object), "\n"))

## moved here from inst/xtraR/plot-psiFun.R; called  plot.psiFun  originally
matplotPsi <- function(x, m.psi, psi, par, main = "full",
			col = c("black", "red3", "blue3", "dark green"),
			leg.loc = "right", lty = 1, ...) {
    ## Original Author: Martin Maechler, Date: 13 Aug 2010, 10:17
    ## Modified by Manuel Koller, Date: 7 Jan 2013
    fExprs <- quote(list(rho(x), psi(x), {psi*minute}(x),
			 w(x) == psi(x)/x, {w*minute}(x)))
    ## build legend
    map <- if (is.null(colnames(m.psi))) {
	1:(ncol(m.psi)+1)
    } else {
	c(1, c(rho=2, psi=3, Dpsi=4, wgt=5, Dwgt=6)[colnames(m.psi)])
    }
    fExprs <- fExprs[map]
    ## ... title
    if(is.character(main)) {
	shortMain <- (main == "short")
	elist <- list(FF = if(shortMain) fExprs[[2]] else fExprs,
		      PSI = psi, PPP = paste(formatC(par), collapse=","))
	tit <- if(shortMain)
	    substitute(FF ~ "etc, with"  ~ psi*"-type" == PSI(PPP), elist)
	else
	    substitute(FF ~~ ~~ " with "~~ psi*"-type" == PSI(PPP), elist)
    } else tit <- NULL
    ## plot
    matplot(x, m.psi, col=col, lty=lty, type="l", main = tit,
	    ylab = quote(f(x)), xlab = quote(x), ...)
    abline(h=0,v=0, lty=3, col="gray30")
    fE <- fExprs; fE[[1]] <- as.name("expression")
    legend(leg.loc, inset=.02, eval(fE), col=col, lty=lty, bty="n")
    invisible(cbind(x=x, m.psi))
}

setMethod("plot", signature(x = "psi_func"),
	  function(x, y, which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
		   main = "full",
		   col = c("black", "red3", "blue3", "dark green", "light green"),
		   leg.loc = "right", ...) {
	      ## x: psi_func object
	      ## y: points to plot at (x-Axis in plot)
	      which <- match.arg(which, several.ok = TRUE)
	      if(missing(y)) y <- seq(-5, 10, length=1501)
              ## For backcompatibility:
              if(!is.null(sm <- list(...)$shortMain)) {
                  if(!missing(main))
                      stop("You cannot specify both 'main' and the deprecated 'shortMain'")
                  warning("'shortMain' is deprecated and will get defunct.\n",
                          "Use 'main = \"short\"' instead of 'shortMain = TRUE'")
                  if(sm) main <- "short"
              }
	      tmp <- lapply(which, function(name) slot(x, name)(y))
	      m.psi <- do.call(cbind, tmp)
	      colnames(m.psi) <- which
	      matplotPsi(y, m.psi, x@name, unlist(formals(x@rho)[-1]),
                         main=main, col=col, leg.loc=leg.loc, ...)
	  })

##-------- TODO: Rather right short  vignette with these formulae

##' \Phi_j(t) := \int_{-\infty}^t  u^j \phi(u) \;du
##' ---------    where \phi(.) (= \code{dnorm()})
##'              is the density of the standard normal distribution  N(0,1).
##' @title "Truncated" Moments of the Gaussian: Int u^j phi(u) du
##' @param t numeric vector
##' @param j an integer (valued scalar), >= 0
##' @return Phi_j(t), i.e. a numeric vector of the same length as t.
##' @author Martin Maechler
PhiI <- function(t, j = 0) {
    stopifnot(j == as.integer(j), length(j) == 1, is.numeric(t))
    if(j >= 4) ## recursion formula
        -t^(j-1)*dnorm(t) + (j-1)* PhiI(t, j-2)
    else
        switch(j+1,
               ## 0:
               pnorm(t),
               ## 1:
               -dnorm(t),
               ## 2:
               pnorm(t) - t*dnorm(t),
               ## 3:
               -(2 + t^2)*dnorm(t))
}

if(FALSE) { ## Checking  PhiI() visually:

    tt <- seq(-4,10, length=64)
    j.max <- 5
    oo <- sfsmisc::mult.fig(j.max+1, main = "Checking PhiI(., j)", marP=-c(1,1,1,0))
    cols <- c("red2", adjustcolor("blue", 0.25))
    for(j in 0:j.max) {
        curve(PhiI(x, j=j), -4, 10, col=cols[1], main = bquote(j == .(j)))
        if(j == j.max %/% 2)
            legend("right", c("PhiI()", "integrate(..)"),
                   col=cols, lwd = c(1,3), lty = c(1,3), inset = 1/40)
        I <- sapply(tt, function(t)
                    integrate(function(u) u^j * dnorm(u), -Inf, t)$value)
        lines(tt, I, col= cols[2], lwd=3, lty = 3)
    }
    par(oo$old.par)

}

## Huber:
huberPsi <- psiFunc(rho =
                  function(x, k) {
                      r <- u <- abs(x); I <- u < k
                      r[ I] <- u[I]^2 / 2
                      r[!I] <- k*(u[!I] - k / 2)
                      r
                  },
                  psi  = function(x, k) pmin.int(k, pmax.int(-k, x)),
                  wgt  = function(x, k) pmin.int(1, k/abs(x)),
                  Dpsi = function(x, k) abs(x) <= k,
                  Erho = function(k) {iP <- pnorm(k, lower=FALSE)
                                      1/2 - iP + k*(dnorm(k) - k*iP)},
                  Epsi2= function(k) ifelse(k < 10,
                  1 - 2*(k*dnorm(k) + (1-k*k)*pnorm(k, lower=FALSE)), 1),
                  EDpsi= function(k) 2*pnorm(k) - 1,
                  name = "Huber",
                  ## the tuning pars and default:
                  k = 1.345)

## Hampel:
hampelPsi <-
    psiFunc(rho = function(x, k)
        {
            u <- abs(x)
            a <- k[1] ; b <- k[2]; r <- k[3]
            Lg <- r <= u
            I <- u < a
            m1 <- !I & (I2 <- u < b)  # a <= u < b : 'constant'
            m2 <- !I2 & !Lg           # b <= u < r : 'descending'
            x[ I] <-  u[I]^2 / 2
            x[m1] <-  a*(a/2 + (u[m1] - a))
            ##x[m2]<- a*(a/2 + (b - a)) + a*(u^2 - b^2)/(2*(r - b))
            ##x[m2]<- a*(b - a/2)       + a*(u^2 - b^2)/(2*(r - b))
            x[m2] <-  a*(b - a/2 + (u[m2] - b)*(r - (b+u[m2])/2)/(r - b))
            ##u=r: a*(b - a/2 + (b + r)/2)
            x[Lg] <-  a/2*(b - a + r)
            x
        },
            psi = function(x, k)
        {
            ## this is "optimized" ==> factors faster than using ifelse()!
            u <- abs(x)
            lrg <- k[3] <= u
            mid <- k[1] < u & !lrg      # constant _and_ descending
            ## x is result for |x| < k[1]
            x[lrg] <- 0
            if(any(mid))
                x[mid] <- k[1] * sign(x[mid])*
                    pmin.int(1, (u[mid] - k[3])/(k[2] - k[3]))
            x
        },
            wgt = function(x, k)
        {
            x <- abs(x)
            lrg <- k[3] <= x
            I <- x < k[1]
            mid <- !I & !lrg            # contains constant and descending
            x[I] <- 1
            x[lrg] <- 0
            if(any(mid))
                x[mid] <- k[1] / x[mid] *
                    pmin.int(1, (x[mid] - k[3])/(k[2] - k[3]))
            x
        },
            Dpsi = function(x, k)
        {
            stopifnot(length(k) == 3, diff(k) >= 0) # for now
            u <- abs(x)
            lrg <- k[3] <= u
            I <- u < k[1]
            m1 <- !I & (I2 <- u < k[2]) # k_1 <= u < k_2: 'constant'
            m2 <- !I2 & !lrg            # k_2 <= u < k_3 : 'descending'
            x[lrg | m1] <- 0
            x[I ] <- 1
            x[m2] <- k[1] / (k[2] - k[3])
            x
        },

            Erho = function(k)
        {
	    names(k) <- c("a","b","r")
	    a <- k[["a"]] ; b <- k[["b"]]; r <- k[["r"]]
	    ph <- dnorm(k)
	    Ph <- pnorm(k)
	    ## rho(x) =	 c0   for  |x| >= r
	    c0 <- a/2*(b - a + r)
	    ## coeff. of rho(x) = a/2(c1 + c2|x| + c2 x^2), for |x| in [b,r]
	    D2 <- r - b
	    c1 <- -(a*r+ b*(b-a)) / D2
	    c2 <- 2*r / D2
	    c3 <- - 1 / D2
	    dPh.rb <- Ph[["r"]] - Ph[["b"]]
	    dph.rb <- ph[["r"]] - ph[["b"]]
	    ## Phi_2(r) - Phi_2(b) :=
	    dPh2.rb <- Ph[["r"]] - Ph[["b"]] - r*ph[["r"]] + b*ph[["b"]]
	    ## E[rho(X)] =
            ## [0,a] : 2* 1/2*(Phi_2(a)  - Phi_2(0))
            (Ph[["a"]]-a*ph[["a"]] - 1/2) +
            ## [a,b] : 2* a*( -a/2*(Phi(b) - Phi(a)) + (Phi_1(b) - Phi_1(a)) )
            2*a*(-a/2*(Ph[["b"]]-Ph[["a"]]) + (ph[["a"]] - ph[["b"]])) +
            ## the upper two can be simplified to
	    ## -1/2 + a*ph[["a"]] + (1+a^2)*Ph[["a"]] -2*a*ph[["b"]] - a^2*Ph[["b"]] +
            ## [b,r] :
		a*(c1*dPh.rb + c2*(-dph.rb) + c3*dPh2.rb) +
            ## [r,Inf] :
		    2*c0*(1 - Ph[["r"]])
        }

            ,
            Epsi2 = function(k) ## E[psi^2]=: 'A' in Hampel et al.(1986), p.150
        {
            names(k) <- c("a","b","r")
            a <- k[["a"]] ; r <- k[["r"]]
            ph <- dnorm(k)
            Ph <- pnorm(k)
            Ph2 <- Ph - k*ph # = Phi_2(k) {see PhiI(.) above}
            2*(Ph2[["a"]] - 1/2 + a^2*(Ph[["b"]] - Ph[["a"]]) +
               (a / (r - k[["b"]]))^2 * (
                   r^2 *(Ph[["r"]] - Ph[["b"]]) -2*r *(ph[["b"]] - ph[["r"]])
                   + Ph2[["r"]] - Ph2[["b"]]))
        },
            EDpsi= function(k) ## E[psi'] =: 'B' in Hampel et al.(1986)
        {
            a <- k[1] ; b <- k[2]; r <- k[3]
            2*(pnorm(a) - 1/2 - a* (pnorm(r) - pnorm(b)) / (r - b))
        },
            name = "Hampel",
            ## the tuning pars and default:
            k = c(2,4,8) / 1.345)# 1/1.345 = 0.7435

## TODO:  Biweight :
## ----   --------  but note that we have
##  (non-S4) ./biweight-funs.R  already {used by lmrob.*()}
##             ~~~~~~~~~~~~~~~
if(FALSE)
tukeyPsi <- c() ##########


## maybe TODO: Optimal tanh() estimator for location



### B.  M-Estimators of Scale --- need chi() and slightly different functionals
### --- ----------------------
###
## one "challenge" is the  a(b)  needed in  chi(x; a,b) = [x^2 -1 -a]_b^b
## for  V-optimal  M-Estimates of scale
## --> but that's solved (!) in ./scale-chi-opt.R
##                              ~~~~~~~~~~~~~~~~~
## Then, I'd also want the optimal chi for s
