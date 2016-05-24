stopifnot(require("robustbase"))
## testing functions:
source(system.file("test-tools-1.R",  package = "Matrix", mustWork=TRUE))# assert.EQ

c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
showProc.time <- local({ ## function + 'pct' variable
    pct <- proc.time()
    function(final="\n") { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	## 'Time ..' *not* to be translated:  tools::Rdiff() skips its lines!
	cat('Time elapsed: ', (pct - ot)[1:3], final)
    }
})

## as long as we don't export these (nor provide an nlrob(., method=.) interface:
nlrob.MM  <- robustbase:::nlrob.MM
nlrob.tau <- robustbase:::nlrob.tau
nlrob.CM  <- robustbase:::nlrob.CM
nlrob.mtl <- robustbase:::nlrob.mtl

(doExtras <- robustbase:::doExtras())
if(doExtras) {
    NP <- 30 ; tol <- 1e-11
} else { ## "fast"
    NP <- 15 ; tol <- 1e-7
}

start.from.true <- !doExtras # (but not necessarily ..)
if(start.from.true) { # population size = NP (random) + 1 (true parameters)
    init_p       <- c(1, 0.2)
    init_p_sigma <- c(1, 0.2, 1)
} else {
    init_p <- init_p_sigma <- NULL
}

if(!dev.interactive(orNone=TRUE))  pdf("nlregrob-tst.pdf")

## Stromberg, Arnold J. (1993).
## Computation of high breakdown nonlinear regression parameters.
## J. Amer. Statist. Assoc. 88(421), 237-244.

## exponential regression
Expo <- function(x, a, b) exp(a + b*x)
set.seed(2345) # for reproducibility
## data without outliers:
d.exp30 <- data.frame(x = sort( runif(30, 0, 10) ), err = rnorm(30))
d.exp30 <- transform(d.exp30, y = Expo(x, 1, 0.2) + err)
## classical (starting at truth .. hmm)
Cfit <- nls(y ~ Expo(x, a, b), data = d.exp30, start = c(a = 1, b = 0.2),
            control = nls.control(tol = 8e-8, printEval = TRUE))
showProc.time()#                        ---- OS X needing 6e-8

## robust
Rfit.MM.S.bisquare <-
    nlrob.MM(y ~ Expo(x, a, b), data = d.exp30,
             lower = c(a = -10, b = -2), upper = c(10, 2),
             NP = NP, tol = tol, add_to_init_pop = init_p )
if(doExtras) {
Rfit.MM.S.lqq        <- update(Rfit.MM.S.bisquare, psi = "lqq")
Rfit.MM.S.optimal    <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.MM.S.hampel     <- update(Rfit.MM.S.bisquare, psi = "hampel")
}
showProc.time()
Rfit.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, init = "lts")
Rfit.MM.lts.lqq      <- update(Rfit.MM.S.bisquare, init = "lts", psi = "lqq")
Rfit.MM.lts.optimal  <- update(Rfit.MM.S.bisquare, init = "lts", psi = "optimal")
Rfit.MM.lts.hampel   <- update(Rfit.MM.S.bisquare, init = "lts", psi = "hampel")
showProc.time()

S.time(Rfit.tau.bisquare <-
    nlrob.tau( y ~ Expo(x, a, b), data = d.exp30,
               lower = c(a = -10, b = -2), upper = c(10, 2),
               NP = NP, add_to_init_pop = init_p ))
S.time(Rfit.tau.optimal <- update(Rfit.tau.bisquare, psi = "optimal"))

S.time(Rfit.CM <- nlrob.CM( y ~ Expo(x, a, b), data = d.exp30,
			    lower = c(a = -10, b = -2, sigma = 0),
			    upper = c(	   10,	    2,	      10),
                            NP = NP, add_to_init_pop = init_p_sigma ))
S.time(Rfit.mtl <- nlrob.mtl(y ~ Expo(x, a, b), data = d.exp30,
			     lower = c(a = -10, b = -2, sigma = 0),
			     upper = c(	    10,	     2,		3),
			     NP = NP+10, # <- higher prob. to get close
                             tol = tol,
                             trace=TRUE, details=TRUE,
                             add_to_init_pop = init_p_sigma ))
showProc.time()

plot(y ~ x, d.exp30, main = "Data = d.exp30")
cTr <- adjustcolor("red4", 0.5)
cLS <- adjustcolor("blue2", 0.5)
cE <- curve(Expo(x, a=1, b=0.2), 0, 10, n=1+2^9, col=cTr, lwd=2, lty=2, add=TRUE)
lines(d.exp30$x, fitted(Cfit), col=cLS, lwd=3)
ll <- length(m1 <- sapply(ls.str(patt="^Rfit"), get, simplify=FALSE))
.tmp <- lapply(m1, function(.) lines(d.exp30$x, fitted(.)))
legend("topleft", c("true", "LS", names(m1)),
       lwd=c(2,3, rep(1,ll)), lty=c(2,1, rep(1,ll)),
       col=c(cTr,cLS, rep(par("fg"),ll)), bty="n", inset=.01)
showProc.time()

## 40% outliers present {use different data name: seen in print(<fitted model>)
d.exp40out <- within(d.exp30, y[15:27] <- y[15:27] + 100)
Cfit.40out  <- update(Cfit, data = d.exp40out,
                      control = nls.control(tol = Cfit$control$tol))
Cfit.no.out <- update(Cfit.40out, subset = -(15:27))

if(doExtras) {
Rf.out.MM.S.bisquare   <- update(Rfit.MM.S.bisquare, data=d.exp40out)
Rf.out.MM.S.lqq        <- update(Rf.out.MM.S.bisquare, psi = "lqq")
Rf.out.MM.S.optimal    <- update(Rf.out.MM.S.bisquare, psi = "optimal")
Rf.out.MM.S.hampel     <- update(Rf.out.MM.S.bisquare, psi = "hampel")
showProc.time()
}
Rf.out.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, data=d.exp40out, init= "lts")
Rf.out.MM.lts.lqq      <- update(Rf.out.MM.lts.bisquare, psi= "lqq") #-----------
Rf.out.MM.lts.optimal  <- update(Rf.out.MM.lts.bisquare, psi= "optimal")
Rf.out.MM.lts.hampel   <- update(Rf.out.MM.lts.bisquare, psi= "hampel")
showProc.time()

Rf.out.tau.bisquare <- update(Rfit.tau.bisquare, data=d.exp40out)
Rf.out.tau.optimal  <- update(Rfit.tau.bisquare, data=d.exp40out, psi = "optimal")
Rf.out.CM  <- update(Rfit.CM,  data=d.exp40out)
Rf.out.mtl <- update(Rfit.mtl, data=d.exp40out)
showProc.time()

plot(y ~ x, d.exp40out, main = "Data = d.exp40out")
cE <- curve(Expo(x, a=1, b=0.2), 0, 10, n=1+2^9, col=cTr, lwd=2, lty=2, add=TRUE)
ll <- length(m1 <- sapply(ls.str(patt="^Rf.out"), get, simplify=FALSE))
.tmp <- lapply(m1, function(.) lines(d.exp40out$x, fitted(.)))
xx <- local({p <- par("usr"); seq(p[1],p[2], length.out=256)})
lines(xx, predict(Cfit.no.out, list(x=xx)), col=cLS, lwd=3)
lines(xx, predict(Cfit.40out , list(x=xx)), col=cLS, lty=2)
legend("topleft", c("true", "LS [w/o outl]", "LS", names(m1)),
       lwd=c(2,3, rep(1,1+ll)), lty=c(2,1,2, rep(1,ll)),
       col=c(cTr,cLS,cLS, rep(par("fg"),ll)), bty="n", inset=.01)
showProc.time()

## presence of high leverage point outliers
d.exp.Hlev <- within(d.exp40out, {
    x[28:30] <- x[28:30] + 10   ## shift  10
    y <- Expo(x, 1, 0.2) + err
    y[28:30] <- y[28:30] + 500
})
Cfit.Hlev <- update(Cfit.40out, data = d.exp.Hlev)
Cfit.no.Hlev <- update(Cfit.Hlev, subset = -(28:30))
showProc.time()

if(doExtras) {
Rf.Hlev.MM.S.bisquare   <- update(Rfit.MM.S.bisquare, data = d.exp.Hlev)
Rf.Hlev.MM.S.lqq        <- update(Rf.Hlev.MM.S.bisquare, psi = "lqq")
Rf.Hlev.MM.S.optimal    <- update(Rf.Hlev.MM.S.bisquare, psi = "optimal")
Rf.Hlev.MM.S.hampel     <- update(Rf.Hlev.MM.S.bisquare, psi = "hampel")
showProc.time()
}
Rf.Hlev.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, data = d.exp.Hlev, init="lts")
Rf.Hlev.MM.lts.lqq      <- update(Rf.Hlev.MM.lts.bisquare, psi= "lqq")
Rf.Hlev.MM.lts.optimal  <- update(Rf.Hlev.MM.lts.bisquare, psi="optimal")
Rf.Hlev.MM.lts.hampel   <- update(Rf.Hlev.MM.lts.bisquare, psi= "hampel")
showProc.time()

Rf.Hlev.tau.bisquare <- update(Rfit.tau.bisquare, data = d.exp.Hlev)
Rf.Hlev.tau.optimal  <- update(Rf.Hlev.tau.bisquare, psi = "optimal")
Rf.Hlev.CM  <- update(Rfit.CM,  data = d.exp.Hlev)
Rf.Hlev.mtl <- update(Rfit.mtl, data = d.exp.Hlev)
showProc.time()

plot(y ~ x, d.exp.Hlev, main = "Data = d.exp.Hlev")
cE <- curve(Expo(x, a=1, b=0.2), 0, par("usr")[2], n=1+2^9, col=cTr, lwd=2, lty=2, add=TRUE)
x.H <- seq(par("usr")[1], par("usr")[2], length.out = 256)
ll <- length(m1 <- sapply(ls.str(patt="^Rf.Hlev"), get, simplify=FALSE))
.tmp <- lapply(m1, function(.) lines(x.H, predict(., list(x=x.H))))
lines(x.H, predict(Cfit.no.Hlev, list(x=x.H)), col=cLS, lwd=3)## L.S.(<good data>)
lines(x.H, predict(Cfit.Hlev,    list(x=x.H)), col=cLS, lty=2)## L.S.
legend("topleft", c("true", "LS [w/o outl]", "LS", names(m1)),
       lwd=c(2,3, rep(1,1+ll)), lty=c(2,1,2, rep(1,ll)),
       col=c(cTr, cLS,cLS, rep(par("fg"),ll)), bty="n", inset=.01)
showProc.time()

				        cfcl <- coef(Cfit)
				        cfcl.n.o <- coef(Cfit.no.out)
				        cfcl.n.H <- coef(Cfit.no.Hlev)
## no outliers present
assert.EQ(coef(Rfit.MM.S.bisquare),	cfcl, tol = 0.01, giveRE=TRUE)
if(doExtras) {
assert.EQ(coef(Rfit.MM.S.lqq),		cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.S.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.S.hampel),	cfcl, tol = 0.01, giveRE=TRUE)
}
assert.EQ(coef(Rfit.MM.lts.bisquare),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.lqq),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.hampel),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.tau.bisquare),	cfcl, tol = 0.02, giveRE=TRUE)# 0.009873
assert.EQ(coef(Rfit.tau.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.CM)[-3],		cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.mtl)[-3],		cfcl, tol = 0.02, giveRE=TRUE)
## 40% outliers present -- compare with L.S.(good.data)
if(doExtras) {
assert.EQ(coef(Rf.out.MM.S.bisquare),	cfcl.n.o, tol = 7e-4, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.lqq),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.optimal),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.hampel),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
}
assert.EQ(coef(Rf.out.MM.lts.bisquare),	cfcl.n.o, tol = 6e-4, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.lqq),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.optimal),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.hampel),	cfcl.n.o, tol = 1e-5, giveRE=TRUE)
assert.EQ(coef(Rf.out.tau.bisquare),	cfcl.n.o, tol = .007, giveRE=TRUE)
assert.EQ(coef(Rf.out.tau.optimal),	cfcl.n.o, tol = .002, giveRE=TRUE)
assert.EQ(coef(Rf.out.CM)[-3],		cfcl.n.o, tol = .012, giveRE=TRUE)# 0.00708,0.01079
assert.EQ(coef(Rf.out.mtl)[-3],		cfcl.n.o, tol = .002, giveRE=TRUE)# better in 64b
## presence of high leverage point outliers -- compare with LS(good.data)
if(doExtras) {
assert.EQ(coef(Rf.Hlev.MM.S.bisquare),	cfcl.n.H, tol = .01,  giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.lqq),	cfcl.n.H, tol = .02,  giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.optimal),	cfcl.n.H, tol = .005, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.hampel),	cfcl.n.H, tol = .02,  giveRE=TRUE)
}
assert.EQ(coef(Rf.Hlev.MM.lts.bisquare),cfcl.n.H, tol = .01,  giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.lqq),	cfcl.n.H, tol = .015, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.optimal), cfcl.n.H, tol = .002, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.hampel),	cfcl.n.H, tol = .02,  giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.tau.bisquare),	cfcl.n.H, tol = .05,  giveRE=TRUE)# 0.0363, 0.0415
assert.EQ(coef(Rf.Hlev.tau.optimal),	cfcl.n.H, tol = .03,  giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.CM)[-3],		cfcl.n.H, tol = .12,  giveRE=TRUE)# 0.032, 0.082
assert.EQ(coef(Rf.Hlev.mtl)[-3],	cfcl.n.H, tol = .08,  giveRE=TRUE)

length(mods <- sapply(ls.str(patt="^Rf"), get, simplify=FALSE)) # 36
is.conv <- sapply(mods, `[[`, "status") == "converged"
prblm <- mods[!is.conv]
if(length(prblm)) {
    cat("\n*** NON-converged model fits:\n")
    print(prblm)
    mods <- mods[is.conv]
} else cat("\n All models converged\n")

## Now, all mods are converged  -----------

dKnd <- as.factor(vapply(mods, function(.m.)
                         as.character(getCall(.m.)[["data"]]), ""))
table(dKnd) ##
(iKnd <- setNames(seq_len(nlevels(dKnd)), levels(dKnd)))

## Coefficients: Some have 'sigma', some not:
pcf <- vapply(lcf <- lapply(mods, coef),  length, 1)
table(pcf) ## 2 and 3
stopifnot(min(pcf) + 1 == max(pcf)) # +1 : those which have 'sigma
pp <- min(pcf)
ccf <- t(simplify2array(lapply(lcf, `[`, 1:max(pcf))))
## take the "Scale" for those that do not have 'sigma' among  coef():
i.n <- is.na(ccf[,"sigma"])
ccf[i.n, "sigma"] <- vapply(mods[i.n], `[[`, 0, "Scale")
    ## not yet: vapply(mods[i.n], sigma, 0.)
ccf
## well, the  'sigma's  are definitely *not* unbiased estimates of
## true sqrt(var(eps))  ...  [FIXME]
## --> indeed, this can be found in the  CM  paper [TODO: write more here]

plot(ccf[,1:2], pch = as.integer(dKnd))## use 'method' to get color
legend("topright", inset=.01, names(iKnd), pch = iKnd)
points(rbind(cfcl.n.H, cfcl, cfcl.n.o), # <- order from iKind
       col=adjustcolor("tomato",.5), cex=2, pch=1:3, lwd=5)
## optional
labs <- sub("^Rfit\\.", '', sub("^Rf\\.[A-Za-z]+\\.", '', rownames(ccf)))
labs <- sub("hampel$", "Ham", sub("optimal$", "opt", sub("bisquare$", "biS", labs)))
labs
text(ccf[,1:2], labs, cex=0.75, col=adjustcolor(1, 0.5),
     adj= -1/5, srt=75, xpd=NA)
points(rbind(cfcl), col=adjustcolor("tomato",.5), cex=2, pch=3, lwd=5)
showProc.time()


###------- Extended Tests for the DNase1 example from >>>> ../man/nlrob-algos.Rd <<<<
###							   =====================
DNase1 <- DNase[DNase$Run == 1,]
form <- density ~ Asym/(1 + exp(( xmid -log(conc) )/scal ))
pnms <- c("Asym", "xmid", "scal")
psNms <- c(pnms, "sigma")
set.seed(47) # as these by default use randomized optimization:

fMM <- robustbase:::nlrob.MM(form, data = DNase1,
           lower = setNames(c(0,0,0), pnms), upper = 3,
           ## call to nlrob.control to pass 'optim.control':
           ctrl = nlrob.control("MM", optim.control = list(trace = 1),
                                optArgs = list(trace = TRUE)))
showProc.time()

## for comparisons, later:
all.eq.mod <- function(m1, m2, sub=FALSE, excl = c("call", "ctrl"), ...) {
    nm1 <- names(m1)
    stopifnot(if(sub) nm1 %in% names(m2) else nm1 == names(m2))
    ni <- if(sub)
	      nm1[is.na(match(nm1, c("call","ctrl")))]
	  else is.na(match(names(m1), excl))## <<- all but those with names in 'excl'
    all.equal(m1[ni], m2[ni], ...)
}
if(doExtras) {## the same, with 'pnames' and unnamed 'lower':
    set.seed(47)
    tools::assertWarning(
        fM2 <- robustbase:::nlrob.MM(form, data = DNase1, pnames = pnms,
                                     lower = 0, upper = 3))
    stopifnot(all.eq.mod(fMM, fM2, tol=1e-15))

    ftau <- robustbase:::nlrob.tau(form, data = DNase1,
                                   lower= setNames(c(0,0,0), pnms), upper= 3,  trace=TRUE)

    fCM  <- robustbase:::nlrob.CM (form, data = DNase1,
                                   lower= setNames(c(0,0,0,0), psNms), upper= 3, trace=TRUE)

    fmtl <- robustbase:::nlrob.mtl(form, data = DNase1,
                                   lower= setNames(c(0,0,0,0), psNms), upper= 3, trace=TRUE)

    mods <- list(MM=fMM, tau=ftau, CM=fCM, MTL=fmtl)
    print(sts <- sapply(mods, `[[`, "status"))
    stopifnot(sts == "converged")

    print(sapply(mods, `[[`, "data"))   # currently 'language' %% FIXME

    print(sapply(mods, `[[`, "coefficients")) # nice matrix
    showProc.time()
}
## Compare with traditional M-estimate, a) started robustly b) psi = Tukey's:
fM <- nlrob(formula(fMM), data=eval(fMM$data), start = coef(fMM),
            psi = .Mwgt.psi1("bisquare"), trace = TRUE)
rbind(M=coef(fM), MM=coef(fMM)) # "similar" ... well, no: the sigma's get much different
## stopifnot(%%____FIXME___
all.equal(coef(fM), coef(fMM), tolerance = 1e-4)
## ) # had 3.26e-5
## FIXME:  nlrob( "M")  should allow to keep specify an initial sigma *and* keep that fixed
showProc.time()


### Now call the above methods via nlrob():
set.seed(47) # (same as above)
## without "sigma"
gMM  <- nlrob(form, data = DNase1, method = "MM",
              lower = setNames(c(0,0,0), pnms), upper = 3)
gtau <- nlrob(form, data = DNase1, method = "tau",
              lower = setNames(c(0,0,0), pnms), upper = 3)
## those with "sigma" -> must be in (lower, upper), too :
gCM  <- nlrob(form, data = DNase1, method = "CM",
              lower = setNames(c(0,0,0,0), psNms), upper = 3)
gmtl  <- nlrob(form, data = DNase1, method = "mtl",
              lower = setNames(c(0,0,0,0), psNms), upper = 3)
showProc.time()

## list {and test print(<nlrob>) for these}:
(mod2 <- list(MM=gMM, tau=gtau, CM=gCM, MTL=gmtl))

if(doExtras) {
    stopifnot(mapply(all.eq.mod, mods, mod2, sub=TRUE))
}
