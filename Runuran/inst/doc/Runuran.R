### R code from vignette source 'Runuran.Rnw'

###################################################
### code chunk number 1: Runuran.Rnw:158-159
###################################################
library("Runuran")


###################################################
### code chunk number 2: specialgenerators.Rnw:40-45
###################################################
## Draw sample of size 10 from standard Gaussian distribution 
urnorm(10)

## Draw sample from truncated non-standard Gaussian distribution 
urnorm(10, mean = 1, sd = 0.5, lb = 2, ub = Inf)


###################################################
### code chunk number 3: specialgenerators.Rnw:52-53
###################################################
urnorm


###################################################
### code chunk number 4: universal.Rnw:24-38
###################################################
## Use method TDR (Transformed Density Rejection) to 
## draw a sample of size 10 from a hyperbolic distribution with PDF
##   f(x) = const * exp(-sqrt(1+x^2))  on domain (-Inf,Inf).

## We first have to define a function that returns the density.
pdf <- function (x) { exp(-sqrt(1+x^2)) }

## Next create the UNU.RAN object.
gen <- tdr.new(pdf=pdf, lb=-Inf, ub=Inf)

## Now we can use this object to draw the sample.
## (Of course we can repeat this step as often as required.)
x <- ur(gen,10)
x


###################################################
### code chunk number 5: universal.Rnw:45-46
###################################################
gen


###################################################
### code chunk number 6: universal.Rnw:54-65
###################################################
## Use method ARS (Adaptive Rejection Sampling) to 
## draw a sample of size 10 from a truncated Gaussian on [100,Inf).

## Define a function that returns the log-density.
lpdf <- function (x) { -0.5*x^2 }

## Create UNU.RAN object.
gen <- ars.new(logpdf=lpdf, lb=100, ub=Inf)

## Draw sample.
ur(gen,10)


###################################################
### code chunk number 7: universal.Rnw:70-78
###################################################
## Draw a sample from Gaussian distribution with 
## mean 2 and standard deviation 5.

## Create UNU.RAN object.
## Use R function 'dnorm(x, mean=2, sd=5, log=TRUE)' as density.
gen <- ars.new(logpdf=dnorm, lb=-Inf, ub=Inf,  mean=2, sd=5, log=TRUE)
## Draw sample.
ur(gen,10)


###################################################
### code chunk number 8: universal.Rnw:84-93
###################################################
## Compute quantiles for hyperbolic distribution with PDF
##   f(x) = const * exp(-sqrt(1+x^2))  on domain (-Inf,Inf).
## Thus we need an inversion method. We choose PINV.

## Create UNU.RAN object.
pdf <- function (x) { exp(-sqrt(1+x^2)) }
gen <- pinv.new(pdf=pdf, lb=0, ub=Inf, uresolution=1e-14)
## Get some quantiles
uq(gen, c(0.005, 0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99, 0.995))


###################################################
### code chunk number 9: universal.Rnw:100-118
###################################################
## Compute density for a given distribution or generator object.
## However, this only works when the density is already stored in 
## the object.

## Use distribution object
distr <- unuran.cont.new(pdf=function(x){exp(-x)}, lb=0,ub=Inf)
x <- ud(distr, 0:5)
x

## Use generator object
gen <- pinvd.new(distr)
x <- ud(gen, 0:5)
x

## Method PINV can also be used to estimate the CDF of the distribution
x <- up(gen, 0:5)
x



###################################################
### code chunk number 10: distributions.Rnw:22-33
###################################################
## Create an object for a gamma distribution with shape parameter 5.
distr <- udgamma(shape=5)

## Create the UNU.RAN generator object. use method PINV (inversion).
gen <- pinvd.new(distr)

## Draw a sample of size 100
x <- ur(gen,100)

## Compute some quantiles for Monte Carlo methods
x <- uq(gen, (1:9)/10)


###################################################
### code chunk number 11: advanced.Rnw:148-172
###################################################
## Use method TDR (Transformed Density Rejection) to 
## draw a sample of size 10 from a hyperbolic distribution with PDF
##   f(x) = const * exp(-sqrt(1+x^2)) 
## restricted to domain [-1,2].

## We first have to define functions that return the log-density and
## its derivative, respectively. (We also could use the density itself.)
lf  <- function (x) { -sqrt(1+x^2) }
dlf <- function (x) { -x/sqrt(1+x^2) }

## Next create the continuous distribution object.
d <- unuran.cont.new(pdf=lf, dpdf=dlf, islog=TRUE, lb=-1, ub=2,
                     name="hyperbolic")

## Create UNU.RAN object. We choose method TDR with 
## immediate acceptance (IA) and parameter c=0.
gen <- unuran.new(distr=d, method="tdr; variant_ia; c=0")

## Now we can use this object to draw the sample.
## (Of course we can repeat this step as often as required.)
ur(gen,10)

## Here is some information about our generator object.
unuran.details(gen)


###################################################
### code chunk number 12: advanced.Rnw:180-195
###################################################
## Use method DGT (Discrete Guide Table method) to 
## draw a sample of size 10 from a Binomial distribution given
## its probability vector.

## Create instances of a discrete distribution object
d <- unuran.discr.new(pv=dbinom(0:100,100,0.4), lb=0, name="binomial(100,0.4)")

## Create UNU.RAN object. We choose method DGT.
gen <- unuran.new(distr=d, method="dgt")

## Now we can use this object to draw the sample.
ur(gen,10)

## Here is some information about our generator object.
unuran.details(gen)


###################################################
### code chunk number 13: advanced.Rnw:204-222
###################################################
## Use method DSROU (Discrete Simple Ratio-Of-Uniforms method) to 
## draw a sample of size 10 from a discrete distribution with
## given PMF, mode, and sum.

## Define functions that return the PMF.
f  <- function (x) { 0.4 * (1-0.4)^x }

## Create the continuous distribution object.
d <- unuran.discr.new(pmf=f, lb=0, ub=Inf, mode=0, sum=1)

## Create UNU.RAN object. We choose method DARI with squeezes.
gen <- unuran.new(distr=d, method="dari; squeeze=on")

## Now we can use this object to draw the sample.
ur(gen,10)

## Here is some information about our generator object.
unuran.details(gen)


###################################################
### code chunk number 14: advanced.Rnw:230-249
###################################################
## Use method VNROU (Multivariate Naive Ratio-Of-Uniforms) to 
## draw a sample of size 5 from a bivariate distribution
## with given PDF, mode and domain.

## Define functions that return the PDF.
f  <- function (x) { exp(-sum(x^4)) }

## Create the continuous distribution object.
d <- unuran.cmv.new(dim=2, pdf=f, mode=c(0,0), ll=c(-1,-1), ur=c(1,1),
                    name="bivariate power-exponential")

## Create UNU.RAN object. We choose method VNROU with parameter r=0.5.
gen <- unuran.new(distr=d, method="vnrou; r=0.5")

## Now we can use this object to draw the sample.
ur(gen,5)

## Here is some information about our generator object.
unuran.details(gen)


###################################################
### code chunk number 15: advanced.Rnw:259-263 (eval = FALSE)
###################################################
## ## Try to use method TDR with missing data.
## lf  <- function (x) { -sqrt(1+x^2) }
## d <- unuran.cont.new(pdf=lf, lb=-Inf, ub=Inf, islog=TRUE)
## gen <- unuran.new(distr=d, method="tdr")


###################################################
### code chunk number 16: fig:rejection
###################################################
x <- c((0:31)/10,3.1416)
y <- round(sin(c((0:31)/10,3.1416)),3)
pdf.sin <- paste("(",paste(x,y,sep=",",collapse=") ("),")",sep="",collapse="")
rm(x,y)


###################################################
### code chunk number 17: pitfalls.Rnw:47-51
###################################################
pdf <- function (x) { x^2 / (1+x^2)^2 }
gen <- pinv.new(pdf=pdf,lb=0,ub=Inf,  center=1 )    ## Add 'center'
x <- ur(gen,10)
x


###################################################
### code chunk number 18: pitfalls.Rnw:74-86
###################################################
## create a unuran object using method 'PINV'
gen <- pinv.new(dnorm,lb=0,ub=Inf)

## such an object can be packed
unuran.packed(gen) <- TRUE

## it can be still used to draw a random sample
x <- ur(gen,10)
x

## we also can check whether a unuran object is packed
unuran.packed(gen)


