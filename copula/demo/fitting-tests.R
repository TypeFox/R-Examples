#### Testing fitCopula
######################################

require(copula)

source(system.file("Rsource", "tstFit-fn.R", package="copula", mustWork=TRUE))
## ../inst/Rsource/tstFit-fn.R

options(nwarnings = 2000)# default 50 - only keeps the 50 last warnings..


## test code
system.time(
rr <- tstFit1cop(normalCopula(), tau.set=seq(0.2,0.8,by=0.2), n.set=c(25,50,100,200), N=200)
)
## ~ 400 seconds
## well, now (2012-09-17, lynne) faster:
##    user  system elapsed 
## 201.006   0.167 203.282 

d <- reshape.tstFit(rr)
plots.tstFit(d)# MM: log-scale desirable but ugly tick labelling

plots.tstFit(d, log=FALSE)


## t-copula instead of normal -- small set for testing here:
set.seed(17)
system.time(
rrt <- tstFit1cop(tCopula(df.fixed=TRUE),
                  tau.set=c(.4, .8),
                  n.set=c(10, 50, 200), N=128)
)
plots.tstFit(reshape.tstFit(rrt))


## Fit a tevCopula() :... still somewhat frequent optim errors(),
## from non-finite loglikCopula() values:
set.seed(3)
rf <- replFitCop(tevCopula(.6, df.fixed=TRUE), 
                 n = 25, N = 40, estimate.variance=FALSE)
warnings() # 11 warnings (out of N = 40)
##  In .local(copula, tau, ...) : tau is out of the range [0, 1]
##  In .local(copula, rho, ...) : rho is out of the range [0, 1]

##
set.seed(321)
system.time(
rtev <- tstFit1cop(tevCopula(df.fixed=TRUE),
                   tau.set= seq(0.2,0.8, by=0.2),
                   n.set = c(25,50,100,200), N=200,
                   estimate.variance = FALSE)##- not implemented, as:
    ## there is no formula for derPdfWrtParam*() for this copula
)
##     user   system  elapsed 
## 1279.605    0.625 1290.210 
## There were 50 or more warnings (use warnings() to see the first 50)
## an all are either  irho or itau related:
## 47: In .local(copula, rho, ...) : rho is out of the range [0, 1]
## 50: In .local(copula, tau, ...) : tau is out of the range [0, 1]
warnings()
plots.tstFit(reshape.tstFit(rtev))
## "mpl" often not plotted {"NA"}

proc.time()
