library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Test of fitting random slopes
#
# Same data set as slope1, refit using variance matrices
#
set.seed(56)
n.subject <- seq(180, by=21, length=9) # number of subjects
slope <- sort(-.5 + rnorm(9, sd=.5))         # true treament effects

inst <- rep(1:9, n.subject)
n <- length(inst)
simdata <- data.frame(id=1:n, inst=inst,
                      trt= rep(0:1, length=n),
                      age= runif(n, 40, 70))
#risk goes up 30%/decade of age
simdata$hazard <- .8* exp(simdata$trt * rep(slope, n.subject) +
                          (simdata$age-55) * .03)

rtime <- function(hazard, censor=c(1,2)) {
    stime <- rexp(length(hazard), rate=hazard)
    ctime <- runif(length(hazard), censor[1], censor[2])
    list(time= pmin(stime, ctime), status=1*(stime <=ctime))
    }
temp <- rtime(simdata$hazard)
simdata$time <- temp$time
simdata$status <- temp$status

#Check fit1
fit1 <- coxme(Surv(time, status) ~ age + trt + (1|inst), simdata)
fit1b <- coxme(Surv(time, status) ~ age + trt + (1|inst), simdata,
               varlist=diag(9))
aeq(fit1$log, fit1b$log)
aeq(as.matrix(fit1$var), as.matrix(fit1b$var))
aeq(fixef(fit1), fixef(fit1b))

# Check fit2
#  To get the same iteration path you need to force the same
#  starting estimates
idlist <- sort(outer(1:9, 0:1, paste, sep='/'))
fit2 <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               varlist=coxmeFull(collapse=TRUE), vinit=c(.1, .1))

mat1 <- matrix(diag(18), 18, dimnames=list(idlist, idlist))
mat2 <-  bdsBlock(idlist, rep(1:9, each=2))
fit2b <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               varlist=list(mat1, mat2), vinit=c(.1, .1))
aeq(fit2$log, fit2b$log)
aeq(as.matrix(fit2$var), as.matrix(fit2b$var), tol=1e-7)
aeq(fixef(fit2), fixef(fit2b))

# Check fit3
# This takes two different paths to the solution, due to use of a
#  slope coefficient rather than nested.  So results are a bit
fit3 <- coxme(Surv(time, status) ~ age + trt + (1|inst) + (trt|inst),simdata,
              vinit=list(.1, .2))
mat3 <- diag(rep(0:1, 9))
dimnames(mat3) <- list(idlist, idlist)
fit3b <-  coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               varlist=coxmeMlist(list(mat2, mat3), rescale=F, pdcheck=F),
               vinit=c(.1, .2))

aeq(fit3$log, fit3b$log)
aeq(fixef(fit3), fixef(fit3b))
all.equal(unlist(VarCorr(fit3)), unlist(VarCorr(fit3b)), tolerance=1e-7, 
          check.attributes=FALSE)
# This function should map between coefficients
map <- matrix(0, 20,20)
for (i in 1:9) {
    map[i, 2*i -1] <- 1
    map[i+9, 2*i -1] <- -1
    map[i+9, 2*i] <- 1
    }
map[19,19] <- 1
map[20,20] <- 1

# Some of the random effects are very close to zero
aeq(unlist(ranef(fit3)), c(map[1:18,1:18] %*% unlist(ranef(fit3b))),
    tol=1e-7)
aeq(as.matrix(fit3$var), map %*% as.matrix(fit3b$var) %*% t(map),
    tol=1e-7)
