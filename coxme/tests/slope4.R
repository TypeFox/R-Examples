#Check fit4
#
# The first model has random effects of 1 + trt | inst.
# The second form is 1 | inst/trt.  The test is to make them
#   the same using variance matrices for the second.
#
# The first has 18 dummy variables i1-i9 and t1-t9, the second has
#   g1-g18, one per group.
# Arithmetic on the dummys shows that g1= i1-t1, g2 = t1, g3=i2-t2, g4=t2,
#   etc.  So 
#    var(g1) =     var(i1) + var(t1) - 2*cov(i1,t1)
#    var(g2) =     var(t1)
#    cov(g1,g2) =  cov(i1,t1) - var(t1)
#
# So if we model g, we get a covariance that is a set of 2x2 blocks
#        a+b-c  c-b
#        c-b    b
# which can be written as a*mat1 + b*mat2 + c*mat3 for 3 constructed matrices
#
library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Same data set as slope1
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

contr.none <- function(n,contrasts=T) {
        if(is.numeric(n) && length(n) == 1.)
                levs <- 1.:n
        else {
                levs <- n
                n <- length(n)
        }
        contr <- array(0., c(n, n), list(levs, levs))
        contr[seq(1., n^2., n + 1.)] <- 1.
        contr
        }
options(contrasts=c('contr.none', 'contr.poly'))
igchol <- function(x) {
    dd <- diag(x)
    ll <- as.matrix(x)
    ll %*% diag(dd) %*% t(ll)
    }

# The basic building blocks
tempx <- matrix(0., nrow(simdata), ncol=18)
for (i in 1:9) {
    tempx[,i] <- 1*(simdata$inst==i)
    tempx[,i+9] <- tempx[,i] * simdata$trt
    }
cox1<- coxph(Surv(time, status) ~ tempx + age + trt, simdata,
             iter=0, x=T)
dt1 <- coxph.detail(cox1)
u1 <- apply(dt1$score, 2, sum)
imat1 <- apply(dt1$imat, 1:2, sum)

group <- strata(simdata$inst, simdata$trt, shortlabel=T, sep='/')
cox2 <- coxph(Surv(time, status) ~ group + age + trt, simdata,
              iter=0, x=T)
dt2 <- coxph.detail(cox2)
u2 <- apply(dt2$score, 2, sum)
imat2 <- apply(dt2$imat, 1:2, sum)

map <- matrix(0, 20,20)
for (i in 1:9) {
    map[i, 2*i -1] <- 1
    map[i+9, 2*i -1] <- -1
    map[i+9, 2*i] <- 1
    }
map[19,19] <- 1
map[20,20] <- 1
imap <- round(solve(map))  #inverse map is integers as well

aeq(u1 %*% map, u2)
aeq(u1, u2 %*% imap)
aeq(t(map) %*% imat1 %*% map, imat2)
aeq(imat1, t(imap) %*% imat2 %*% imap)


# Now for a fit of the first form of the model
fit0a <- coxme(Surv(time, status) ~ age + trt + (1 +trt |inst), simdata,
               iter=0, vfixed=c(.2, .1, .3))

aeq(fit0a$u, u1)

pen1 <- matrix(0., 20, 20)
pen1[cbind(1:9, 1:9)] <- .2
pen1[cbind(10:18, 10:18)] <- .3
pen1[cbind(1:9, 10:18)] <- .1*sqrt(.2* .3)
pen1[cbind(10:18, 1:9)] <- .1*sqrt(.2* .3)
ipen1 <- solve(gchol(pen1)) #generalized inverse

aeq(imat1 + ipen1, igchol(fit0a$hmat))
step1 <- solve(fit0a$hmat, fit0a$u)

# iteration 1
fit1a <- coxme(Surv(time, status) ~ age + trt + (1 +trt |inst), simdata,
               iter=1, vfixed=c(.2, .1, .3))
aeq(step1, c(unlist(ranef(fit1a)), fixef(fit1a)))

cox1.1<- coxph(Surv(time, status) ~ tempx + age + trt, simdata,
             iter=0, x=T, init=step1)

dt1.1 <- coxph.detail(cox1.1)
aeq(apply(dt1.1$score, 2, sum) - step1 %*% ipen1, fit1a$u)
aeq(apply(dt1.1$imat, 1:2,sum) + ipen1, igchol(fit1a$hmat))
step2 <- solve(fit1a$hmat, fit1a$u)

#iteration 2
fit2a <- coxme(Surv(time, status) ~ age + trt + (1 +trt |inst), simdata,
               iter=2, vfixed=c(.2, .1, .3))
aeq(step1+step2, c(unlist(ranef(fit2a)), fixef(fit2a)))



# Now try the solution using method #2
vfix <- c(.2, .3, .1*sqrt(.2*.3))  #var, var, covar
gname <- levels(group)

mat1 <- bdsmatrix(rep(c(1,0,0,0), 9), blocksize=rep(2,9),
                  dimnames=list(gname, gname))
mat2 <- bdsmatrix(rep(c(1,-1,-1,1), 9), blocksize=rep(2,9),
                  dimnames=list(gname, gname))
mat3 <- bdsmatrix(rep(c(-2,1,1,0), 9), blocksize=rep(2,9),
                  dimnames=list(gname, gname))
mat1 <- as.matrix(mat1)  #force non-sparse

pen2 <- matrix(0, 20, 20)
pen2[1:18, 1:18]  <- vfix[1]*mat1 + vfix[2]*as.matrix(mat2) + 
                    vfix[3]*as.matrix(mat3)
ipen2 <- solve(gchol(pen2))
aeq( t(map) %*% pen1 %*% map, pen2)
aeq( ipen1, map %*% ipen2 %*% t(map))

fit0b <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               iter=0, vfixed=vfix,
              varlist=coxmeMlist(list(mat1,mat2,mat3), pdcheck=F, rescale=F))
               
aeq(u2, fit0b$u)
aeq(imat2 + ipen2, igchol(fit0b$hmat))
step1b <- solve(fit0b$hmat, fit0b$u)

# Note that step1b != step1 %*% map, or t(map) or imap or t(imap),
#  since U and imat transform by map, but ipen by inverse-map = imap

# Iteration 1
fit1b <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               iter=1, vfixed=vfix,
              varlist=coxmeMlist(list(mat1,mat2,mat3), pdcheck=F, rescale=F))
aeq(step1b, c(unlist(ranef(fit1b)), fixef(fit1b)))


# And now the full fit
#  After more looks at this data, I realized that the likelihood
# wrt the correlation term is very very flat. Tiny changes in
# starting point or solution path move it a lot, but not the final
# loglik.
fita <- coxme(Surv(time, status) ~ age + trt + (1+trt | inst), simdata)
 
# fitb has a hard time, by way of wandering into bad solutions
# vtemp <- VarCorr(fita)[[1]][c(1,4,2)]
# vtemp[3] <- vtemp[3] * sqrt(vtemp[1] * vtemp[2])
#fitb <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
#              vinit=vtemp,
#              varlist=coxmeMlist(list(mat1,mat2,mat3), pdcheck=F, rescale=F,
#                                  positive=F))

# So create our own variance function, which expects v1, v2, cor
#  and matrices that are already in the right form
myvar <- function(varlist) {
    varlist <- varlist
    init <- function(vinit, fixed, intercept, G, X, sparse, ...) {
        ngroup <- length(G)
        n <- nrow(G)
    
        if (length(vinit) ==3) theta <- vinit 
        else theta <- c(.1, .2, 0)
        theta[3] <- (1+theta[3])/(1-theta[3])
        
    
        G <- rev(expand.nested(G))
        imap <- as.matrix(as.numeric(G[,1]))
        rname <- names(G)[1]
        bname <- levels(G[[1]])

        list(theta=log(theta), imap=imap, X=NULL, xmap=NULL, 
             parms=list(varlist=varlist, bname=bname, rname=rname,
                        vname=names(varlist)))
        }
    generate <- function(newtheta, parms) {
        theta <- exp(newtheta)
        theta[3] <- (theta[3]-1)/(theta[3] +1)  #correlation
        theta[3] <- theta[3] * sqrt(theta[1]*theta[2]) #covar
        
        theta[1]*parms$varlist[[1]] + theta[2]*parms$varlist[[2]] +
            theta[3] * parms$varlist[[3]]
        }
    wrapup <- function(newtheta, b, parms) {
        theta <- exp(newtheta)
        theta[3] <- (theta[3]-1)/(theta[3] +1)  #correlation
        
        mname <- c("Intercept","Slope")
        tmat <- matrix(theta[c(1,3,3,2)], 2,
                       dimnames=list(mname, mname))
        tmat <- list(tmat)
        names(tmat) <- "Inst"
        list(theta=tmat, b=b)
        }
    out <- list(initialize=init, generate=generate, wrapup=wrapup)
    class(out) <- 'coxmevar'
    out
    }

fitc <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
              varlist=myvar(list(mat1, mat2, mat3)))

aeq(fitc$log, fita$log, tol=1e-5)
aeq(fixef(fita), fixef(fitc), tol=1e-4)

vtemp <- VarCorr(fita)[[1]][c(1,4,2)]
fitc2 <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
              vfixed=vtemp, varlist=myvar(list(mat1, mat2, mat3)))
aeq(unlist(ranef(fita)),  map[1:18, 1:18] %*% unlist(ranef(fitc2)), tol=1e5)
aeq(fitc2$log, fita$log, tol=1e-5)
