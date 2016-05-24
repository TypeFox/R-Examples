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


#Hammer through fit3 an iteration at a time
# Iteration 0, first form
cox1 <- coxph(Surv(time, status) ~ factor(inst)*trt + age, simdata, 
              iter=0, x=T)
indx <- c(1:9, 12:20, 11,10) #order of variables in coxph
cx1 <- cox1$x[,indx]
cox1 <- coxph(Surv(time, status) ~ cx1, simdata, iter=0, x=T)
dt1 <- coxph.detail(cox1)
u1 <- apply(dt1$score, 2, sum)
imat1 <- apply(dt1$imat,1:2, sum)

fit3a <- coxme(Surv(time, status) ~ age + trt + (1|inst) + (trt|inst),
               simdata, vfixed=list(.1,.3), iter=0)
aeq(u1, fit3a$u)
pen1 <- diag(c(rep(1/.1,9), rep(1/.3,9), 0,0))
aeq(imat1+pen1, as.matrix(igchol(fit3a$hmat)))

# Iteration 0, second form
idlist <- sort(outer(1:9, 0:1, paste, sep='/'))
mat1 <- matrix(diag(18), 18, dimnames=list(idlist, idlist))
mat2 <-  bdsBlock(idlist, rep(1:9, each=2))
mat3 <- diag(rep(0:1, 9))
dimnames(mat3) <- list(idlist, idlist)
fit3b <-  coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               varlist=coxmeMlist(list(mat2, mat3), rescale=F, pdcheck=F),
               vfixed=c(.1,.3), iter=0)

group <- strata(simdata$inst, simdata$trt, shortlabel=T, sep='/')
cox2 <- coxph(Surv(time, status) ~ group + age + trt, simdata,
              iter=0, x=T)
dt2 <- coxph.detail(cox2)
u2 <- apply(dt2$score, 2, sum)
aeq(fit3b$u, u2)

imat2 <- apply(dt2$imat, 1:2, sum)
pen2 <- matrix(0., 20,20)
pen2[1:18, 1:18] <- solve(.1*mat2 + .3*mat3)
aeq(imat2 + pen2, as.matrix(igchol(fit3b$hmat)))


# This function should map between coefficients
map <- matrix(0, 20,20)
for (i in 1:9) {
    map[i, 2*i -1] <- 1
    map[i+9, 2*i -1] <- -1
    map[i+9, 2*i] <- 1
    }
map[19,19] <- 1
map[20,20] <- 1
aeq(u1 %*%map, u2)     #scores map
aeq(pen2, t(map) %*% pen1 %*% map)  #penalty matrices map
aeq(imat2, t(map) %*% imat1 %*% map) # Cox variances map

# Now iteration 1
step1 <- solve(imat1 + pen1, u1)
step2 <- solve(imat2 + pen2, u2)
aeq(solve(fit3a$hmat, fit3a$u), step1)
aeq(solve(fit3b$hmat, fit3b$u), step2)


fit3a.1 <- coxme(Surv(time, status) ~ age + trt + (1|inst) + (trt|inst),
               simdata, vfixed=list(.1,.3), iter=1)
aeq(c(unlist(ranef(fit3a.1)), fixef(fit3a.1)), step1)

fit3b.1 <-  coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               varlist=coxmeMlist(list(mat2, mat3), rescale=F, pdcheck=F),
               vfixed=c(.1,.3), iter=1)
aeq(c(unlist(ranef(fit3b.1)), fixef(fit3b.1)), step2)
