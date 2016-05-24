###################################################
### chunk number 1: directories
###################################################
anadir <- "/home/komari/win/work/papers/bayesaft/CGDdata/"
dirsim1 <- paste(anadir, "anapaper1b/chain1", sep="")
dirsim2 <- paste(anadir, "anapaper1b/chain2", sep="")


###################################################
### chunk number 2: initop1
###################################################
library(bayesSurv)
data(cgd)
print(cgd[1:6,])


###################################################
### chunk number 3: initop2
###################################################
cgd$trtmt <- -(cgd$trtmt - 2)                                 
cgd$gender <- cgd$gender - 1                                  
cgd$inherit <- cgd$inherit - 1                                
cgd$cortico <- -(cgd$cortico - 2)                             
cgd$prophy <- -(cgd$prophy - 2)                               
cgd$gender <- factor(cgd$gender, labels = c("male", "female"))
cgd$inherit <- factor(cgd$inherit, labels = c("X-l", "AuRec"))
cgd$hcat <- factor(cgd$hcat, labels = c("US-NIH", "US-other", "EU-Am", "EU-other"))

cgd$event <- -(cgd$event - 2) 


###################################################
### chunk number 4: initop3
###################################################
cgd$time <- cgd$T1 - cgd$T2
npatient <- length(unique(cgd$ID))
nobs <- dim(cgd)[1]
print(cgd[1:6,])


###################################################
### chunk number 5: fitinit1
###################################################
ifit <- survreg(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat + frailty(ID, dist = "gaussian"),
                dist = "lognormal", data = cgd)
resid <- ifit$y[, 1] - ifit$linear.predictors 
R <- max(resid) - min(resid)
ifit2 <- survreg(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat,
                 dist = "lognormal", data = cgd)


###################################################
### chunk number 6: fitinit2
###################################################
summary(ifit)
print(R)


###################################################
### chunk number 7: fitinit3
###################################################
summary(ifit2)


###################################################
### chunk number 8: prior1
###################################################
X <- bayessurvreg1(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat + cluster(ID),
                   random = ~1, data = cgd, onlyX = TRUE)
nregres <- dim(X)[2]
X[1:3,]


###################################################
### chunk number 9: mixprior1
###################################################
prior <- list()


###################################################
### chunk number 10: mixprior2
###################################################
prior$kmax <- 30
prior$k.prior <- "poisson"
prior$poisson.k <- 5


###################################################
### chunk number 11: mixprior3
###################################################
prior$dirichlet.w <- 1


###################################################
### chunk number 12: mixprior4
###################################################
prior$mean.mu <- 3.66
prior$var.mu <- 5^2


###################################################
### chunk number 13: mixprior5
###################################################
prior$shape.invsig2 <- 2.0 
prior$shape.hyper.invsig2 <- 0.2
prior$rate.hyper.invsig2 <- 0.1


###################################################
### chunk number 14: mixprior6
###################################################
prior$pi.split <- c(1, rep(0.5, prior$kmax - 2), 0)


###################################################
### chunk number 15: mixprior7
###################################################
prior$pi.birth <- c(1, rep(0.5, prior$kmax - 2), 0)


###################################################
### chunk number 16: mixprior8
###################################################
prior$Eb0.depend.mix <- FALSE


###################################################
### chunk number 17: mixprior9
###################################################
print(prior)


###################################################
### chunk number 18: betaprior1
###################################################
prior.beta.gibbs <- list()
prior.beta.gibbs$mean.prior <- rep(0, nregres)
prior.beta.gibbs$var.prior <- rep(1e3, nregres)


###################################################
### chunk number 19: betaprior2
###################################################
print(prior.beta.gibbs)


###################################################
### chunk number 20: betaprior3
###################################################
prior.beta.mh1 <- list()
prior.beta.mh1$mean.prior <- rep(0, nregres)
prior.beta.mh1$var.prior <- rep(1e3, nregres)


###################################################
### chunk number 21: betaprior4
###################################################
prior.beta.mh1$blocks <- list()
prior.beta.mh1$blocks$ind.block <- list()


###################################################
### chunk number 22: betaprior5
###################################################
prior.beta.mh1$blocks$ind.block[[1]] <- 1:9            
nblock <- length(prior.beta.mh1$blocks$ind.block)


###################################################
### chunk number 23: betaprior6
###################################################
vars <- c(0.15, 0.20, 0.0003, 1.3, 0.08, 0.25, 0.1, 0.35, 0.35)
cors <- c(1,  0.1,  0,     0.1,  0.15, 0,    0.4, 0.1, 0.2,
          1, -0.15, 0.15, -0.2, -0.3,  0.2, -0.1, 0,
          1, -0.2,  0.15,  0.2,  0.3,  0.2,  0.1,
          1,  0.2, -0.5,   0.2, -0.4,  0.4,
          1,  0.15, 0.5,   0.3,  0.4,
          1,  0.15, 0.15,  0,
          1,  0.35, 0.65,
          1,  0.2,
          1)
corsm <- diag(9)
corsm[lower.tri(corsm, diag = TRUE)] <- cors
corsm[upper.tri(corsm, diag = FALSE)] <- t(corsm)[upper.tri(t(corsm), diag = FALSE)]
covm <- diag(sqrt(vars)) %*% corsm %*% diag(sqrt(vars))


###################################################
### chunk number 24: betaprior7
###################################################
print(corsm)


###################################################
### chunk number 25: betaprior8
###################################################
print(round(covm, digits=3))


###################################################
### chunk number 26: betaprior9
###################################################
prior.beta.mh1$blocks$cov.prop <- list()
prior.beta.mh1$blocks$cov.prop[[1]] <- covm[lower.tri(covm, diag = TRUE)]


###################################################
### chunk number 27: betaprior10
###################################################
prior.beta.mh1$type.upd <- rep("random.walk.metropolis", nblock)


###################################################
### chunk number 28: betaprior11
###################################################
prior.beta.mh1$weight.unif <- rep(0.05, nblock)   


###################################################
### chunk number 29: betaprior12
###################################################
prior.beta.mh1$half.range.unif <- c(0.25, 0.25, 0.01, 1.0, 0.15, 0.25, 0.3, 1.0, 1.0) 


###################################################
### chunk number 30: betaprior13
###################################################
print(prior.beta.mh1)


###################################################
### chunk number 31: betaprior14
###################################################
prior.beta.mh2 <- list()
prior.beta.mh2$mean.prior <- rep(0, nregres)
prior.beta.mh2$var.prior <- rep(1e3, nregres)


###################################################
### chunk number 32: betaprior15
###################################################
prior.beta.mh2$blocks <- list()
prior.beta.mh2$blocks$ind.block <- list()
prior.beta.mh2$blocks$ind.block[[1]] <- 1:6 
prior.beta.mh2$blocks$ind.block[[2]] <- 7:9 
nblock <- length(prior.beta.mh2$blocks$ind.block)


###################################################
### chunk number 33: betaprior16
###################################################
vars <- c(0.1, 0.35, 0.35)
cors <- c(1, 0.9, 0.9,
          1, 0.9,
          1)
corsm <- diag(3)
corsm[lower.tri(corsm, diag = TRUE)] <- cors
corsm[upper.tri(corsm, diag = FALSE)] <- t(corsm)[upper.tri(t(corsm), diag = FALSE)]
covm <- diag(sqrt(vars)) %*% corsm %*% diag(sqrt(vars))


###################################################
### chunk number 34: betaprior17
###################################################
print(corsm)


###################################################
### chunk number 35: betaprior18
###################################################
print(covm)


###################################################
### chunk number 36: betaprior19
###################################################
prior.beta.mh2$blocks$cov.prop <- list()
prior.beta.mh2$blocks$cov.prop[[1]] <- NULL
prior.beta.mh2$blocks$cov.prop[[2]] <- covm[lower.tri(covm, diag = TRUE)]


###################################################
### chunk number 37: betaprior20
###################################################
prior.beta.mh2$type.upd <- c("gibbs", "random.walk.metropolis")


###################################################
### chunk number 38: betaprior21
###################################################
prior.beta.mh2$weight.unif <- c(0.05, 0.05)   


###################################################
### chunk number 39: betaprior22
###################################################
prior.beta.mh2$half.range.unif <- c(0.25, 0.25, 0.01, 1.0, 0.15, 0.25, 0.3, 1.0, 1.0) 


###################################################
### chunk number 40: betaprior23
###################################################
print(prior.beta.mh2)


###################################################
### chunk number 41: bprior1
###################################################
prior.b.gamma <- list()


###################################################
### chunk number 42: bprior2
###################################################
prior.b.gamma$prior.D <- "inv.wishart"
prior.b.gamma$df.D <- 0.002                    # tau = degrees of freedom of an inverse-Wishart distribution
prior.b.gamma$scale.D <- 0.002                 # sb = scale of an inverse-Wishart distribution


###################################################
### chunk number 43: bprior3
###################################################
prior.b.gamma$type.upd <- "gibbs"


###################################################
### chunk number 44: bprior4
###################################################
print(prior.b.gamma)


###################################################
### chunk number 45: bprior5
###################################################
prior.b.unif <- list()
prior.b.unif$prior.D <- "sduniform"


###################################################
### chunk number 46: bprior6
###################################################
prior.b.unif$scale.D <- 100


###################################################
### chunk number 47: bprior6b
###################################################
prior.b.unif$type.upd <- "gibbs"


###################################################
### chunk number 48: bprior7
###################################################
print(prior.b.unif)


###################################################
### chunk number 49: revjump1
###################################################
prop.revjump <- list()


###################################################
### chunk number 50: revjump2
###################################################
prop.revjump$algorithm <- "correlated.av"


###################################################
### chunk number 51: revjump3
###################################################
prop.revjump$moody.ring <- c(0.1, 0.05)


###################################################
### chunk number 52: revjump4
###################################################
prop.revjump$transform.split.combine <- "brooks"
prop.revjump$transform.split.combine.parms <- c(2, 2, 2, 2, 1, 1)


###################################################
### chunk number 53: revjump5
###################################################
prop.revjump$transform.birth.death <- "richardson.green"


###################################################
### chunk number 54: revjump6
###################################################
print(prop.revjump)


###################################################
### chunk number 55: init11
###################################################
init1 <- list()


###################################################
### chunk number 56: init12
###################################################
init1$iter <- 0


###################################################
### chunk number 57: init13
###################################################
init1$mixture <- c(1,
                   1, rep(0, prior$kmax - 1),
                   3.9, rep(0, prior$kmax - 1),
                   1.2, rep(0, prior$kmax - 1))


###################################################
### chunk number 58: init14
###################################################
init1$beta <- c(1.1, -0.66, 0.04, -1.76, 0.94, 1.03, 0.37, 1.22, 0.82)


###################################################
### chunk number 59: init15
###################################################
init1$D <- 0.16


###################################################
### chunk number 60: init16
###################################################
init1$b <- rep(0, npatient)    


###################################################
### chunk number 61: init17
###################################################
init1$y <- NULL


###################################################
### chunk number 62: init18
###################################################
init1$r <- rep(1, nobs)


###################################################
### chunk number 63: init19
###################################################
init1$otherp <- rgamma(1, shape = prior$shape.hyper.invsig2, rate = prior$rate.hyper.invsig2)


###################################################
### chunk number 64: init110
###################################################
init1$u <- c(runif(1), 0, 0, runif(3*(prior$kmax - 1)))


###################################################
### chunk number 65: init111
###################################################
print(init1)


###################################################
### chunk number 66: init21
###################################################
init2 <- list()


###################################################
### chunk number 67: init22
###################################################
init2$iter <- 0


###################################################
### chunk number 68: init23
###################################################
init2$mixture <- c(2,
                   0.5, 0.5, rep(0, prior$kmax - 2),
                   2.5, 5.5, rep(0, prior$kmax - 2),
                   1, 1, rep(0, prior$kmax - 2))


###################################################
### chunk number 69: init24
###################################################
init2$beta <- rep(0, nregres)


###################################################
### chunk number 70: init25
###################################################
init2$D <- 0.05


###################################################
### chunk number 71: init26
###################################################
init2$b <- rnorm(npatient, 0, sqrt(init2$D))


###################################################
### chunk number 72: init27
###################################################
init2$y <- NULL


###################################################
### chunk number 73: init28
###################################################
init2$r <- c(rep(1, 102), rep(2, 101))


###################################################
### chunk number 74: init29
###################################################
init2$otherp <- rgamma(1, shape = prior$shape.hyper.invsig2, rate = prior$rate.hyper.invsig2)


###################################################
### chunk number 75: init210
###################################################
init2$u <- c(runif(1), 0, 0, runif(3*(prior$kmax - 1)))


###################################################
### chunk number 76: forrun1
###################################################
store <- list(y = FALSE, r = FALSE, u = FALSE, b = TRUE, MHb = FALSE, regresres = FALSE)


###################################################
### chunk number 77: forrun2
###################################################
nsimul <- list(niter = 1000, nthin = 3, nburn = 500, nnoadapt = 0, nwrite = 500)


###################################################
### chunk number 78: forrun3
###################################################
dir.create("cgdchain1test")
dir.create("cgdchain2test")
dirsim1test <- paste(getwd(), "/cgdchain1test", sep = "")
dirsim2test <- paste(getwd(), "/cgdchain2test", sep = "")


###################################################
### chunk number 79: runsimChain1
###################################################
simul1 <- bayessurvreg1(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat + cluster(ID),
                        random = ~1,
                        data = cgd, dir = dirsim1test, nsimul = nsimul,
                        prior = prior, prior.beta = prior.beta.gibbs, prior.b = prior.b.unif, prop.revjump = prop.revjump,
                        init = init1, store = store)


###################################################
### chunk number 80: runsimChain2
###################################################
simul2 <- bayessurvreg1(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat + cluster(ID),
                        random = ~1,
                        data = cgd, dir = dirsim2test, nsimul = nsimul,
                        prior = prior, prior.beta = prior.beta.gibbs, prior.b = prior.b.unif, prop.revjump = prop.revjump,
                        init = init2, store = store)


###################################################
### chunk number 81: predict1
###################################################
nnewpat <- 8
nID <- 1:nnewpat
ntrtmt <-          c(0, 1, 0, 1, 0, 1, 0, 1)
ninherit <- factor(c(0, 0, 1, 1, 0, 0, 1, 1), levels = 0:1, labels = c("X-l", "AuRec"))
nage <- rep(14.6, nnewpat)
ncortico <- rep(0, nnewpat)
nprophy <- rep(1, nnewpat)
ngender <-  factor(c(0, 0, 0, 0, 1, 1, 1, 1), levels = 0:1, labels = c("male", "female"))
nhcat <- factor(rep(2, nnewpat), levels = 1:4, labels = c("US-NIH", "US-other", "EU-Am", "EU-other"))

ntime <- rep(1, nnewpat)
nevent <- rep(0, nnewpat)


###################################################
### chunk number 82: predict2
###################################################
preddata <- data.frame(ID = nID, trtmt = ntrtmt, inherit = ninherit, age = nage, cortico = ncortico, prophy = nprophy,
                       gender = ngender, hcat = nhcat, time = ntime, event = nevent)
print(preddata)


###################################################
### chunk number 83: predict3
###################################################
predict <- list(Et = TRUE, t = FALSE, Surv = TRUE, hazard = TRUE, cum.hazard = FALSE)
store <- list(Et = FALSE, t = FALSE, Surv = FALSE, hazard = FALSE, cum.hazard = FALSE)


###################################################
### chunk number 84: predict4
###################################################
grid <- seq(1, 401, by = 2.5)


###################################################
### chunk number 85: runPredict
###################################################
simulp <- predictive(Surv(time, event) ~ trtmt + inherit + age + cortico + prophy + gender + hcat + cluster(ID),
                     random = ~1, data = preddata, dir = dirsim1test,
                     quantile = c(0, 0.025, 0.5, 0.975, 1),
                     skip = 0, by = 1, predict = predict, store = store, grid = grid,
                     Eb0.depend.mix = FALSE, type = "mixture")


###################################################
### chunk number 86: fig1
###################################################
labels <- c("Male, X-l, plcb", "Male, X-l, trt", "Male, AR, plcb", "Male, AR, trt",
            "Female, X-l, plcb", "Female, X-l, trt", "Female, AR, plcb", "Female, AR, trt")
par(mfrow = c(4, 2))
for(i in 1:8){
  gridS <- scan(paste(dirsim1, "/quantS", i, ".sim", sep = ""), nlines = 1) 
  Sfun <- read.table(paste(dirsim1, "/quantS", i, ".sim", sep = ""), header = TRUE)
  rownames(Sfun) <- c("0%", "2.5%", "50%", "97.5%", "100%", "mean")
  plot(gridS, Sfun["mean", ], type = "l", lty = 1, ylim = c(0, 1), xlab = "Time (days)", ylab = "Survivor", bty = "n")  
  lines(gridS, Sfun["2.5%", ], lty = 2)
  lines(gridS, Sfun["97.5%", ], lty = 2)
  title(main = labels[i])
}  


###################################################
### chunk number 87: fig2
###################################################
labels <- c("Male, X-l, plcb", "Male, X-l, trt", "Male, AR, plcb", "Male, AR, trt",
            "Female, X-l, plcb", "Female, X-l, trt", "Female, AR, plcb", "Female, AR, trt")
par(mfrow = c(4, 2))
for(i in 1:8){
  gridhaz <- scan(paste(dirsim1, "/quanthazard", i, ".sim", sep = ""), nlines = 1) 
  hfun <- read.table(paste(dirsim1, "/quanthazard", i, ".sim", sep = ""), header = TRUE)
  rownames(hfun) <- c("0%", "2.5%", "50%", "97.5%", "100%", "mean")
  plot(gridhaz, hfun["97.5%", ], type = "l", lty = 2, xlab = "Time (days)", ylab = "Hazard", bty = "n")  
  lines(gridhaz, hfun["mean", ], lty = 1)
  lines(gridhaz, hfun["2.5%", ], lty = 2)
  title(main = labels[i])
}  


###################################################
### chunk number 88: fig3
###################################################
gg <- c("Male", "Female")
par(mfrow = c(2, 1))
for (j in 1:2){
  for(i in 1:4){
    gridS <- scan(paste(dirsim1, "/quantS", (j-1)*4+i, ".sim", sep = ""), nlines = 1) 
    Sfun <- read.table(paste(dirsim1, "/quantS", (j-1)*4+i, ".sim", sep = ""), header = TRUE)
    rownames(Sfun) <- c("0%", "2.5%", "50%", "97.5%", "100%", "mean")
    if (i == 1) plot(gridS, Sfun["mean", ], type = "l", lty = 1, ylim = c(0, 1), xlab = "Time (days)", ylab = "Survivor", bty = "n")
    else        lines(gridS, Sfun["mean", ], lty = i)
    legend(0, 0.6, legend = c("trtmt, X-l", "trtmt, AuRec", "placebo, X-l", "placebo, AuRec"), lty = c(2, 4, 1, 3), bty = "n")
  }
  title(main = gg[j])
}  


###################################################
### chunk number 89: fig4
###################################################
gg <- c("Male", "Female")
par(mfrow = c(2, 1))
leg <- c(0.0080, 0.0040)
ylim <- c(0.0100, 0.0040)
for (j in 1:2){
  for(i in 1:4){
    gridhaz <- scan(paste(dirsim1, "/quanthazard", (j-1)*4+i, ".sim", sep = ""), nlines = 1) 
    hfun <- read.table(paste(dirsim1, "/quanthazard", (j-1)*4+i, ".sim", sep = ""), header = TRUE)
    rownames(hfun) <- c("0%", "2.5%", "50%", "97.5%", "100%", "mean")
    if (i == 1) plot(gridhaz, hfun["mean", ], type = "l", lty = 1, ylim = c(0, ylim[j]), xlab = "Time (days)", ylab = "Hazard", bty = "n")
    else        lines(gridhaz, hfun["mean", ], lty = i)
    legend(200, leg[j], legend = c("trtmt, X-l", "trtmt, AuRec", "placebo, X-l", "placebo, AuRec"), lty = c(2, 4, 1, 3), bty = "n")
  }
  title(main = gg[j])
}  


###################################################
### chunk number 90: predError1
###################################################
dgrid <- seq(-2, 9, length = 100)
dgrids <- seq(-3, 3, length = 100)
dens <- list()
for (ch in 1:2){
  dens[[ch]] <- bayesDensity(dir = get(paste("dirsim", ch, sep = "")), grid = dgrid, stgrid = dgrids)
}  


###################################################
### chunk number 91: fig5
###################################################
par(bty = "n", mfrow = c(2, 2))
for (ch in 1:2){
  xlim <- c(-2, 8); xleg <- -2; yleg <- 0.3; ylim <- c(0, 0.3)
  plot(dens[[ch]], k.cond = 0:9, standard = FALSE, dim.plot = FALSE, xlim = xlim, ylim = ylim, xleg = xleg, yleg = yleg, main = "")
  title(main = paste("Unstandardized, Chain ", ch, sep = ""))
}
for (ch in 1:2){
  xlim <- c(-2.5, 2.5); xleg <- -2.5; yleg <- 0.7; ylim <- c(0, 0.7)
  plot(dens[[ch]], k.cond = 0:9, standard = TRUE, dim.plot = FALSE, xlim = xlim, ylim = ylim, xleg = xleg, yleg = yleg, main = "")
  title(main = paste("Standardized, Chain ", ch, sep = ""))
}


###################################################
### chunk number 92: predb1
###################################################
ids <- unique(cgd$ID)
bb <- list()
for (ch in 1:2){
  bb[[ch]] <- matrix(scan(paste(get(paste("dirsim", ch, sep = "")), "/b.sim", sep = ""), skip = 1), ncol = 128, byrow = TRUE)
  colnames(bb[[ch]]) <- ids
} 
bbs <- rbind(bb[[1]], bb[[2]])


###################################################
### chunk number 93: predb2
###################################################
b.mean <- apply(bbs, 2, mean)
b.median <- apply(bbs, 2, quantile, 0.50)
b.low <- apply(bbs, 2, quantile, 0.025)
b.up <- apply(bbs, 2, quantile, 0.975)


###################################################
### chunk number 94: predb3
###################################################
n <- dim(cgd)[1]
id1 <- cgd$ID[1:(n-1)]; id2 <- cgd$ID[2:n]; difid <- c(1, id2 - id1)
first <- difid > 0           
frevent <- table(cgd$ID)
freqv <- as.numeric(frevent)
frval <- data.frame(ID = cgd$ID[first], trtmt = cgd$trtmt[first], freq = as.numeric(frevent),
                    b.mean, b.median, b.low, b.up, nevent = freqv)
frval <- frval[order(frval$trtmt), ]
frval <- frval[order(frval$freq), ]


###################################################
### chunk number 95: fig15
###################################################
par(bty="n")
plot(frval$b.mean, 1:128, type = "p", pch = 20, xlim = c(-3, 2.5), bty = "n", ylab = "Patient", xlab = "Random effect")
lines(frval$b.mean, 1:128)
points(frval$b.low, 1:128, pch = 15, cex = 0.5)
points(frval$b.up, 1:128, pch = 15, cex = 0.5)
for (pat in 1:n){
  lines(c(frval$b.low[pat], frval$b.up[pat]), c(pat, pat), lty = 1)
}
title(main = "Individual random effects")
abline(h = 84.5, lty = 2)
abline(h = 112.5, lty = 2)
abline(h = 120.5, lty = 2)
abline(v = 0, lty = 1)
text(-3, 40, "1 event", pos = 4)
text(-3, 30, "/patient", pos = 4)
text(-3, 100,  "2 events", pos = 4)
text(-3, 90, "/patient", pos = 4)
text(-3, 115, "3", pos = 4)
text(-3, 127, ">= 4", pos = 4)  


###################################################
### chunk number 96: summ1
###################################################
library(coda)
nchains <- 2


###################################################
### chunk number 97: summ2
###################################################
sdb <- list()
logscale <- list()
for (ch in 1:nchains){
  sdb[[ch]] <- read.table(paste(get(paste("dirsim", ch, sep = "")), "/D.sim", sep = ""), header = TRUE)
  sdb[[ch]] <- data.frame(sdb = sqrt(sdb[[ch]][,2]))
  logscale[[ch]] <- read.table(paste(get(paste("dirsim", ch, sep = "")), "/mixmoment.sim", sep = ""), header = TRUE)
  logscale[[ch]] <- data.frame(logscale = sqrt(logscale[[ch]][,2]))
}  


###################################################
### chunk number 98: summ3
###################################################
pars <- list()
for (ch in 1:nchains){
  pars[[ch]] <- files2coda(files = c("beta.sim", "mixmoment.sim"), data.frames = c("sdb", "logscale"),
                           thin = 1, dir = paste(get(paste("dirsim", ch, sep = ""))), chain = ch)
}


###################################################
### chunk number 99: summ4
###################################################
parsls <- mcmc.list(pars[[1]], pars[[2]])
rm(list = c("pars", "sdb", "logscale"))


###################################################
### chunk number 100: summ5
###################################################
dimnames(parsls[[1]])[[2]]


###################################################
### chunk number 101: summ6
###################################################
quant <- c(0, 0.025, 0.5, 0.75, 0.975, 1)
means <- list();  quantiles <- list();  summ <- list()
for (ch in 1:nchains){
  means[[ch]] <- apply(parsls[[ch]], 2, mean)
  quantiles[[ch]] <- apply(parsls[[ch]], 2, quantile, quant)
  summ[[ch]] <- rbind(means[[ch]], quantiles[[ch]])
  rownames(summ[[ch]])[1] <- "mean"
}  
names(summ) <- paste("Chain ", 1:nchains, sep = "")
print(summ)


###################################################
### chunk number 102: fig6
###################################################
par(bty = "n")
par(mfrow = c(2, 2))
kall <- c(parsls[[1]][, "k"], parsls[[2]][, "k"])
hist(kall, xlab = "k", prob = TRUE, main = "Both chains", breaks = 0:22)
plot.new()
hist(parsls[[1]][, "k"], xlab = "k", prob = TRUE, main = "Chain 1", breaks = 0:22)
hist(parsls[[2]][, "k"], xlab = "k", prob = TRUE, main = "Chain 2", breaks = 0:22)


###################################################
### chunk number 103: fig7
###################################################
ch <- 1
par(bty = "n")
par(mfrow = c(3, 3))
for (i in 1:9){
   densplot(parsls[[ch]][,i], show.obs = FALSE, bty = "n")
   title(main = paste(attr(parsls[[ch]], "dimnames")[[2]][i], ", chain ", ch, sep = ""))
}


###################################################
### chunk number 104: fig8
###################################################
ch <- 1
par(mfrow = c(2, 2))
densplot(parsls[[ch]][, "k"], show.obs = FALSE, bty = "n")
title(main = paste("k, chain", ch, sep = ""))
densplot(parsls[[ch]][, "sdb"], show.obs = FALSE, bty = "n")
title(main = paste("Std. Dev. of b, chain", ch, sep = ""))
densplot(parsls[[ch]][, "Intercept"], show.obs = FALSE, bty = "n")
title(main = paste("Intercept, chain ", ch, sep = ""))
densplot(parsls[[ch]][, "Scale"], show.obs = FALSE, bty = "n")
title(main = paste("Error Scale, chain ", ch, sep = ""))


###################################################
### chunk number 105: fig9
###################################################
ch <- 1
par(bty = "n")
par(mfrow = c(3, 3))
autocorr.plot(parsls[[ch]][, 1:9], ask = FALSE, sub = paste("Chain ", ch, sep = ""))


###################################################
### chunk number 106: fig10
###################################################
par(bty = "n")
par(mfrow = c(2, 2))
plot.new()
autocorr.plot(parsls[[ch]][, 13], auto.layout = FALSE, ask = FALSE, sub = paste("Chain ", ch, sep = ""), main = "sdb")
autocorr.plot(parsls[[ch]][, 11:12], auto.layout = FALSE, ask = FALSE, sub = paste("Chain ", ch, sep = ""))


###################################################
### chunk number 107: crosscorr
###################################################
croscor <- lapply(parsls, crosscorr)
croscor <- lapply(croscor, round, digits = 2)
names(croscor) <- paste("Chain ", 1:nchains, sep = "")
print(croscor)


###################################################
### chunk number 108: GelmanRubin
###################################################
gelm <- gelman.diag(parsls)
rownames(gelm$psrf) <- dimnames(parsls[[1]])[[2]]
print(gelm)


###################################################
### chunk number 109: fig11
###################################################
ch <- 1
par(bty = "n")
par(mfrow = c(2, 2));  traceplot2(parsls[[ch]], chains = 1:3, sub = paste("Chain ", ch, sep = "")); 


###################################################
### chunk number 110: fig12
###################################################
par(bty = "n")
par(mfrow = c(2, 2));  traceplot2(parsls[[ch]], chains = 4:6, sub = paste("Chain ", ch, sep = "")); 


###################################################
### chunk number 111: fig13
###################################################
par(bty = "n")
par(mfrow = c(2, 2));  traceplot2(parsls[[ch]], chains = 7:9, sub = paste("Chain ", ch, sep = "")); 


###################################################
### chunk number 112: fig14
###################################################
par(bty = "n")
par(mfrow = c(2, 2))
traceplot2(parsls[[ch]], chains = 10, sub = paste("Chain ", ch, sep = ""))
traceplot2(parsls[[ch]], chains = 13, sub = paste("Chain ", ch, sep = ""))
traceplot2(parsls[[ch]], chains = 11:12, sub = paste("Chain ", ch, sep = ""))


###################################################
### chunk number 113: performRevJump
###################################################
mh <- list() 
averMH <- list()
for (ch in 1:nchains){
  mh[[ch]] <- files2coda(files = c("MHinfo.sim"), start = 1, thin = 1, dir = get(paste("dirsim", ch, sep = "")))
  averMH[[ch]] <- apply(mh[[ch]][, c(1, 3)], 2, mean)
}  
for (ch in 1:nchains){
  cat("Chain ", ch, ":\n", sep = "")
  print(averMH[[ch]])
}  


###################################################
### chunk number 114: cleaning
###################################################
files1 <- dir("./cgdchain1test")
files2 <- dir("./cgdchain2test")
file.remove(paste("./cgdchain1test/", files1, sep = ""))
file.remove(paste("./cgdchain2test/", files2, sep = ""))
file.remove("cgdchain1test")
file.remove("cgdchain2test")


