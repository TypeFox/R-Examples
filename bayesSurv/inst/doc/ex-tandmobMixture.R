###################################################
### chunk number 1: directories
###################################################
anadir <- "/home/komari/win/work/papers/bayesaft/Tandmob/tana1/"
docdir <- "/home/komari/win/work/papers/bayesaft/RforCRAN/tandmob/"
dirsim <- character()
dirsim[1] <- paste(anadir, "chain1small", sep="")
dirsim[2] <- paste(anadir, "chain2small", sep="")


###################################################
### chunk number 2: initop1
###################################################
library(bayesSurv)
data(tandmob2)


###################################################
### chunk number 3: initop2
###################################################
pteeth <- c(14, 24, 34, 44, 15, 25, 35, 45)
dteeth <- pteeth + 40
nteeth <- length(pteeth)


###################################################
### chunk number 4: initop3
###################################################
upper4 <- c(1, 1, 0, 0, 0, 0, 0, 0)
lower4 <- c(0, 0, 1, 1, 0, 0, 0, 0)
upper5 <- c(0, 0, 0, 0, 1, 1, 0, 0)
lower5 <- c(0, 0, 0, 0, 0, 0, 1, 1)


###################################################
### chunk number 5: initop4
###################################################
childvars <- c("IDNR", "GENDER", "GENDERNum", "DOB", "PROVINCE", "EDUC")
timevars <- paste(c("EBEG.", "EEND."), rep(pteeth, rep(2, nteeth)), sep="")
dmfvars <- paste("T", dteeth, ".DMF", sep="")
vars <- c(childvars, timevars, dmfvars)
print(vars)


###################################################
### chunk number 6: initop5
###################################################
tandmob2 <- tandmob2[, vars]


###################################################
### chunk number 7: init5b
###################################################
names(tandmob2)[names(tandmob2) == "GENDERNum"] <- "GIRL"
childvars[childvars == "GENDERNum"] <- "GIRL"
vars[vars == "GENDERNum"] <- "GIRL"


###################################################
### chunk number 8: initop6
###################################################
is.dmf <- !is.na(tandmob2$T54.DMF) & !is.na(tandmob2$T64.DMF) & !is.na(tandmob2$T74.DMF) & !is.na(tandmob2$T84.DMF) &
          !is.na(tandmob2$T55.DMF) & !is.na(tandmob2$T65.DMF) & !is.na(tandmob2$T75.DMF) & !is.na(tandmob2$T85.DMF) 
tandmob2 <- tandmob2[is.dmf,]


###################################################
### chunk number 9: initop7
###################################################
startage <- 5.0
tandmob2[, timevars] <- tandmob2[, timevars] - startage


###################################################
### chunk number 10: initop8
###################################################
nFromGender <- 50
ids <- read.table(paste(docdir, "IDsampled", nFromGender, ".dat", sep=""), header=TRUE)[,1]
wanna <- match(tandmob2$IDNR, ids, nomatch=0)
wanna <- as.logical(wanna)
wanna[wanna > 1] <- TRUE
tandmob2 <- tandmob2[wanna,]


###################################################
### chunk number 11: initop5a
###################################################
print(tandmob2[1:10, ])


###################################################
### chunk number 12: initop9
###################################################
nchild <- dim(tandmob2)[1]
longdata <- list()
for (i in 1:length(childvars)){
  longdata[[i]] <- rep(tandmob2[[childvars[i]]], nteeth)
}
names(longdata) <- childvars
longdata <- as.data.frame(longdata)
longdata$TOOTH <- as.factor(rep(pteeth, rep(nchild, nteeth)))
longdata$UPPER4 <- rep(upper4, rep(nchild, nteeth))
longdata$LOWER4 <- rep(lower4, rep(nchild, nteeth))
longdata$UPPER5 <- rep(upper5, rep(nchild, nteeth))
longdata$LOWER5 <- rep(lower5, rep(nchild, nteeth))
dmf <- numeric(); for (i in 1:nteeth) dmf <- c(dmf, tandmob2[[dmfvars[i]]])
ebeg <- numeric(); for (i in 1:nteeth) ebeg <- c(ebeg, tandmob2[[paste("EBEG.", pteeth[i], sep="")]])
eend <- numeric(); for (i in 1:nteeth) eend <- c(eend, tandmob2[[paste("EEND.", pteeth[i], sep="")]])
longdata$DMF <- dmf
longdata$EBEG <- ebeg
longdata$EEND <- eend
rm(list=c("dmf", "ebeg", "eend"))
longdata <- longdata[order(longdata$IDNR), ]


###################################################
### chunk number 13: initop10
###################################################
print(longdata[1:16,])


###################################################
### chunk number 14: initval1
###################################################
ifit <- survreg(Surv(EBEG, EEND, type="interval2") ~ GIRL*DMF + (DMF+GIRL)*(LOWER4 + UPPER5 + LOWER5),
                 dist = "lognormal", data = longdata)
summary(ifit)


###################################################
### chunk number 15: prior1
###################################################
X <- bayessurvreg1(Surv(EBEG, EEND, type="interval2") ~ GIRL*DMF + (DMF+GIRL)*(LOWER4 + UPPER5 + LOWER5) + cluster(IDNR),
                   random =~ LOWER4 + UPPER5 + LOWER5, data = longdata[1:80,], onlyX = TRUE)


###################################################
### chunk number 16: prior2
###################################################
print(X[1:5,])
nregres <- dim(X)[2]
nobs <- dim(longdata)[1]


###################################################
### chunk number 17: priormix1
###################################################
prior <- list()


###################################################
### chunk number 18: priormix2
###################################################
prior$kmax <- 30
prior$k.prior <- "poisson"
prior$poisson.k <- 5


###################################################
### chunk number 19: priormix3
###################################################
prior$dirichlet.w <- 1


###################################################
### chunk number 20: priormix4
###################################################
prior$mean.mu <- 1.8
prior$var.mu <- 0.75^2


###################################################
### chunk number 21: priormix5
###################################################
prior$shape.invsig2 <- 2.0
prior$shape.hyper.invsig2 <- 0.2
prior$rate.hyper.invsig2 <- 0.1


###################################################
### chunk number 22: priormix6
###################################################
prior$pi.split <- c(1, rep(0.5, prior$kmax - 2), 0)


###################################################
### chunk number 23: priormix7
###################################################
prior$pi.birth <- c(1, rep(0.5, prior$kmax - 2), 0)


###################################################
### chunk number 24: priormix8
###################################################
prior$Eb0.depend.mix <- FALSE


###################################################
### chunk number 25: priormix9
###################################################
print(prior)


###################################################
### chunk number 26: priorbeta1
###################################################
prior.beta <- list()


###################################################
### chunk number 27: priorbeta2
###################################################
prior.beta$mean.prior <- rep(0, nregres)
prior.beta$var.prior <- rep(1e2, nregres)


###################################################
### chunk number 28: priorbeta3
###################################################
print(prior.beta)


###################################################
### chunk number 29: priorb1
###################################################
prior.b <- list()


###################################################
### chunk number 30: priorb2
###################################################
prior.b$df.D <- 4
prior.b$scale.D <- 0.002*c(1, 0, 0, 0, 1, 0, 0, 1, 0, 1)
prior.b$prior.D <- "inv.wishart"


###################################################
### chunk number 31: priorb3
###################################################
prior.b$type.upd <- "gibbs"


###################################################
### chunk number 32: priorb4
###################################################
print(prior.b)


###################################################
### chunk number 33: revjump1
###################################################
prop.revjump <- list()


###################################################
### chunk number 34: revjump2
###################################################
prop.revjump$algorithm <- "correlated.av"


###################################################
### chunk number 35: revjump3
###################################################
prop.revjump$moody.ring <- c(0.1, 0.05)


###################################################
### chunk number 36: revjump4
###################################################
prop.revjump$transform.split.combine <- "brooks"
prop.revjump$transform.split.combine.parms <- c(2, 2, 2, 2, 1, 1)


###################################################
### chunk number 37: revjump5
###################################################
prop.revjump$transform.birth.death <- "richardson.green"


###################################################
### chunk number 38: revjump6
###################################################
print(prop.revjump)


###################################################
### chunk number 39: init11
###################################################
init1 <- list()


###################################################
### chunk number 40: init12
###################################################
init1$iter <- 0


###################################################
### chunk number 41: init13
###################################################
init1$mixture <- c(1,
                   1, rep(0, prior$kmax - 1),
                   1.8, rep(0, prior$kmax - 1),
                   0.25^2, rep(0, prior$kmax - 1))


###################################################
### chunk number 42: init14
###################################################
init1$beta <- c(-0.09, -0.11, -0.01, 0.16, 0.17, 0.05, 0.02, 0.01, 0.03, -0.02, 0.01, 0)  


###################################################
### chunk number 43: init15
###################################################
init1$D <- c(1, 0, 0, 0, 1, 0, 0, 1, 0, 1)


###################################################
### chunk number 44: init16
###################################################
b0 <- rep(0, nchild)
b1 <- rnorm(nchild, 0, 0.3)
b2 <- rnorm(nchild, 0.16, 0.3)
b3 <- rnorm(nchild, 0.16, 0.3)
init1$b <- as.numeric(rbind(b0, b1, b2, b3))


###################################################
### chunk number 45: init17
###################################################
init1$y <- NULL


###################################################
### chunk number 46: init18
###################################################
init1$r <- rep(1, nobs)


###################################################
### chunk number 47: init19
###################################################
init1$otherp <- rgamma(1, shape = prior$shape.hyper.invsig2, rate = prior$rate.hyper.invsig2)


###################################################
### chunk number 48: init110
###################################################
init1$u <- c(runif(1), 0, 0, runif(3*(prior$kmax - 1)))


###################################################
### chunk number 49: init21
###################################################
init2 <- list()


###################################################
### chunk number 50: init22
###################################################
init2$iter <- 0


###################################################
### chunk number 51: init23
###################################################
init2$mixture <- c(1, 
                   1, rep(0, prior$kmax - 1),
                   1.7, rep(0, prior$kmax - 1),
                   0.20^2, rep(0, prior$kmax - 1))


###################################################
### chunk number 52: init24
###################################################
init2$beta <- rep(0, nregres)


###################################################
### chunk number 53: init25
###################################################
varb <- c(0.01, 0.01, 0.01, 0.01)
corb <- c(0.2, 0.2, 0.2,  0.2, 0.2,  0.2)
Corb <- diag(4)
Corb[lower.tri(Corb, diag=FALSE)] <- corb
Corb[upper.tri(Corb, diag=FALSE)] <- t(Corb)[upper.tri(t(Corb), diag=FALSE)]
Covb <- diag(sqrt(varb)) %*% Corb %*% diag(sqrt(varb))
init2$D <- Covb[lower.tri(Covb, diag=TRUE)]


###################################################
### chunk number 54: init26
###################################################
b0 <- rnorm(nchild, 0, sqrt(init2$D[1])) 
b1 <- rnorm(nchild, 0, sqrt(init2$D[5]))
b2 <- rnorm(nchild, 0.16, sqrt(init2$D[8])) 
b3 <- rnorm(nchild, 0.16, sqrt(init2$D[10]))
init2$b <- as.numeric(rbind(b0, b1, b2, b3))


###################################################
### chunk number 55: init27
###################################################
init2$y <- NULL


###################################################
### chunk number 56: init28
###################################################
init2$r <- rep(1, nobs)


###################################################
### chunk number 57: init29
###################################################
init2$otherp <- rgamma(1, shape = prior$shape.hyper.invsig2, rate = prior$rate.hyper.invsig2)


###################################################
### chunk number 58: init210
###################################################
init2$u <- c(runif(1), 0, 0, runif(3*(prior$kmax - 1)))


###################################################
### chunk number 59: run1
###################################################
store <- list(y = FALSE, r = FALSE, u = FALSE, b = FALSE, MHb = FALSE, regresres = FALSE)


###################################################
### chunk number 60: run2
###################################################
nsimul <- list(niter = 1000, nthin = 3, nburn = 500, nwrite = 500)


###################################################
### chunk number 61: run3
###################################################
dir.create("tandchain1test")
dir.create("tandchain2test")
dirsimtest <- character()
dirsimtest[1] <- paste(getwd(), "/tandchain1test", sep = "")
dirsimtest[2] <- paste(getwd(), "/tandchain2test", sep = "")


###################################################
### chunk number 62: runChain1
###################################################
sim1 <- bayessurvreg1(Surv(EBEG, EEND, type="interval2") ~ GIRL*DMF + (DMF+GIRL)*(LOWER4 + UPPER5 + LOWER5) + cluster(IDNR),
                random=~ LOWER4 + UPPER5 + LOWER5, data=longdata, dir=dirsimtest[1], nsimul=nsimul,
                prior=prior, init=init1, prop.revjump=prop.revjump, prior.beta=prior.beta, prior.b=prior.b,
                store=store)


###################################################
### chunk number 63: runChain2
###################################################
sim2 <- bayessurvreg1(Surv(EBEG, EEND, type="interval2") ~ GIRL*DMF + (DMF+GIRL)*(LOWER4 + UPPER5 + LOWER5) + cluster(IDNR),
                random=~ LOWER4 + UPPER5 + LOWER5, data=longdata, dir=dirsimtest[2], nsimul=nsimul,
                prior=prior, init=init2, prop.revjump=prop.revjump, prior.beta=prior.beta, prior.b=prior.b,
                store=store)


###################################################
### chunk number 64: pred1
###################################################
newIDNR    = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
newGIRL    = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
newUPPER4  = c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0)
newLOWER4  = c(0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0)
newUPPER5  = c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0)
newLOWER5  = c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
newDMF     = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1)
newEBEG    = rep(1, 16)
newEEND    = rep(NA, 16)                      
preddata <- data.frame(IDNR    = newIDNR,
                       GIRL    = newGIRL,
                       UPPER4  = newUPPER4,
                       LOWER4  = newLOWER4,
                       UPPER5  = newUPPER5,
                       LOWER5  = newLOWER5,
                       DMF     = newDMF,
                       EBEG    = newEBEG,
                       EEND    = newEEND)              
print(preddata)


###################################################
### chunk number 65: pred2
###################################################
predict <- list(Et = TRUE, t = FALSE, Surv = TRUE, hazard = TRUE, cum.hazard = FALSE)
store <- list(Et = FALSE, t = FALSE, Surv = FALSE, hazard = FALSE, cum.hazard = FALSE)


###################################################
### chunk number 66: pred3
###################################################
grid <- seq(6.5, 13, by=0.1)


###################################################
### chunk number 67: runPredict
###################################################
simulp1 <- predictive(Surv(EBEG, EEND, type="interval2") ~ GIRL*DMF + (DMF+GIRL)*(LOWER4 + UPPER5 + LOWER5) + cluster(IDNR),
                     random=~ LOWER4 + UPPER5 + LOWER5, time0=startage,
                     data = preddata, dir = dirsimtest[1], skip = 0, by = 1, predict = predict, store = store, grid = grid,
                     type = "mixture")


###################################################
### chunk number 68: emerc1
###################################################
ch <- 1
line <- numeric(13); line[1] <- 1; line[5] <- 2; line[9] <- 1; line[13] <- 2
col <- character(13); col[1] <- "blue"; col[5] <- "blue"; col[9] <- "red"; col[13] <- "red"
title <- c("Maxilla 4", "Mandible 4", "Maxilla 5", "Mandible 5")
lwd <- 2
xlim <- c(7, 12)
par(mfrow = c(2, 2))
for (i in 1:4){
  for(j in c(1,5,9,13)){
    gridS <- scan(paste(dirsim[ch], "/quantS", (i-1)+j, ".sim", sep = ""), nlines = 1) 
    Sfun <- read.table(paste(dirsim[ch], "/quantS", (i-1)+j, ".sim", sep = ""), header = TRUE)
    rownames(Sfun) <- c("0%", "2.5%", "50%", "97.5%", "100%", "mean")
    if (j == 1) plot(gridS, 1-Sfun["mean", ], type="l", lty=line[1], ylim=c(0, 1), xlab="Age", ylab="Emergence", bty="n", lwd=lwd, col=col[1], xlim=xlim)
    else        lines(gridS, 1-Sfun["mean", ], lty=line[j], lwd=lwd, col=col[j])
    legend(7, 1, legend = c("Girl, dmf>0", "Girl, sound", "Boy, dmf>0", "Boy, sound"), lty=line[c(13,9,5,1)], col=col[c(13,9,5,1)], lwd=lwd, bty = "n")
  }
  title(main = title[i])
}  


###################################################
### chunk number 69: summ1
###################################################
nchains <- length(dirsim)
print(nchains)


###################################################
### chunk number 70: summ2
###################################################
covb <- list()
covt <- list()
for (ch in 1:nchains){
  dd <- read.table(paste(dirsim[ch], "/D.sim", sep=""), header=TRUE)

  varb <- cbind(dd$D.1.1, dd$D.2.2, dd$D.3.3, dd$D.4.4)
  sdb <- sqrt(varb)
  corb12 <- dd$D.2.1/sqrt(dd$D.1.1*dd$D.2.2)
  corb13 <- dd$D.3.1/sqrt(dd$D.1.1*dd$D.3.3)
  corb14 <- dd$D.4.1/sqrt(dd$D.1.1*dd$D.4.4)
  corb23 <- dd$D.3.2/sqrt(dd$D.2.2*dd$D.3.3)
  corb24 <- dd$D.4.2/sqrt(dd$D.2.2*dd$D.4.4)
  corb34 <- dd$D.4.3/sqrt(dd$D.3.3*dd$D.4.4)
  covb[[ch]] <- data.frame(sdb, corb12, corb13, corb14, corb23, corb24, corb34)
  colnames(covb[[ch]]) <- c(paste("sd", 1:4, sep=""), paste("cor", c(12, 13, 14, 23, 24, 34), sep=""))
  
  sdt1 <- sqrt(dd$D.1.1)
  sdt2 <- sqrt(dd$D.1.1 + dd$D.2.2 + 2*dd$D.2.1)
  sdt3 <- sqrt(dd$D.1.1 + dd$D.3.3 + 2*dd$D.3.1)
  sdt4 <- sqrt(dd$D.1.1 + dd$D.4.4 + 2*dd$D.4.1)
  cort12 <- (dd$D.1.1 + dd$D.2.1)/(sdt1*sdt2)
  cort13 <- (dd$D.1.1 + dd$D.3.1)/(sdt1*sdt3)
  cort14 <- (dd$D.1.1 + dd$D.4.1)/(sdt1*sdt4)
  cort23 <- (dd$D.1.1 + dd$D.2.1 + dd$D.3.1 + dd$D.3.2)/(sdt2*sdt3)
  cort24 <- (dd$D.1.1 + dd$D.2.1 + dd$D.4.1 + dd$D.4.2)/(sdt2*sdt4)
  cort34 <- (dd$D.1.1 + dd$D.3.1 + dd$D.4.1 + dd$D.4.3)/(sdt3*sdt4)
  covt[[ch]] <- data.frame(sdt1, sdt2, sdt3, sdt4, cort12, cort13, cort14, cort23, cort24, cort34)  
}


###################################################
### chunk number 71: summ3
###################################################
beta <- list()
for(ch in 1:nchains){
  beta[[ch]] <- read.table(paste(dirsim[ch], "/beta.sim", sep=""), header=TRUE)
}  


###################################################
### chunk number 72: summ4
###################################################
mixture <- list()
for(ch in 1:nchains){
  mixture[[ch]] <- read.table(paste(dirsim[ch], "/mixmoment.sim", sep=""), header=TRUE)
}  


###################################################
### chunk number 73: summ5
###################################################
library(coda)
mcovb <- list(); mcovt <- list(); mbeta <- list(); mmixture <- list()
for (ch in 1:nchains){
  mcovb[[ch]] <- files2coda(data.frames="covb", thin=1, chain=ch)
  mcovt[[ch]] <- files2coda(data.frames="covt", thin=1, chain=ch)
  mbeta[[ch]] <- files2coda(data.frames="beta", thin=1, chain=ch)  
  mmixture[[ch]] <- files2coda(data.frames="mixture", thin=1, chain=ch)  
}
mcovb <- mcmc.list(mcovb[[1]], mcovb[[2]])
mcovt <- mcmc.list(mcovt[[1]], mcovt[[2]])
mbeta <- mcmc.list(mbeta[[1]], mbeta[[2]])
mmixture <- mcmc.list(mmixture[[1]], mmixture[[2]])


###################################################
### chunk number 74: summ6
###################################################
quant <- c(0, 0.025, 0.5, 0.75, 0.975, 1)
qnames <- paste(quant*100, "\\%", sep = "")

meancovb <- apply(rbind(covb[[1]], covb[[2]]), 2, mean)
meancovt <- apply(rbind(covt[[1]], covt[[2]]), 2, mean)
meanbeta <- apply(rbind(beta[[1]], beta[[2]]), 2, mean)
meanmixture <- apply(rbind(mixture[[1]], mixture[[2]]), 2, mean)

quantcovb <- apply(rbind(covb[[1]], covb[[2]]), 2, quantile, probs=quant)
quantcovt <- apply(rbind(covt[[1]], covt[[2]]), 2, quantile, probs=quant)
quantbeta <- apply(rbind(beta[[1]], beta[[2]]), 2, quantile, probs=quant)
quantmixture <- apply(rbind(mixture[[1]], mixture[[2]]), 2, quantile, probs=quant)

sumcovb <- rbind(meancovb, quantcovb[c("50%", "2.5%", "97.5%"), ]); rownames(sumcovb)[1] <- "Mean"
sumcovt <- rbind(meancovt, quantcovt[c("50%", "2.5%", "97.5%"), ]); rownames(sumcovt)[1] <- "Mean"
sumbeta <- rbind(meanbeta, quantbeta[c("50%", "2.5%", "97.5%"), ]); rownames(sumbeta)[1] <- "Mean"
summixture <- rbind(meanmixture, quantmixture[c("50%", "2.5%", "97.5%"), ]); rownames(summixture)[1] <- "Mean"

sumwant <- cbind(sumcovt, sumbeta, summixture)
print(sumwant)


###################################################
### chunk number 75: summ7
###################################################
regres <- list(); pval <- numeric(9)
for (i in 1:12){
  regres[[i]] <- c(beta[[1]][,i], beta[[2]][,i])
  M <- length(regres[[i]])
  p1 <- sum(regres[[i]] < 0)/M
  p2 <- sum(regres[[i]] > 0)/M
  pval[[i]] <- 2*min(p1, p2)
  names(pval)[i] <- colnames(beta[[1]])[i]
}  
pval <- round(pval, 3)
rm(list="regres")
print(pval)


###################################################
### chunk number 76: summ8
###################################################
betaall <- rbind(beta[[1]], beta[[2]])
dmfeff <- data.frame(
            girl.max4 = betaall$DMF + betaall$GIRL.DMF,
            boy.max4 = betaall$DMF,
            girl.max5 = betaall$DMF + betaall$GIRL.DMF + betaall$DMF.UPPER5,
            boy.max5 = betaall$DMF + betaall$DMF.UPPER5,                     
            girl.man4 = betaall$DMF + betaall$GIRL.DMF + betaall$DMF.LOWER4,
            boy.man4 = betaall$DMF + betaall$DMF.LOWER4,
            girl.man5 = betaall$DMF + betaall$GIRL.DMF + betaall$DMF.LOWER5,
            boy.man5 = betaall$DMF + betaall$DMF.LOWER5                     
          )                     
meandmfeff <- apply(dmfeff, 2, mean)
quantdmfeff <- apply(dmfeff, 2, quantile, probs=quant)
sumdmfeff <- rbind(meandmfeff, quantdmfeff[c("50%", "2.5%", "97.5%"), ]); rownames(sumdmfeff)[1] <- "Mean"
dmfpval <- numeric()
for (i in 1:8){
  M <- length(dmfeff[[i]])
  p1 <- sum(dmfeff[[i]] < 0)/M
  p2 <- sum(dmfeff[[i]] > 0)/M
  dmfpval[[i]] <- 2*min(p1, p2)
  names(dmfpval)[i] <- colnames(dmfeff)[i]
}
print(sumdmfeff)
print(dmfpval)


###################################################
### chunk number 77: dens1
###################################################
ch <- 1
par(mfrow=c(3, 4)); for(i in 1:12){ densplot2(mcmc(beta[[ch]][,i]), main=colnames(beta[[ch]])[i]) }


###################################################
### chunk number 78: dens2
###################################################
par(mfrow=c(3, 4)); for(i in 1:10){ densplot2(mcmc(covt[[ch]][,i]), main=colnames(covt[[ch]])[i]) }


###################################################
### chunk number 79: dens3
###################################################
par(mfrow=c(2, 2)); for(i in 1:3){ densplot2(mcmc(mixture[[ch]][,i]), main=colnames(mixture[[ch]])[i]) }


###################################################
### chunk number 80: acor1
###################################################
par(mfrow=c(3, 4)); autocorr.plot(mcmc(beta[[ch]]), auto.layout=FALSE, ask=FALSE, bty="n")


###################################################
### chunk number 81: acor2
###################################################
par(mfrow=c(3, 4)); autocorr.plot(mcmc(covt[[ch]]), auto.layout=FALSE, ask=FALSE, bty="n")


###################################################
### chunk number 82: acor3
###################################################
par(mfrow=c(2, 2)); autocorr.plot(mcmc(mixture[[ch]]), auto.layout=FALSE, ask=FALSE, bty="n")


###################################################
### chunk number 83: trace1
###################################################
par(mfrow=c(3, 4)); traceplot2(mcmc(beta[[ch]]))


###################################################
### chunk number 84: trace2
###################################################
par(mfrow=c(3, 4)); traceplot2(mcmc(covt[[ch]]))


###################################################
### chunk number 85: trace3
###################################################
par(mfrow=c(2, 2)); traceplot2(mcmc(mixture[[ch]]))


###################################################
### chunk number 86: gelmanRubin
###################################################
gelm.beta <- gelman.diag(mbeta)
gelm.covt <- gelman.diag(mcovt)
gelm.mixture <- gelman.diag(mmixture)
rownames(gelm.beta$psrf) <- dimnames(mbeta[[1]])[[2]]
rownames(gelm.covt$psrf) <- dimnames(mcovt[[1]])[[2]]
rownames(gelm.mixture$psrf) <- dimnames(mmixture[[1]])[[2]]


###################################################
### chunk number 87: pgelman1
###################################################
print(gelm.beta)


###################################################
### chunk number 88: pgelman2
###################################################
print(gelm.covt)


###################################################
### chunk number 89: pgelman3
###################################################
print(gelm.mixture)


###################################################
### chunk number 90: errorDensity1
###################################################
dens <- list()
dens[[1]] <- bayesDensity(dirsim[[1]], grid=seq(1.4, 2, by=0.005), stgrid=seq(-3, 2, by=0.05))
dens[[2]] <- bayesDensity(dirsim[[2]], grid=seq(1.4, 2, by=0.005), stgrid=seq(-3, 2, by=0.05))


###################################################
### chunk number 91: errdens
###################################################
par(mfrow = c(2, 2), bty="n")
plot(dens[[1]], standard=TRUE, dim.plot=FALSE, main="Standardized, chain 1", k.cond=0:10)
plot(dens[[1]], standard=FALSE, dim.plot=FALSE, main="Unstandardized, chain 1", k.cond=0:10)
plot(dens[[2]], standard=TRUE, dim.plot=FALSE, main="Standardized, chain 2", k.cond=0:10)
plot(dens[[2]], standard=FALSE, dim.plot=FALSE, main="Unstandardized, chain 2", k.cond=0:10)


###################################################
### chunk number 92: cleaning
###################################################
files1 <- dir("./tandchain1test")
files2 <- dir("./tandchain2test")
file.remove(paste("./tandchain1test/", files1, sep = ""))
file.remove(paste("./tandchain2test/", files2, sep = ""))
file.remove("tandchain1test")
file.remove("tandchain2test")


