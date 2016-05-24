### R code from vignette source 'UsingSpikeSlabGAM.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE,
        digits=5, mc.cores=2)
clrs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
        "#A65628", "#F781BF", "#999999")


###################################################
### code chunk number 2: loadLibs
###################################################
library("spikeSlabGAM")
library("mboost")
library("ggplot2")
library("gtable")


###################################################
### code chunk number 3: mkData
###################################################
set.seed(1312424)
n <- 200
snr <- 3

sm1 <- runif(n)
fsm1 <- dbeta(sm1, 7, 3)/2

sm2 <- runif(n, 0, 1)
f <- gl(3, n/3)
ff <- as.numeric(f)/2
fsm2f <- ff + ff * sm2 +
        ((f == 1) * -dbeta(sm2, 6, 4) + (f == 2) * dbeta(sm2, 6, 9) +
           (f == 3) * dbeta(sm2, 9, 6))/2

lin <- matrix(rnorm(n * 3), n, 3)
colnames(lin) <- paste("lin", 1:3, sep = "")

noise1 <- sm1 + rnorm(n)
noise2 <- runif(n)
noise3 <- runif(n)
noise4 <- sample(gl(4, n/4))

eta <- drop(fsm1 + fsm2f + lin %*% c(.1, .2, .3))
y <- eta + sd(eta)/snr * rt(n, df = 5)
d <- data.frame(y, sm1, sm2, f, lin, noise1, noise2, noise3, noise4)


###################################################
### code chunk number 4: defForm1
###################################################
f1 <- y ~ (sm1 + sm2 + f + lin1)^2 + lin2 +lin3 + noise1 + noise2 + noise3 + noise4


###################################################
### code chunk number 5: fit1
###################################################
m <- spikeSlabGAM(formula = f1, data = d)


###################################################
### code chunk number 6: summary1.1Fake (eval = FALSE)
###################################################
## summary(m)


###################################################
### code chunk number 7: summary1.1Prep
###################################################
summ <- summary(m)


###################################################
### code chunk number 8: summary1.1
###################################################
print(summ, printModels=F)


###################################################
### code chunk number 9: summary1.2
###################################################
#this is just ctrl-c-v from summary.spikeSlabGAM s.t. we don't get too much
#redundant info about model formula etc. roughly the same as: print(sum,
#printPGamma=F, printModels=T)
cat("\nPosterior model probabilities (inclusion threshold =",summ$thresh,"):\n")
modelTable <- {
    #make ("x","")-vectors out of model names
    models <- sapply(names(summ$modelTable), function(x){
                sapply(strsplit(gsub("1", "x", x),""),
                        function(y) gsub("0", "", y))
            })
    models <- rbind(round(summ$modelTable,3), models,
                    round(cumsum(summ$modelTable), 3))
    rownames(models) <- c("prob.:",
                          names(summ$predvars[!grepl("u(", names(summ$predvars),
                                                    fixed=TRUE)]),
                          "cumulative:")
    models <- data.frame(models)
    showModels <- 8
    models <- models[,1:showModels, drop=F]
    colnames(models) <- 1:NCOL(models)
    models
}
print(modelTable)


###################################################
### code chunk number 10: plotm1Fake (eval = FALSE)
###################################################
## plot(m)


###################################################
### code chunk number 11: plotm1prep
###################################################
pl <- plot(m, nrow=4, ncol=3,
  ggElems=list(theme(legend.key.size = grid:::unit(0.01, "npc"),
    legend.text = element_text(size = 7),
    legend.title=element_text(size=8, hjust = 0),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8, hjust=1),
    plot.margin=grid:::unit(c(0.5, 0, 0, 0), "lines"))))


###################################################
### code chunk number 12: plotm1
###################################################
do.call(gridExtra:::grid.arrange, c(pl, nrow=4, ncol=3))


###################################################
### code chunk number 13: plotm1sm1Fake (eval = FALSE)
###################################################
## plot(m, labels = c("lin(sm1)", "sm(sm1)"), cumulative = FALSE)
## trueFsm1 <- data.frame(truth = fsm1 - mean(fsm1), sm1 = sm1)
## plot(m, labels = "sm(sm1)", ggElems = list(geom_line(aes(x = sm1, y = truth),
##                                                      data = trueFsm1, linetype = 2)))


###################################################
### code chunk number 14: plotm1sm1Prep
###################################################
pSm11 <- plot(m, labels="sm(sm1)",
  ggElems=list(theme(axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9, hjust=1))))
pSm1Sep <-  plot(m, labels=c("lin(sm1)", "sm(sm1)"), cumulative=FALSE,
  ggElems=list(theme(axis.text.x = element_text(size=9),
    axis.text.y = element_text(size=9, hjust=1))),
  nrow=1, ncol=2) #,
# JSS version only:        ggElems=list(list(ylab("$\\eta$"))))


###################################################
### code chunk number 15: plotm1sm1
###################################################
trueFsm1 <-data.frame(truth=fsm1-mean(fsm1), sm1=sm1)
pSm1Cum <- pSm11[[1]] + geom_line(data=trueFsm1, aes(x=sm1, y=truth), lty=2)  +
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9, hjust=1))# + ylab("$\\eta$")
gridExtra:::grid.arrange(pSm1Sep[[1]], pSm1Sep[[2]], pSm1Cum, nrow=1)


###################################################
### code chunk number 16: plotm1sm2fFake (eval = FALSE)
###################################################
## trueFsm2f <-data.frame(truth = fsm2f - mean(fsm2f), sm2 = sm2, f = f)
## plot(m, labels = "sm(sm2):fct(f)",
##      ggElems = list(geom_line(aes(x = sm2, y = truth, colour = f),
##                             data = trueFsm2f, linetype = 2)))


###################################################
### code chunk number 17: plotm1sm2fprep
###################################################
trueFsm2f <-data.frame(truth=fsm2f-mean(fsm2f), sm2=sm2, f=f)
#JSS only: plot(m, labels="sm(sm2):fct(f)", ggElems=list(list(ylab("$\\eta$")), geom_line(aes(x=sm2, y=truth, colour=f), data=trueFsm2f, linetype=2)))
pl2 <- plot(m, labels="sm(sm2):fct(f)",
  ggElems=list(geom_line(aes(x=sm2, y=truth, colour=f), data=trueFsm2f,
    linetype=2)))


###################################################
### code chunk number 18: plotm1sm2f
###################################################
pl2[[1]]


###################################################
### code chunk number 19: pimaDataPrep
###################################################
data("PimaIndiansDiabetes2", package = "mlbench")
pimaDiab <- na.omit(PimaIndiansDiabetes2[, -c(4, 5)])
pimaDiab <- within(pimaDiab,{
            diabetes <- 1*(diabetes == "pos")
        })
set.seed(1109712439)
testInd <- sample(1:nrow(pimaDiab), 200)
pimaDiabTrain <- pimaDiab[-testInd,]


###################################################
### code chunk number 20: pimaM0
###################################################
mcmc <- list(nChains = 8, chainLength = 1000, burnin = 500, thin = 5)
m0 <- spikeSlabGAM(diabetes ~ pregnant + glucose + pressure + mass + pedigree + age,
        family = "binomial", data = pimaDiabTrain, mcmc = mcmc)


###################################################
### code chunk number 21: pimaM0Summary
###################################################
pr0 <- predict(m0, newdata = pimaDiab[testInd,])
print(summary(m0), printModels = FALSE)


###################################################
### code chunk number 22: pimaB0
###################################################
b <- gamboost(as.factor(diabetes) ~ pregnant + glucose + pressure + mass + pedigree + age,
        family = Binomial(), data = pimaDiabTrain)[300]
aic <- AIC(b, method = "classical")
prB <- predict(b[mstop(aic)], newdata = pimaDiab[testInd,])


###################################################
### code chunk number 23: pimaB0Sum
###################################################
summary(b[mstop(aic)])$selprob


###################################################
### code chunk number 24: pimaM0Pred
###################################################
dev <- function(y, p){
    -2*sum(dbinom(x = y, size = 1, prob = p, log = T))
}
c(spikeSlabGAM = dev(pimaDiab[testInd, "diabetes"], pr0),
    gamboost   = dev(pimaDiab[testInd, "diabetes"], plogis(prB)))


###################################################
### code chunk number 25: pimaM1
###################################################
hyper1 <- list(gamma = c(v0 = 0.005))
m1 <- spikeSlabGAM(diabetes ~ pregnant + glucose + pressure + mass + pedigree + age,
                   family = "binomial", data = pimaDiabTrain, mcmc = mcmc,
                   hyperparameters = hyper1)
pr1 <- predict(m1, newdata = pimaDiab[testInd, ])


###################################################
### code chunk number 26: pimaM1Sum
###################################################
print(summary(m1), printModels = FALSE)
(dev(pimaDiab[testInd, "diabetes"], pr1))


###################################################
### code chunk number 27: sessionInfo
###################################################
	sessionInfo()


