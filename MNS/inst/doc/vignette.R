### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignette.Rnw:395-404 (eval = FALSE)
###################################################
## library("MNS")
## set.seed(1)
## N=3
## Net = gen.Network(method = "cohort", p = 20, 
##                        Nsub = N, sparsity = .2, 
##                        REsize = 20, REprob = .5,
##                        REnoise = 1)
## # plot simulated networks:
## plot(Net, view="sub")


###################################################
### code chunk number 2: vignette.Rnw:447-452 (eval = FALSE)
###################################################
## library("MNS")
## set.seed(1)
## Net = gen.Network(method = "Danaher", p = 20)
## # plot simulated networks:
## plot(Net, view="sub")


###################################################
### code chunk number 3: vignette.Rnw:468-478 (eval = FALSE)
###################################################
## library("MNS")
## set.seed(1)
## N=3
## Net = gen.Network(method = "cohort", p = 20, 
##                        Nobs = 500,
##                        Nsub = N, sparsity = .2, 
##                        REsize = 20, REprob = .5,
##                        REnoise = 1)
## # plot simulated networks:
## plot.ts(Net$Data[[1]][, c(1,2,3,4,5, 6)], main="")


###################################################
### code chunk number 4: vignette.Rnw:585-600 (eval = FALSE)
###################################################
## set.seed(1)
## N=10
## Net = gen.Network(method = "cohort", p = 10, 
##                         Nsub = N, sparsity = .2, 
##                         REsize = 10, REprob = .75,
##                         REnoise = 1, Nobs = 75)
## # run cross-validation 
## CV = cv.MNS(dat = Net$Data, 
##             l1range = seq(.05, .075, length.out = 5), 
##             alpharange = seq(.25, .75, length.out = 3),
##             K=5, parallel=TRUE)
## # fit MNS model:
## mns = MNS(dat = Net$Data, 
##           lambda_pop = CV$l1 * (1-CV$alpha), 
##           lambda_random = CV$l1 * (CV$alpha))


###################################################
### code chunk number 5: vignette.Rnw:620-621 (eval = FALSE)
###################################################
## plot(mns, view="pop")


###################################################
### code chunk number 6: vignette.Rnw:635-636 (eval = FALSE)
###################################################
## plot(mns, view="var")


###################################################
### code chunk number 7: vignette.Rnw:675-676 (eval = FALSE)
###################################################
## plot(mns, view="sub", subID=c(2,4,6))


