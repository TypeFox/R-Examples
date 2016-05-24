## ## commented out for time-reasons

## ################################################################################
## ## test 1: "validation" using fitMod and vcov.DRMod and predict.DRMod
## model <- "emax"
## sigma <- 1
## n <- c(100,50,50,50,100)
## doses <- c(0,10,20,40,50)
## cf <- c(0,1,10)
## V <- DoseFinding:::aprCov(doses, model, cf, S=diag(1/n))
## DoseFinding:::getPredVar(model, cf, V=V, pDose=50)

## ## now validation using the formulas in fitMod
## doseVec <- rep(doses, n)
## respVec <- rnorm(length(doseVec))
## dd <- fitMod(doseVec, respVec, model="emax")
## ## now change to achieve desired values
## dd$coefs <- cf
## dd$df <- dd$RSS <- sum(n)
## vcov(dd)
## predict(dd, predType = "effect-curve", doseSeq=50, se.fit=TRUE)$se.fit^2

## ################################################################################
## ## test 2: "validation" using simulation
## model <- "emax"
## sigma <- 1
## ## select very large sample size (to validate asymptotics)
## n <- c(100000, 50000, 50000, 50000, 100000)
## n <- c(100, 50, 50, 50, 100)
## doses <- c(0,10,20,40,50)
## cf <- c(0,0.4,10)
## Delta <- 0.2

## V <- DoseFinding:::aprCov(doses, model, cf, S=diag(0.3^2/n))
## tdvar <- DoseFinding:::getTDVar(model, cf, V=V, Delta=Delta, scale = "unrestricted")
## edvar <- DoseFinding:::getEDVar(model, cf, V=V, p=0.5, maxD=50, scale = "unrestricted")
## pavar <-  DoseFinding:::getPredVar(model, cf, V=V, pDose=50)
## tdvart <- DoseFinding:::getTDVar(model, cf, V=V,scale = "log", Delta=Delta)
## edvart <- DoseFinding:::getEDVar(model, cf, V=V,scale = "logit", p=0.5, maxD=50)

## ## simulation
## mn <- emax(doses, cf[1], cf[2], cf[3])
## doseVec <- rep(doses, n)
## mnVec <- rep(mn, n)
## td <- pl <- ed <- numeric(10000)
## for(i in 1:10000){
##   respVec <- mnVec + rnorm(length(mnVec),0,0.3)
##   ff <- fitMod(doseVec, respVec, model="emax", bnds = c(0.05, 75))
##   ed[i] <- ED(ff, p=0.5)
##   td[i] <- TD(ff, Delta = Delta)
##   pl[i] <- predict(ff, doseSeq=50, predType = "effect-curve")
##   pb <- txtProgressBar(min=1, max=1000, char="*", width = 20, style = 3)
##   setTxtProgressBar(pb, i)
## }
## cat("\n")


## mm <- Mods(emax=cf[3], doses=doses, placEff=0, maxEff=emax(50,cf[1],cf[2],cf[3]))
## edt <- ED(mm, p=0.5)
## edt7 <- ED(mm, p=0.7)
## edt3 <- ED(mm, p=0.3)
## tdt <- TD(mm, Delta=Delta)
## hist(td[td < 100], freq=FALSE, breaks = 21)
## curve(dnorm(x, tdt, sqrt(tdvar)), add=TRUE)
## hist(ed, freq=FALSE, breaks = 101)
## curve(dnorm(x, edt, sqrt(edvar)), add=TRUE)
## hist(pl, freq=FALSE, breaks = 21)
## curve(dnorm(x, emax(50,cf[1],cf[2],cf[3]), sqrt(pavar)), add=TRUE)

## hist(td[td < 100], freq=FALSE, breaks = 101)
## curve(dlnorm(x, log(tdt), sqrt(tdvart)), add=TRUE)
## hist(ed, freq=FALSE, breaks = 101)
## ## plot against logit-normal distribution

## mean(ed < edt7 & ed > edt3)
## pnorm(edt7, edt, sqrt(edvar))-pnorm(edt3, edt, sqrt(edvar))

## ################################################################################
## ## test 3: study example
## nSim <- 100
## doses <- c(0,1,3,10,30,50,75,150,300,450)
## n <- c(100,rep(38,8),100)
## sigma <- 380
## mm <- Mods(sigEmax = rbind(c(100,6),c(170,4), c(80,3), c(290,5)),
##            emax = c(5,20,50,120), linear=NULL, doses=doses,
##            placEff=0, maxEff=150)

## model <- "sigEmax"
## pp <- planMod(model, mm, n, sigma, doses=doses,
##               simulation = TRUE, cores = 4,
##               alpha = 0.025, nSim = nSim,
##               p = 0.5, pLB = 0.25, pUB = 0.75)
## print(pp)
## summary(pp, Delta = 130, p = 0.5)
## plot(pp)
## plot(pp, type="ED", 0.5)
## plot(pp, type="TD", Delta = 130, direction = "increasing")

## model <- "emax"
## pp <- planMod(model, mm, n, sigma, doses=doses,
##               simulation = TRUE, cores = 4,
##               alpha = 0.025, 
##               p = 0.5, pLB = 0.25, pUB = 0.75)

## model <- "linear"
## pp <- planMod(model, mm, n, sigma, doses=doses,
##               simulation = TRUE, cores = 4,
##               alpha = 0.025, nSim = nSim,
##               p = 0.5, pLB = 0.25, pUB = 0.75)

## ## now model selection
## model <- c("sigEmax", "emax")
## pp1 <- planMod(model, mm, n, sigma, doses=doses, asyApprox = FALSE,
##               simulation = TRUE,  cores = 4,
##               alpha = 0.025, nSim = nSim,
##               p = 0.5, pLB = 0.25, pUB = 0.75)
## print(pp1)
## summary(pp1)
## plot(pp1)
## plot(pp1, type="ED", 0.5)
## plot(pp1, type="TD", Delta = 130, direction = "increasing")


## ## ################################################################################
## ## ## test 4: study example
## doses <- c(0,10,25,50,100,150)
## fmodels <- Mods(linear = NULL, emax = 25,
##                 logistic = c(50, 10.88111), exponential= 85,
##                 betaMod=rbind(c(0.33,2.31),c(1.39,1.39)),
##                 doses = doses, addArgs=list(scal = 200),
##                 placEff = 0, maxEff = 0.4)
## sigma <- 1
## n <- rep(62, 6)*2

## ## use all models not used previously
## model <- c("linear", "quadratic", "exponential", "betaMod", "logistic", "linlog")
## altb <- defBnds(200)
## pp <- planMod(model, fmodels, n, sigma, doses=doses,
##               asyApprox = FALSE, simulation = TRUE, cores = 4,
##               alpha = 0.025, nSim = nSim, bnds=altb,
##               p = 0.5, pLB = 0.25, pUB = 0.75)
## pp
## summary(pp, p = 0.5, Delta = 0.3)
## plot(pp)
## plot(pp, type = "TD", Delta=0.3, direction = "i")
## plot(pp, type = "ED", p = 0.5)
