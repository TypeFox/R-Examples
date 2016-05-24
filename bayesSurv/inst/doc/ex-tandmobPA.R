###################################################
### chunk number 1: directories
###################################################
anadir <- "/home/komari/win/work/papers/gsplineTand/anaDMFsymm/"
chaindir.PA <- paste(anadir, "chains/modelPA", c("1646", "2636"), sep="")
initdir.PA <- paste(anadir, "chains/end_PA", c("1646", "2636"), sep="")
plotdir.PA <- paste(anadir, "summPlots/modelPA/", sep="")
resultdir <- paste(anadir, "results/", sep="")
predCurvesdir <- paste(resultdir, "predCurves/", sep="")
figuredir <- "/home/komari/win/work/papers/gsplineTand/RforCRAN/figuresPA/"


###################################################
### chunk number 2: loadLibData
###################################################
library(bayesSurv)
data(tandmobRoos)


###################################################
### chunk number 3: data01
###################################################
pair <- list(pair1=c(16, 46), pair2=c(26, 36))
dmfp <- list()
for (pii in 1:length(pair)){
  tt1 <- pair[[pii]][1]
  tt2 <- pair[[pii]][2]
  remove <- (is.na(tandmobRoos[, paste("FBEG.", tt1, sep="")])  |
            is.na(tandmobRoos[, paste("FBEG.", tt2, sep="")])  |
            (tandmobRoos[, paste("TOOTH.", tt1, sep="")] == 0)  |
            (tandmobRoos[, paste("TOOTH.", tt2, sep="")] == 0)  |
            is.na(tandmobRoos[, paste("T", tt1+39, "d", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt1+39, "m", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt1+39, "f", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt1+39, "s", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt2+39, "d", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt2+39, "m", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt2+39, "f", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt2+39, "s", sep="")])  |  
            is.na(tandmobRoos[, paste("SEAL.", tt1, sep="")])  |
            is.na(tandmobRoos[, paste("SEAL.", tt2, sep="")])  |
            is.na(tandmobRoos[, "FREQ.BR"])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt1, ".1", sep="")])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt1, ".2", sep="")])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt2, ".1", sep="")])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt2, ".2", sep="")]))
  dmfp[[pii]] <- tandmobRoos[!remove,]
}
names(dmfp) <- c("1646", "2636")
n.sample <- sapply(dmfp, nrow)
print(n.sample)


###################################################
### chunk number 4: data02
###################################################
time.0 <- 5


###################################################
### chunk number 5: data03
###################################################
Ebeg <- list()
Eend <- list()
for (pii in 1:length(pair)){
  tt1 <- pair[[pii]][1]
  tt2 <- pair[[pii]][2]  
  Ebeg[[pii]] <- cbind(dmfp[[pii]][, paste("EBEG.", tt1, sep="")], dmfp[[pii]][, paste("EBEG.", tt2, sep="")])
  Eend[[pii]] <- cbind(dmfp[[pii]][, paste("EEND.", tt1, sep="")], dmfp[[pii]][, paste("EEND.", tt2, sep="")])
  Ebeg[[pii]] <- as.vector(t(Ebeg[[pii]])) - time.0
  Eend[[pii]] <- as.vector(t(Eend[[pii]])) - time.0
  Ebeg[[pii]][Ebeg[[pii]] <= 0] <- NA
}


###################################################
### chunk number 6: data04
###################################################
Fbeg <- list()
Fend <- list()
for (pii in 1:length(pair)){
  tt1 <- pair[[pii]][1]
  tt2 <- pair[[pii]][2]  
  Fbeg[[pii]] <- cbind(dmfp[[pii]][, paste("FBEG.", tt1, sep="")], dmfp[[pii]][, paste("FBEG.", tt2, sep="")])
  Fend[[pii]] <- cbind(dmfp[[pii]][, paste("FEND.", tt1, sep="")], dmfp[[pii]][, paste("FEND.", tt2, sep="")])
  Fbeg[[pii]] <- as.vector(t(Fbeg[[pii]])) - time.0
  Fend[[pii]] <- as.vector(t(Fend[[pii]])) - time.0
}


###################################################
### chunk number 7: data05
###################################################
Idnr <- list()
Maxilla <- list()
Girl <- list()
Gender <- list()
Brush <- list()
fBrush <- list()
Plaque <- list()
Plaque.1 <- list()
Plaque.2 <- list()
fPlaque <- list()
Seal <- list()
fSeal <- list()
Prim5 <- list()
Prim5d <- list()
Prim5m <- list()
Prim5f <- list()
fPrim5 <- list()
for (pii in 1:length(pair)){
  tt1 <- pair[[pii]][1]
  tt2 <- pair[[pii]][2]  
  Idnr[[pii]] <- rep(dmfp[[pii]][, "IDNR"], rep(2, dim(dmfp[[pii]])[1]))
  Maxilla[[pii]]  <- rep(c(1, 0), dim(dmfp[[pii]])[1])

  ## dummies
  Girl[[pii]]     <- rep(dmfp[[pii]][, "GIRL"], rep(2, dim(dmfp[[pii]])[1]))
  Brush[[pii]]    <- rep(dmfp[[pii]][, "FREQ.BR"], rep(2, dim(dmfp[[pii]])[1])) 
  Plaque.1[[pii]] <- as.vector(t(cbind(dmfp[[pii]][, paste("PLAQUE.", tt1, ".1", sep="")], 
                                       dmfp[[pii]][, paste("PLAQUE.", tt2, ".1", sep="")])))
  Plaque.2[[pii]] <- as.vector(t(cbind(dmfp[[pii]][, paste("PLAQUE.", tt1, ".2", sep="")], 
                                       dmfp[[pii]][, paste("PLAQUE.", tt2, ".2", sep="")])))
  Seal[[pii]]     <- as.vector(t(cbind(dmfp[[pii]][, paste("SEAL.", tt1, sep="")],
                                       dmfp[[pii]][, paste("SEAL.", tt2, sep="")])))
  Prim5d[[pii]]   <- as.vector(t(cbind(dmfp[[pii]][, paste("T", tt1+39, "d", sep="")], 
                                       dmfp[[pii]][, paste("T", tt2+39, "d", sep="")])))
  Prim5m[[pii]]   <- as.vector(t(cbind(dmfp[[pii]][, paste("T", tt1+39, "m", sep="")], 
                                       dmfp[[pii]][, paste("T", tt2+39, "m", sep="")])))
  Prim5f[[pii]]   <- as.vector(t(cbind(dmfp[[pii]][, paste("T", tt1+39, "f", sep="")], 
                                       dmfp[[pii]][, paste("T", tt2+39, "f", sep="")])))

  ## factors
  Gender[[pii]] <- factor(1*(Girl[[pii]]==1), levels=0:1, labels=c("boy", "girl"))
  fBrush[[pii]] <- factor(1*(Brush[[pii]]==1), levels=0:1, labels=c("not.daily", "daily"))  
  fPlaque[[pii]] <- ordered(0*(Plaque.1[[pii]]==0 & Plaque.2[[pii]]==0) +
                            1*(Plaque.1[[pii]]==1 & Plaque.2[[pii]]==0) +
                            2*(Plaque.1[[pii]]==0 & Plaque.2[[pii]]==1), levels=0:2, labels=c("none", "pits.fiss", "total"))
  fSeal[[pii]] <- factor(1*(Seal[[pii]]==1), levels=0:1, labels=c("no", "yes"))  
  fPrim5[[pii]] <- factor(0*(Prim5d[[pii]]==0 & Prim5f[[pii]]==0 & Prim5m[[pii]]==0) +
                          1*(Prim5d[[pii]]==1 & Prim5f[[pii]]==0 & Prim5m[[pii]]==0) +
                          2*(Prim5d[[pii]]==0 & Prim5f[[pii]]==1 & Prim5m[[pii]]==0) +
                          3*(Prim5d[[pii]]==0 & Prim5f[[pii]]==0 & Prim5m[[pii]]==1), levels=0:3,
                          labels=c("sound", "decayed", "filled", "missing"))

  ## simplified covariates for Plaque and Prim5
  Plaque[[pii]] <- 0*(fPlaque[[pii]]=="none") + 1*(fPlaque[[pii]]=="pits.fiss" | fPlaque[[pii]]=="total")
  Prim5[[pii]] <- 0*(fPrim5[[pii]]=="sound") + 1*(fPrim5[[pii]]=="decayed" | fPrim5[[pii]]=="filled" | fPrim5[[pii]]=="missing")    
}


###################################################
### chunk number 8: data06
###################################################
data.PA <- list()
for (pii in 1:length(pair)){
  data.PA[[pii]] <- data.frame(Idnr=Idnr[[pii]],
                               Ebeg=Ebeg[[pii]],
                               Eend=Eend[[pii]],
                               Fbeg=Fbeg[[pii]],
                               Fend=Fend[[pii]], 
                               Maxilla=Maxilla[[pii]],
                               Girl=Girl[[pii]],
                               Gender=Gender[[pii]],
                               Brush=Brush[[pii]],
                               fBrush=fBrush[[pii]],
                               Plaque=Plaque[[pii]],
                               Plaque.1=Plaque.1[[pii]],
                               Plaque.2=Plaque.2[[pii]],
                               fPlaque=fPlaque[[pii]],
                               Seal=Seal[[pii]],
                               fSeal=fSeal[[pii]],
                               Prim5=Prim5[[pii]],
                               Prim5d=Prim5d[[pii]],
                               Prim5m=Prim5m[[pii]],
                               Prim5f=Prim5f[[pii]],
                               fPrim5=fPrim5[[pii]]                               
                               )
}
names(data.PA) <- c("1646", "2636")
rm(list=c("dmfp", "Idnr", "Ebeg", "Eend", "Fbeg", "Fend", "Maxilla", "Girl", "Gender", "Brush", "fBrush",
          "Plaque", "Plaque.1", "Plaque.2", "fPlaque", "Seal", "fSeal",
          "Prim5", "Prim5d", "Prim5m", "Prim5f", "fPrim5"))

rm(list=c("pii", "remove", "tt1", "tt2"))


###################################################
### chunk number 9: data07
###################################################
print(data.PA$"1646"[1:10,])


###################################################
### chunk number 10: initAFT01
###################################################
fit.aft <- function(data){
  vbeg <- data[, "Ebeg"]
  vend <- data[, "Eend"]
  vbeg[vbeg <= 0] <- NA               ## this is left-censored

  tfit <- survreg(Surv(vbeg, vend, type="interval2")~Maxilla+Girl, dist="loglogistic", data=data)
  return(tfit)
}
emerg.PA.aft <- list()
emerg.PA.aft$"1646" <- fit.aft(data=data.PA$"1646")
emerg.PA.aft$"2636" <- fit.aft(data=data.PA$"2636")


###################################################
### chunk number 11: initAFT02
###################################################
## To get some impression on the distribution of caries we ad hod impute emergence times as 
## either midpoints of intervals or a~midpoint between 6 and left-censored time or equal to the right-censored time
fit.aft2 <- function(data){
  vebeg <- data[, "Ebeg"]
  veend <- data[, "Eend"]
  vfbeg <- data[, "Fbeg"]
  vfend <- data[, "Fend"]
  vebeg[vebeg <= 0] <- NA         ## this is left censored

  left <- is.na(vebeg)
  interv <- !is.na(vebeg) & !is.na(veend)
  right <- is.na(veend)

  onset <- vebeg
  onset[left] <- 0.5*(0 + veend[left])
  onset[right] <- vebeg[right]
  onset[interv] <- 0.5*(vebeg[interv] + veend[interv])
  vfbeg2 <- vfbeg - onset
  vfbeg2[vfbeg2 < 0] <-  NA          ### -> both emergence and caries were in one interval, make it left censored now
  vfend2 <- vfend - onset

  tfit <- survreg(Surv(vfbeg2, vfend2, type="interval2")~Maxilla+Girl+Brush+Plaque+Seal+Prim5,
                  dist="loglogistic", data=data)
  return(tfit)
}
caries.PA.aft <- list()
caries.PA.aft$"1646" <- fit.aft2(data=data.PA$"1646")
caries.PA.aft$"2636" <- fit.aft2(data=data.PA$"2636")


###################################################
### chunk number 12: initAFT03
###################################################
lapply(emerg.PA.aft, summary)


###################################################
### chunk number 13: initAFT04
###################################################
lapply(caries.PA.aft, summary)


###################################################
### chunk number 14: prior01
###################################################
prior.PA.gspl.emerg <- list(
  specification   = 2, 
  K               = c(15, 15),  
  izero           = c(0, 0),
  neighbor.system = "uniCAR",
  order           = 3,
  equal.lambda    = FALSE,
  c4delta         = c(1.5, 1.5),
  prior.lambda    = c("gamma", "gamma"),
  prior.intercept = c("normal", "normal"),
  prior.scale     = c("gamma", "gamma"),
  prior.gamma     = c("fixed", "fixed"),
  prior.sigma     = c("fixed", "fixed"),
  shape.lambda    = c(1, 1),
  rate.lambda     = c(0.005, 0.005),
  mean.intercept  = c(0, 0),
  var.intercept   = c(100, 100),
  shape.scale     = c(1, 1),
  rate.scale      = c(0.005, 0.005)
)
prior.PA.gspl.caries <- prior.PA.gspl.emerg


###################################################
### chunk number 15: prior02
###################################################
prior.PA.beta.emerg <- list(mean.prior=0, var.prior=1e2)
prior.PA.beta.caries <- list(mean.prior=rep(0, 5), var.prior=rep(1e2, 5))


###################################################
### chunk number 16: init01
###################################################
iinit.PA.emerg <- list(
  lambda    = c(3000, 3000),
  intercept = c(0.40, 0.42),
  scale     = rep(0.10, 2),
  gamma     = c(0, 0),
  sigma     = c(0.2, 0.2),
  beta      = -0.01
)
init.PA.emerg <- list(iinit.PA.emerg, iinit.PA.emerg)
names(init.PA.emerg) <- c("1646", "2636")


###################################################
### chunk number 17: init02
###################################################
iinit.PA.caries <- list(
  lambda    = c(3000, 3000),
  intercept = c(2.20, 2.20),
  scale     = c(0.55, 0.55),
  gamma     = c(0, 0),
  sigma     = c(0.2, 0.2),
  beta      = c(-0.07, 0.32, -0.20, 0.04, -0.61)
)
init.PA.caries <- list(iinit.PA.caries, iinit.PA.caries)
names(init.PA.caries) <- c("1646", "2636")


###################################################
### chunk number 18: init03
###################################################
init.PA.emerg <- list()
init.PA.caries <- list()

for (k in 1:2){  
  init.PA.emerg[[k]] <- list()
    init.PA.emerg[[k]]$iter      <- scan(paste(initdir.PA[k], "/iteration.sim", sep=""), skip=1)
    init.PA.emerg[[k]]$lambda    <- scan(paste(initdir.PA[k], "/lambda.sim", sep=""), skip=1)
    init.PA.gspline         <- scan(paste(initdir.PA[k], "/gspline.sim", sep=""), skip=1)
    init.PA.emerg[[k]]$gamma     <- init.PA.gspline[1:2]
    init.PA.emerg[[k]]$sigma     <- init.PA.gspline[3:4]
    init.PA.emerg[[k]]$intercept <- init.PA.gspline[7:8]
    init.PA.emerg[[k]]$scale     <- init.PA.gspline[9:10]
    init.PA.emerg[[k]]$beta      <- scan(paste(initdir.PA[k], "/beta.sim", sep=""), skip=1)
    init.PA.emerg[[k]]$a         <- scan(paste(initdir.PA[k], "/mlogweight.sim", sep=""), skip=0)
    init.PA.emerg[[k]]$y         <- matrix(scan(paste(initdir.PA[k], "/Y.sim", sep=""), skip=0), ncol=2, byrow=TRUE)
    init.PA.r <- scan(paste(initdir.PA[k], "/r.sim", sep=""), skip=0)
    init.PA.emerg[[k]]$r <- vecr2matr(init.PA.r, prior.PA.gspl.emerg$K)

  init.PA.caries[[k]] <- list()
    init.PA.caries[[k]]$iter      <- scan(paste(initdir.PA[k], "/iteration.sim", sep=""), skip=1)
    init.PA.caries[[k]]$lambda    <- scan(paste(initdir.PA[k], "/lambda_2.sim", sep=""), skip=1)
    init.PA.gspline         <- scan(paste(initdir.PA[k], "/gspline_2.sim", sep=""), skip=1)
    init.PA.caries[[k]]$gamma     <- init.PA.gspline[1:2]
    init.PA.caries[[k]]$sigma     <- init.PA.gspline[3:4]
    init.PA.caries[[k]]$intercept <- init.PA.gspline[7:8]
    init.PA.caries[[k]]$scale     <- init.PA.gspline[9:10]
    init.PA.caries[[k]]$beta      <- scan(paste(initdir.PA[k], "/beta_2.sim", sep=""), skip=1)
    init.PA.caries[[k]]$a         <- scan(paste(initdir.PA[k], "/mlogweight_2.sim", sep=""), skip=0)
    init.PA.caries[[k]]$y         <- matrix(scan(paste(initdir.PA[k], "/Y_2.sim", sep=""), skip=0), ncol=2, byrow=TRUE)
    init.PA.r <- scan(paste(initdir.PA[k], "/r_2.sim", sep=""), skip=0)
    init.PA.caries[[k]]$r <- vecr2matr(init.PA.r, prior.PA.gspl.caries$K)
}
names(init.PA.emerg) <- c("1646", "2636")
names(init.PA.caries) <- c("1646", "2636")  


###################################################
### chunk number 19: init04
###################################################
str(init.PA.emerg$"1646")


###################################################
### chunk number 20: init05
###################################################
str(init.PA.caries$"1646")


###################################################
### chunk number 21: sample01
###################################################
nsimul.PA <- list(niter=250000, nthin=3, nburn=225000, nwrite=100)


###################################################
### chunk number 22: sample02 eval=FALSE
###################################################
## ##dyn.load("/home/komari/Rlib/bayesSurvival/bayessurvreg02/src/bayessurvreg.so")
## ##source("/home/komari/Rlib/bayesSurvival/bayessurvreg02/Rextra/toSource.R")
## sides <- c("RIGHT", "LEFT")
## sample.PA <- list()
## for (k in 1:2){
##   cat("\nPerforming ", sides[k], " part of the mouth\n", sep="")
##   cat("===========================================\n")
##   
##   data.now <- data.PA[[k]]
##   sample.PA[[k]] <- bayesBisurvreg(
##     formula=Surv(Ebeg, Eend, type="interval2")~ Girl + cluster(Idnr),
##     formula2=Surv(Fbeg, Fend, type="interval2")~ Girl + Brush + Plaque + Seal + Prim5+ cluster(Idnr),
##     onlyX=FALSE,
##     dir=chaindir.PA[k], nsimul=nsimul.PA,
##     prior=prior.PA.gspl.emerg, prior2=prior.PA.gspl.caries,
##     prior.beta=prior.PA.beta.emerg, prior.beta2=prior.PA.beta.caries,
##     init=init.PA.emerg[[k]], init2=init.PA.caries[[k]],
##     store=list(a=FALSE, a2=FALSE),                                   
##     data=data.now)
## }


###################################################
### chunk number 23: predict01
###################################################
skip <- 0
nwrite <- 5000
pred.grid <- c(seq(0.1, 1.5, by=0.05), seq(1.65, 6, by=0.15))
quants <- c(0.025, 0.5, 0.975)
nquant <- length(quants)


###################################################
### chunk number 24: predict02
###################################################
pGender <- factor(c("boy", "girl"))
pBrush <- factor(c("not.daily", "daily"))
pSeal <- factor(c("no", "yes"))
pPlaque <- factor(c("none", "present"))
pPrim5 <- factor(c("sound", "dmf"))
pdata <- expand.grid(ffPrim5=pPrim5, ffPlaque=pPlaque, fSeal=pSeal, fBrush=pBrush, Gender=pGender)
pdata$Prim5 <- 1*(pdata$ffPrim5=="dmf")
pdata$Plaque <- 1*(pdata$ffPlaque=="present")
pdata$Seal <- 1*(pdata$fSeal=="yes")
pdata$Brush <- 1*(pdata$fBrush=="daily")
pdata$Girl <- 1*(pdata$Gender=="girl")
ncov <- nrow(pdata)

pred.data <- data.frame(Idnr=1:(2*ncov), Maxilla=c(rep(1, ncov), rep(0, ncov)), Dimension=c(rep(1, ncov), rep(2, ncov)))
pred.data <- cbind(pred.data, rbind(pdata, pdata))


###################################################
### chunk number 25: predict03
###################################################
print(pred.data)


###################################################
### chunk number 26: predict04
###################################################
side <- c("RIGHT", "LEFT")
jaw <- c("MAXILLARY", "MANDIBULAR")
start <- list(maxilla=seq(1, ncov, by=2), mandible=seq(ncov+1, 2*ncov, by=2))
end <- list(maxilla=seq(2, ncov, by=2), mandible=seq(ncov+2, 2*ncov, by=2))


###################################################
### chunk number 27: predict05 eval=FALSE
###################################################
## pred.PA <- list()
## for (k in 1:2){    ## loop over right and left
##   cat("Performing ", side[k], " side\n", sep="")
##   
##   for (jj in 1:2){  ## loop over maxilla and mandible
##     cat("Performing ", jaw[jj], " jaw\n", sep="")    
##     pred.PA[[(k-1)*2 + jj]] <- list()
##     
##     for (ii in 1:length(start[[jj]])){  ## loop over sets of covariates      
##       cat("Performing the covariate set number ", ii, "\n", sep="")
##       pdata.now <- pred.data[start[[jj]][ii]:end[[jj]][ii], ]
##       nr <- nrow(pdata.now)
##       pred.PA[[(k-1)*2+jj]][[ii]] <- 
##           predictive2(Surv(rep(1, nr), rep(1, nr)) ~ Girl+Brush+Plaque+Seal+Prim5+ cluster(Idnr),
##                       obs.dim=pdata.now$Dimension, grid=pred.grid, data=pdata.now,
##                       Gspline=list(dim=2, K=c(15, 15)), quantile=quants,
##                       skip=skip, by=1, nwrite=nwrite, only.aver=FALSE,
##                       predict=list(density=TRUE, Surv=TRUE, hazard=TRUE),
##                       dir=chaindir.PA[k], extens="_2")
##     }
##   }
## }  
## tnames <- c("16", "46", "26", "36")
## names(pred.PA) <- tnames


###################################################
### chunk number 28: predict06 eval=FALSE
###################################################
## predd.PA <- list()
## for (tt in 1:4){    ## loop over teeth
##   predd.PA[[tt]] <- list()
##   predd.PA[[tt]]$grid <- pred.PA[[tt]][[1]]$grid
##   
##   predd.PA[[tt]]$Surv <- pred.PA[[tt]][[1]]$Surv
##   predd.PA[[tt]]$hazard <- pred.PA[[tt]][[1]]$hazard
##   predd.PA[[tt]]$density <- pred.PA[[tt]][[1]]$density
## 
##   predd.PA[[tt]]$quant.Surv <- pred.PA[[tt]][[1]]$quant.Surv
##   predd.PA[[tt]]$quant.hazard <- pred.PA[[tt]][[1]]$quant.hazard
##   predd.PA[[tt]]$quant.density <- pred.PA[[tt]][[1]]$quant.density
##   
##   for (ii in 2:length(start[[1]])){
##     predd.PA[[tt]]$Surv <- rbind(predd.PA[[tt]]$Surv, pred.PA[[tt]][[ii]]$Surv)
##     predd.PA[[tt]]$hazard <- rbind(predd.PA[[tt]]$hazard, pred.PA[[tt]][[ii]]$hazard)    
##     predd.PA[[tt]]$density <- rbind(predd.PA[[tt]]$density, pred.PA[[tt]][[ii]]$density)
## 
##     templ <- length(predd.PA[[tt]]$quant.Surv)
##     for (ij in 1:length(pred.PA[[tt]][[ii]]$quant.Surv)){    
##       predd.PA[[tt]]$quant.Surv[[templ + ij]] <- pred.PA[[tt]][[ii]]$quant.Surv[[ij]]
##       predd.PA[[tt]]$quant.hazard[[templ + ij]] <- pred.PA[[tt]][[ii]]$quant.hazard[[ij]]
##       predd.PA[[tt]]$quant.density[[templ + ij]] <- pred.PA[[tt]][[ii]]$quant.density[[ij]]      
##     }
##     
##   }  
## }  
## names(predd.PA) <- tnames
## pred.PA <- predd.PA
## rm(list="predd.PA")


###################################################
### chunk number 29: predict07 eval=FALSE
###################################################
## for (tt in 1:4){
##   sink(paste(predCurvesdir, "predDensCaries.", tnames[tt], ".PA.", "dat", sep=""))
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")
##   for (j in 1:ncov) cat(pred.PA[[tt]]$density[j,], "\n", sep="   ")
##   sink()
## 
##   sink(paste(predCurvesdir, "predSurvCaries.", tnames[tt], ".PA.", "dat", sep=""))
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")
##   for (j in 1:ncov) cat(pred.PA[[tt]]$Surv[j,], "\n", sep="   ")
##   sink()
## 
##   sink(paste(predCurvesdir, "predhazardCaries.", tnames[tt], ".PA.", "dat", sep=""))  
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")
##   for (j in 1:ncov) cat(pred.PA[[tt]]$hazard[j,], "\n", sep="   ")
##   sink()
## }


###################################################
### chunk number 30: predict08 eval=FALSE
###################################################
## for (tt in 1:4){
##   sink(paste(predCurvesdir, "quantDensCaries.", tnames[tt], ".PA.", "dat", sep=""))
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")  
##   for (j in 1:ncov){
##     for (jj in 1:nquant) cat(pred.PA[[tt]]$quant.density[[j]][jj,], "\n", sep="   ")
##   }
##   sink()
## 
##   sink(paste(predCurvesdir, "quantSurvCaries.", tnames[tt], ".PA.", "dat", sep=""))  
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")  
##   for (j in 1:ncov){
##     for (jj in 1:nquant) cat(pred.PA[[tt]]$quant.Surv[[j]][jj,], "\n", sep="   ")
##   }
##   sink()
## 
##   sink(paste(predCurvesdir, "quanthazardCaries.", tnames[tt], ".PA.", "dat", sep=""))    
##   cat(pred.PA[[tt]]$grid, "\n", sep="   ")  
##   for (j in 1:ncov){
##     for (jj in 1:nquant) cat(pred.PA[[tt]]$quant.hazard[[j]][jj,], "\n", sep="   ")
##   }
##   sink()
## }


###################################################
### chunk number 31: predict09
###################################################
tnames <- c("16", "26", "36", "46")
qnames <- c("2.5%", "50%", "97.5%")
nquant <- length(qnames)
ncov <- 32


###################################################
### chunk number 32: predict10
###################################################
pred.PA <- list()
for (tt in 1:4){
  pred.PA[[tt]] <- list()
  pred.PA[[tt]]$grid <- scan(paste(predCurvesdir, "predDensCaries.", tnames[tt], ".PA.dat", sep=""), nlines=1)
  ngrid <- length(pred.PA[[tt]]$grid)
  pred.PA[[tt]]$density <- matrix(scan(paste(predCurvesdir, "predDensCaries.", tnames[tt], ".PA.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
  pred.PA[[tt]]$Surv    <- matrix(scan(paste(predCurvesdir, "predSurvCaries.", tnames[tt], ".PA.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
  pred.PA[[tt]]$hazard  <- matrix(scan(paste(predCurvesdir, "predhazardCaries.", tnames[tt], ".PA.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
}
names(pred.PA) <- tnames


###################################################
### chunk number 33: predict11
###################################################
for (tt in 1:4){
  ngrid <- length(pred.PA[[tt]]$grid)
  pred.PA[[tt]]$quant.density <- list()
  pred.PA[[tt]]$quant.Surv    <- list()
  pred.PA[[tt]]$quant.hazard  <- list()
  for (j in 1:ncov){
    pred.PA[[tt]]$quant.density[[j]] <- matrix(scan(paste(predCurvesdir, "quantDensCaries.", tnames[tt], ".PA.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    pred.PA[[tt]]$quant.Surv[[j]]    <- matrix(scan(paste(predCurvesdir, "quantSurvCaries.", tnames[tt], ".PA.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    pred.PA[[tt]]$quant.hazard[[j]]  <- matrix(scan(paste(predCurvesdir, "quanthazardCaries.", tnames[tt], ".PA.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    rownames(pred.PA[[tt]]$quant.density[[j]]) <- qnames
    rownames(pred.PA[[tt]]$quant.Surv[[j]]) <- qnames
    rownames(pred.PA[[tt]]$quant.hazard[[j]]) <- qnames    
  }
}


###################################################
### chunk number 34: predict12
###################################################
pdfwidth <- 6
pdfheight <- 9
teeth <- c(16, 26, 36, 46)


###################################################
### chunk number 35: predict13
###################################################
pGender <- factor(c("boy", "girl"))
pBrush <- factor(c("not.daily", "daily"))
pSeal <- factor(c("no", "yes"))
pPlaque <- factor(c("none", "present"))
pPrim5 <- factor(c("sound", "dmf"))
pdata <- expand.grid(ffPrim5=pPrim5, ffPlaque=pPlaque, fSeal=pSeal, fBrush=pBrush, Gender=pGender)
pdata$Prim5 <- 1*(pdata$ffPrim5=="dmf")
pdata$Plaque <- 1*(pdata$ffPlaque=="present")
pdata$Seal <- 1*(pdata$fSeal=="yes")
pdata$Brush <- 1*(pdata$fBrush=="daily")
pdata$Girl <- 1*(pdata$Gender=="girl")
ncov <- nrow(pdata)
print(pdata)


###################################################
### chunk number 36: predict14
###################################################
sets <- list(set1=1:8, set2=9:16, set3=17:24, set4=25:32)
label <- paste(pdata$Gender, ", brush: ", pdata$fBrush, ", seal: ", pdata$fSeal, ", plaq: ", pdata$ffPlaque, ", ",
               pdata$ffPrim5, sep="")

pdf(paste(figuredir, "predSurvCI_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
for (tt in 1:4){
  for (ss in 1:length(sets)){
    par(mfrow=c(4, 2), bty="n", lwd=1.2)
    for (ii in sets[[ss]]){
      plot(pred.PA[[tt]]$grid, pred.PA[[tt]]$Surv[ii,], type="l", lty=1, xlab="Age (years)", ylab="Caries free", ylim=c(0,1))
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$quant.Surv[[ii]]["2.5%", ], lty=2)
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$quant.Surv[[ii]]["97.5%", ], lty=2)
      title(main = paste("Tooth ", teeth[tt], sep=""), sub=label[ii])
    }  
  }
}
dev.off()


###################################################
### chunk number 37: predict15
###################################################
pdf(paste(figuredir, "predhazardCI_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
ylim <- NULL
for (tt in 1:4){
  for (ss in 1:length(sets)){
    par(mfrow=c(4, 2), bty="n", lwd=1.2)
    for (ii in sets[[ss]]){
      plot(pred.PA[[tt]]$grid, pred.PA[[tt]]$quant.hazard[[ii]]["97.5%", ], type="l", lty=2,
           xlab="Age (years)", ylab="Hazard for caries", ylim=ylim) 
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$hazard[ii,], lty=1)
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$quant.hazard[[ii]]["2.5%", ], lty=2)
      title(main = paste("Tooth ", teeth[tt], sep=""), sub=label[ii])
    }  
  }
}
dev.off()


###################################################
### chunk number 38: predict16
###################################################
col <- c(rep(c("orange", "red"), 4), rep(c("darkblue", "blue"), 4))
lty <- rep(c(2, 2, 4, 4, 1, 1, 3, 3), 2)
glabel <- c("Boy", "Girl")
ncov.gender <- ncov/2

pdf(paste(figuredir, "predSurv_PA.pdf", sep=""), width=pdfheight, height=pdfwidth)
for (tt in 1:4){
  for (gg in 1:2){
    par(mfrow=c(1, 1), bty="n", lwd=1.2)
    plot(pred.PA[[tt]]$grid, pred.PA[[tt]]$Surv[(gg-1)*ncov.gender+1,],
         type="l", lty=lty[1], col=col[1], xlab="Time since emergence (years)", ylab="Caries free", ylim=c(0,1))    
    for (i in 2:ncov.gender){
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$Surv[(gg-1)*ncov.gender+i,], lty=lty[i], col=col[i])
    }
    title(main = paste("Tooth ", teeth[tt], ", ", glabel[gg], sep=""))
    legend(0, 0.7, legend=c("No plaque, sealed", "No plaque, not sealed", "Plaque, sealed", "Plaque, not sealed"),
           lty=1:4, bty="n", y.intersp=1.2, col="black")
    legend(0, 0.3, legend=c("Daily brush, sound primary", "Daily brush, dmf primary",
                            "Not daily brush, sound primary", "Not daily brush, dmf primary"),
           lty=1, bty="n", y.intersp=1.2, col=c("darkblue", "blue", "orange", "red"))
  }    
}  
dev.off()


###################################################
### chunk number 39: predict17
###################################################
pdf(paste(figuredir, "predhazard_PA.pdf", sep=""), width=pdfheight, height=pdfwidth)
ylim <- c(0, 0.25)
for (tt in 1:4){
  for (gg in 1:2){
    par(mfrow=c(1, 1), bty="n", lwd=1.2)
    plot(pred.PA[[tt]]$grid, pred.PA[[tt]]$hazard[(gg-1)*ncov.gender+1,],
         type="l", lty=lty[1], col=col[1], xlab="Time since emergence (years)", ylab="Hazard for caries", ylim=ylim)    
    for (i in 2:ncov.gender){
      lines(pred.PA[[tt]]$grid, pred.PA[[tt]]$hazard[(gg-1)*ncov.gender+i,], lty=lty[i], col=col[i])
    }
    title(main = paste("Tooth ", teeth[tt], ", ", glabel[gg], sep=""))
    legend(2.5, 0.27, legend=c("No plaque, sealed", "No plaque, not sealed", "Plaque, sealed", "Plaque, not sealed"),
           lty=1:4, bty="n", y.intersp=1.2, col="black")
    legend(2.5, 0.20, legend=c("Daily brush, sound primary", "Daily brush, dmf primary",
                            "Not daily brush, sound primary", "Not daily brush, dmf primary"),
           lty=1, bty="n", y.intersp=1.2, col=c("darkblue", "blue", "orange", "red"))    
  }    
}  
dev.off()


###################################################
### chunk number 40: errDens01
###################################################
skip <- 0
nwrite <- 5000


###################################################
### chunk number 41: errDens02 eval=FALSE
###################################################
## grid1.emerg <- seq(0, 0.8, length=50)
## grid2.emerg <- seq(0, 0.8, length=50)
## dens.emerg <- list()
## for (k in 1:2){
##   dens.emerg[[k]] <- bayesGspline(chaindir.PA[k], grid1=grid1.emerg, grid2=grid2.emerg, skip=skip, by=1, nwrite=nwrite, only.aver=TRUE)
##   sink(paste(chaindir.PA[k], "/densEmerg.dat", sep=""), append=FALSE)
##   cat(dens.emerg[[k]]$grid1, "\n", sep="  ")
##   cat(dens.emerg[[k]]$grid2, "\n", sep="  ")
##   cat(dens.emerg[[k]]$average, "\n", sep="  ")
##   sink()
## }


###################################################
### chunk number 42: errDens03 eval=FALSE
###################################################
## grid1.caries <- seq(-0.5, 8.0, length=50)
## grid2.caries <- seq(-0.5, 8.0, length=50)
## dens.caries <- list()
## for (k in 1:2){
##   dens.caries[[k]] <- bayesGspline(chaindir.PA[k], extens="_2",
##                                    grid1=grid1.caries, grid2=grid2.caries, skip=0, by=1, nwrite=100, only.aver=TRUE)
##   sink(paste(chaindir.PA[k], "/densCaries.dat", sep=""), append=FALSE)
##   cat(dens.caries[[k]]$grid1, "\n", sep="  ")
##   cat(dens.caries[[k]]$grid2, "\n", sep="  ")
##   cat(dens.caries[[k]]$average, "\n", sep="  ")
##   sink()
## }


###################################################
### chunk number 43: plotErrDens01
###################################################
pdfwidth <- 9
pdfheight <- 6
tname <- c(1646, 2636)
tname2 <- c("16, 46", "26, 36")
col <- "lightblue"


###################################################
### chunk number 44: plotErrDens02
###################################################
dens.emerg <- list()
for (k in 1:2){
  dens.emerg[[k]] <- list()
  dens.emerg[[k]]$grid1 <- scan(paste(chaindir.PA[k], "/densEmerg.dat", sep=""), nlines=1)
  dens.emerg[[k]]$grid2 <- scan(paste(chaindir.PA[k], "/densEmerg.dat", sep=""), skip=1, nlines=1)
  ngrid1 <- length(dens.emerg[[k]]$grid1)
  ngrid2 <- length(dens.emerg[[k]]$grid2)  
  if (ngrid1 != ngrid2) stop("Different grid sizes read???")
  ngrid <- ngrid1
  
  dens.emerg[[k]]$average <- matrix(scan(paste(chaindir.PA[k], "/densEmerg.dat", sep=""), skip=2, nlines=1), ncol=ngrid)
}
names(dens.emerg) <- c("right", "left")


###################################################
### chunk number 45: plotErrDens03
###################################################
dens.caries <- list()
for (k in 1:2){
  dens.caries[[k]] <- list()
  dens.caries[[k]]$grid1 <- scan(paste(chaindir.PA[k], "/densCaries.dat", sep=""), nlines=1)
  dens.caries[[k]]$grid2 <- scan(paste(chaindir.PA[k], "/densCaries.dat", sep=""), skip=1, nlines=1)
  ngrid1 <- length(dens.caries[[k]]$grid1)
  ngrid2 <- length(dens.caries[[k]]$grid2)  
  if (ngrid1 != ngrid2) stop("Different grid sizes read???")
  ngrid <- ngrid1
  
  dens.caries[[k]]$average <- matrix(scan(paste(chaindir.PA[k], "/densCaries.dat", sep=""), skip=2, nlines=1), ncol=ngrid)
}
names(dens.caries) <- c("right", "left")


###################################################
### chunk number 46: plotErrDens04
###################################################
theta <- c(30, 300)
phi <- c(30, 30)
for (k in 1:2){
  pdf(paste(figuredir, "errDensEmerg", tname[k], "_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
  par(bty="n", mfrow=c(2, 2), mar=c(4, 4, 4, 2)+0.1)
  contour(dens.emerg[[k]]$grid1, dens.emerg[[k]]$grid2, dens.emerg[[k]]$average, lty=1, 
          main=paste("Emergence, teeth ", tname2[k], sep=""),
          xlim=c(0, 1), ylim=c(0, 1))
  plot.new()
  par(mar=c(0, 0, 0, 0)+0.1)
  for (i in 1:2){
    persp(dens.emerg[[k]]$grid1, dens.emerg[[k]]$grid2, dens.emerg[[k]]$average, theta=theta[i], phi=phi[i], box=FALSE,
          xlab="y1", ylab="y2", zlab="g(y1, y2)", col=col)
  }
  dev.off()
}


###################################################
### chunk number 47: plotErrDens05
###################################################
theta <- c(30, 300)
phi <- c(30, 30)
for (k in 1:2){
  pdf(paste(figuredir, "errDensCaries", tname[k], "_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
  par(bty="n", mfrow=c(2, 2), mar=c(4, 4, 4, 2) + 0.1)
  contour(dens.caries[[k]]$grid1, dens.caries[[k]]$grid2, dens.caries[[k]]$average, lty=1, 
          main=paste("Caries, teeth ", tname2[k], sep=""),
          xlim=c(-1, 9), ylim=c(-1, 9))
  plot.new()
  par(mar=c(0, 0, 0, 0) + 0.1)
  for (i in 1:2){
    persp(dens.caries[[k]]$grid1, dens.caries[[k]]$grid2, dens.caries[[k]]$average, theta=theta[i], phi=phi[i], box=FALSE,
          xlab="y1", ylab="y2", zlab="g(y1, y2)", col=col)
  }
  dev.off()
}


###################################################
### chunk number 48: marg.errDens01
###################################################
skip <- 0
nwrite <- 1000
last.iter <- 1000


###################################################
### chunk number 49: marg.errDens02
###################################################
KK <- prior.PA.gspl.emerg$K
grid1.emerg <- seq(0, 0.8, length=100)
grid2.emerg <- seq(0, 0.8, length=100)
marg.dens.emerg <- list()
for (k in 1:2){
  marg.dens.emerg[[k]] <- marginal.bayesGspline(chaindir.PA[k], grid1=grid1.emerg, grid2=grid2.emerg, K=KK, skip=skip, by=1, nwrite=nwrite, last.iter=last.iter, only.aver=TRUE)
}


###################################################
### chunk number 50: marg.errDens03
###################################################
KK <- prior.PA.gspl.caries$K
grid1.caries <- seq(-0.5, 8.0, length=100)
grid2.caries <- seq(-0.5, 8.0, length=100)
marg.dens.caries <- list()
for (k in 1:2){
  marg.dens.caries[[k]] <- marginal.bayesGspline(chaindir.PA[k], extens="_2", grid1=grid1.caries, grid2=grid2.caries, K=KK, skip=skip, by=1, nwrite=nwrite, last.iter=last.iter, only.aver=TRUE)
}


###################################################
### chunk number 51: plot.marg.ErrDens01
###################################################
pdfwidth <- 9
pdfheight <- 6


###################################################
### chunk number 52: plot.marg.ErrDens02
###################################################
pmain <- paste("Error density, EMERGENCE, tooth ", c(16, 26, 36, 46), sep="")
pdf(paste(figuredir, "margErrDensEmerg_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(bty="n", mfrow=c(2, 2), mar=c(4, 4, 4, 2)+0.1)
plot(marg.dens.emerg[[1]]$margin1$grid, marg.dens.emerg[[1]]$margin1$average, type="l", xlab="y", ylab="g(y)", main=pmain[1])
plot(marg.dens.emerg[[2]]$margin1$grid, marg.dens.emerg[[2]]$margin1$average, type="l", xlab="y", ylab="g(y)", main=pmain[2])
plot(marg.dens.emerg[[2]]$margin2$grid, marg.dens.emerg[[2]]$margin2$average, type="l", xlab="y", ylab="g(y)", main=pmain[3])
plot(marg.dens.emerg[[1]]$margin2$grid, marg.dens.emerg[[1]]$margin2$average, type="l", xlab="y", ylab="g(y)", main=pmain[4])
dev.off()


###################################################
### chunk number 53: plot.marg.ErrDens03
###################################################
pmain <- paste("Error density, CARIES, tooth ", c(16, 26, 36, 46), sep="")
pdf(paste(figuredir, "margErrDensCaries_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(bty="n", mfrow=c(2, 2), mar=c(4, 4, 4, 2)+0.1)
plot(marg.dens.caries[[1]]$margin1$grid, marg.dens.caries[[1]]$margin1$average, type="l", xlab="y", ylab="g(y)", main=pmain[1])
plot(marg.dens.caries[[2]]$margin1$grid, marg.dens.caries[[2]]$margin1$average, type="l", xlab="y", ylab="g(y)", main=pmain[2])
plot(marg.dens.caries[[2]]$margin2$grid, marg.dens.caries[[2]]$margin2$average, type="l", xlab="y", ylab="g(y)", main=pmain[3])
plot(marg.dens.caries[[1]]$margin2$grid, marg.dens.caries[[1]]$margin2$average, type="l", xlab="y", ylab="g(y)", main=pmain[4])
dev.off()


###################################################
### chunk number 54: kendall.tau01
###################################################
skip <- 0
nwrite <- 1000
last.iter <- 1000


###################################################
### chunk number 55: kendall.tau02
###################################################
KK <- prior.PA.gspl.emerg$K
ktau.emerg <- list()
for (k in 1:2){
  ktau.emerg[[k]] <- sampled.kendall.tau(chaindir.PA[k], K=KK, skip=skip, nwrite=nwrite, last.iter=last.iter)
}
names(ktau.emerg) <- c("RIGHT", "LEFT")


###################################################
### chunk number 56: kendall.tau03
###################################################
KK <- prior.PA.gspl.caries$K
ktau.caries <- list()
for (k in 1:2){
  ktau.caries[[k]] <- sampled.kendall.tau(chaindir.PA[k], extens="_2", K=KK, skip=skip, nwrite=nwrite, last.iter=last.iter)
}
names(ktau.caries) <- c("RIGHT", "LEFT")


###################################################
### chunk number 57: kendall.tau04
###################################################
lapply(ktau.emerg, give.summary)
lapply(ktau.caries, give.summary)


###################################################
### chunk number 58: kendall.tau05
###################################################
pdfwidth <- 9
pdfheight <- 6
pdf(paste(figuredir, "histKendallTau_PA.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(bty="n", mfrow=c(2, 2), mar=c(4, 4, 4, 2)+0.1)
hist(ktau.emerg[[1]], xlab="Kendall tau", main="Emergence, teeth 16-46")
hist(ktau.emerg[[2]], xlab="Kendall tau", main="Emergence, teeth 26-36")
hist(ktau.caries[[1]], xlab="Kendall tau", main="Caries, teeth 16-46")
hist(ktau.caries[[2]], xlab="Kendall tau", main="Caries, teeth 26-36")
dev.off()


###################################################
### chunk number 59: chains01
###################################################
nbeta.emerg.PA <- 1
nbeta.caries.PA <- 5


###################################################
### chunk number 60: chains02
###################################################
itersInd.1 <- read.table(paste(chaindir.PA[1], "/iteration.sim", sep=""), header=TRUE)[,1]
itersInd.2 <- read.table(paste(chaindir.PA[2], "/iteration.sim", sep=""), header=TRUE)[,1]
if (sum(itersInd.1 - itersInd.2) > 0) stop("Different simulation lengths for the left and right chain")
itersInd <- itersInd.1
rm(list=c("itersInd.1", "itersInd.2"))
nlines <- length(itersInd)
str(itersInd)


###################################################
### chunk number 61: chains03
###################################################
mixmoment.emerg.PA <- list()
error.emerg.PA <- list()
lambda.emerg.PA <- list()
logposter.emerg.PA <- list()
beta.emerg.PA <- list()
for (k in 1:2){
  mixmoment.emerg.PA[[k]] <- matrix(scan(paste(chaindir.PA[k], "/mixmoment.sim", sep=""), skip=1, nlines=nlines), ncol=6, byrow=TRUE)
  lambda.emerg.PA[[k]]    <- matrix(scan(paste(chaindir.PA[k], "/lambda.sim", sep=""), skip=1, nlines=nlines), ncol=2, byrow=TRUE)
  logposter.emerg.PA[[k]] <- matrix(scan(paste(chaindir.PA[k], "/logposter.sim", sep=""), skip=1, nlines=nlines), ncol=4 , byrow=TRUE)
  beta.emerg.PA[[k]]      <- matrix(scan(paste(chaindir.PA[k], "/beta.sim", sep=""), skip=1, nlines=nlines), ncol=nbeta.emerg.PA, byrow=TRUE)

  colnames(mixmoment.emerg.PA[[k]]) <- scan(paste(chaindir.PA[k], "/mixmoment.sim", sep=""), what="character", nlines=1)
  colnames(lambda.emerg.PA[[k]])    <- scan(paste(chaindir.PA[k], "/lambda.sim", sep=""), what="character", nlines=1)
  colnames(logposter.emerg.PA[[k]]) <- scan(paste(chaindir.PA[k], "/logposter.sim", sep=""), what="character", nlines=1)
  colnames(beta.emerg.PA[[k]])      <- scan(paste(chaindir.PA[k], "/beta.sim", sep=""), what="character", nlines=1) 

  scale1 <- sqrt(mixmoment.emerg.PA[[k]][, "D.1.1"])
  scale2 <- sqrt(mixmoment.emerg.PA[[k]][, "D.2.2"])
  correl <- mixmoment.emerg.PA[[k]][,"D.2.1"]/(scale1*scale2)
  error.emerg.PA[[k]] <- data.frame(Mean1=mixmoment.emerg.PA[[k]][,"Mean.1"],
                                    Mean2=mixmoment.emerg.PA[[k]][,"Mean.2"],                                 
                                    Scale1=scale1,
                                    Scale2=scale2,
                                    Covar=mixmoment.emerg.PA[[k]][,"D.2.1"],
                                    Correl=correl)
  
  beta.emerg.PA[[k]] <- data.frame(Maxilla=error.emerg.PA[[k]]$Mean1 - error.emerg.PA[[k]]$Mean2,
                                   Girl=beta.emerg.PA[[k]][,"Girl"])
  logposter.emerg.PA[[k]] <- as.data.frame(logposter.emerg.PA[[k]])
  
  rm(list=c("scale1", "scale2", "correl"))
}
rm(list=c("mixmoment.emerg.PA"))
names(error.emerg.PA) <- c("right", "left")
names(beta.emerg.PA) <- c("right", "left")
names(logposter.emerg.PA) <- c("right", "left")
lambda.emerg.PA <- data.frame(lambda16=lambda.emerg.PA[[1]][,1], lambda26=lambda.emerg.PA[[2]][,1],
                              lambda36=lambda.emerg.PA[[2]][,2], lambda46=lambda.emerg.PA[[1]][,2])


###################################################
### chunk number 62: chains04
###################################################
str(error.emerg.PA$right)
str(beta.emerg.PA$right)
str(logposter.emerg.PA$right)
str(lambda.emerg.PA)


###################################################
### chunk number 63: chain05
###################################################
mixmoment.caries.PA <- list()
error.caries.PA <- list()
lambda.caries.PA <- list()
logposter.caries.PA <- list()
beta.caries.PA <- list()
for (k in 1:2){
  mixmoment.caries.PA[[k]] <- matrix(scan(paste(chaindir.PA[k], "/mixmoment_2.sim", sep=""), skip=1, nlines=nlines), ncol=6, byrow=TRUE)
  lambda.caries.PA[[k]]    <- matrix(scan(paste(chaindir.PA[k], "/lambda_2.sim", sep=""), skip=1, nlines=nlines), ncol=2, byrow=TRUE)
  logposter.caries.PA[[k]] <- matrix(scan(paste(chaindir.PA[k], "/logposter_2.sim", sep=""), skip=1, nlines=nlines), ncol=4 , byrow=TRUE)
  beta.caries.PA[[k]]      <- matrix(scan(paste(chaindir.PA[k], "/beta_2.sim", sep=""), skip=1, nlines=nlines), ncol=nbeta.caries.PA, byrow=TRUE)

  colnames(mixmoment.caries.PA[[k]]) <- scan(paste(chaindir.PA[k], "/mixmoment_2.sim", sep=""), what="character", nlines=1)
  colnames(lambda.caries.PA[[k]])    <- scan(paste(chaindir.PA[k], "/lambda_2.sim", sep=""), what="character", nlines=1)
  colnames(logposter.caries.PA[[k]]) <- scan(paste(chaindir.PA[k], "/logposter_2.sim", sep=""), what="character", nlines=1)
  colnames(beta.caries.PA[[k]])      <- scan(paste(chaindir.PA[k], "/beta_2.sim", sep=""), what="character", nlines=1) 

  scale1 <- sqrt(mixmoment.caries.PA[[k]][, "D.1.1"])
  scale2 <- sqrt(mixmoment.caries.PA[[k]][, "D.2.2"])
  correl <- mixmoment.caries.PA[[k]][,"D.2.1"]/(scale1*scale2)
  error.caries.PA[[k]] <- data.frame(Mean1=mixmoment.caries.PA[[k]][,"Mean.1"],
                                     Mean2=mixmoment.caries.PA[[k]][,"Mean.2"],                                 
                                     Scale1=scale1,
                                     Scale2=scale2,
                                     Covar=mixmoment.caries.PA[[k]][,"D.2.1"],
                                     Correl=correl)
  beta.caries.PA[[k]] <- data.frame(Maxilla=error.caries.PA[[k]]$Mean1 - error.caries.PA[[k]]$Mean2,
                                    Girl=beta.caries.PA[[k]][,"Girl"],
                                    Brush=beta.caries.PA[[k]][,"Brush"],
                                    Plaque=beta.caries.PA[[k]][,"Plaque"],
                                    Seal=beta.caries.PA[[k]][,"Seal"],
                                    Prim5=beta.caries.PA[[k]][,"Prim5"])
  logposter.caries.PA[[k]] <- as.data.frame(logposter.caries.PA[[k]])  
  
  rm(list=c("scale1", "scale2", "correl"))
}
rm(list=c("mixmoment.caries.PA"))
names(error.caries.PA) <- c("right", "left")
names(beta.caries.PA) <- c("right", "left")
names(logposter.caries.PA) <- c("right", "left")
lambda.caries.PA <- data.frame(lambda16=lambda.caries.PA[[1]][,1], lambda26=lambda.caries.PA[[2]][,1],
                               lambda36=lambda.caries.PA[[2]][,2], lambda46=lambda.caries.PA[[1]][,2])



###################################################
### chunk number 64: chains06
###################################################
str(error.caries.PA$right)
str(beta.caries.PA$right)
str(logposter.caries.PA$right)
str(lambda.caries.PA)


###################################################
### chunk number 65: summary01
###################################################
summ.error.emerg.PA <- lapply(error.emerg.PA, give.summary)
summ.error.caries.PA <- lapply(error.caries.PA, give.summary)
summ.beta.emerg.PA <- lapply(beta.emerg.PA, give.summary)
summ.beta.caries.PA <- lapply(beta.caries.PA, give.summary)
summ.lambda.emerg.PA <- give.summary(lambda.emerg.PA)
summ.lambda.caries.PA <- give.summary(lambda.caries.PA)


###################################################
### chunk number 66: summary02
###################################################
rsumm.error.emerg.PA <- lapply(summ.error.emerg.PA, round, dig=4)
rsumm.error.caries.PA <- lapply(summ.error.caries.PA, round, dig=4)
rsumm.beta.emerg.PA <- lapply(summ.beta.emerg.PA, round, dig=4)
rsumm.beta.caries.PA <- lapply(summ.beta.caries.PA, round, dig=4)
rsumm.lambda.emerg.PA <- round(summ.lambda.emerg.PA, dig=4)
rsumm.lambda.caries.PA <- round(summ.lambda.caries.PA, dig=4)


###################################################
### chunk number 67: summary03
###################################################
print(rsumm.beta.emerg.PA)
print(rsumm.error.emerg.PA)
print(rsumm.lambda.emerg.PA)


###################################################
### chunk number 68: summary04
###################################################
print(rsumm.beta.caries.PA)
print(rsumm.error.caries.PA)
print(rsumm.lambda.caries.PA)


###################################################
### chunk number 69: postDens01
###################################################
plfun <- "densplot2"
plname <- "dens"

pdfwidth <- 9
pdfheight <- 6
tt <- c(16, 26, 36, 46)


###################################################
### chunk number 70: postDens02
###################################################
pdf(paste(figuredir, plname, "_emergBeta.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.emerg.PA[[2]][, "Maxilla"]), main="Maxilla (emerg), left"))
  title(sub=paste("Maxilla = ", rsumm.beta.emerg.PA[[2]][1, "Maxilla"], 
        " (", rsumm.beta.emerg.PA[[2]][3, "Maxilla"], ", ", rsumm.beta.emerg.PA[[2]][4, "Maxilla"], ")", sep=""))
  eval(call(plfun, mcmc(beta.emerg.PA[[1]][, "Maxilla"]), main="Maxilla (emerg), right"))
  title(sub=paste("Maxilla = ", rsumm.beta.emerg.PA[[1]][1, "Maxilla"], 
        " (", rsumm.beta.emerg.PA[[1]][3, "Maxilla"], ", ", rsumm.beta.emerg.PA[[1]][4, "Maxilla"], ")", sep=""))

  eval(call(plfun, mcmc(beta.emerg.PA[[2]][, "Girl"]), main="Girl (emerg), left"))
  title(sub=paste("Girl = ", rsumm.beta.emerg.PA[[2]][1, "Girl"], 
        " (", rsumm.beta.emerg.PA[[2]][3, "Girl"], ", ", rsumm.beta.emerg.PA[[2]][4, "Girl"], ")", sep=""))
  eval(call(plfun, mcmc(beta.emerg.PA[[1]][, "Girl"]), main="Girl (emerg), right"))
  title(sub=paste("Girl = ", rsumm.beta.emerg.PA[[1]][1, "Girl"], 
        " (", rsumm.beta.emerg.PA[[1]][3, "Girl"], ", ", rsumm.beta.emerg.PA[[1]][4, "Girl"], ")", sep=""))
dev.off()


###################################################
### chunk number 71: postDens03
###################################################
pdf(paste(figuredir, plname, "_emergIntcpt.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Mean1"]), main="Intcpt (emerg), tooth 16"))
  title(sub=paste("Intercept = ", rsumm.error.emerg.PA[[1]][1, "Mean1"], 
        " (", rsumm.error.emerg.PA[[1]][3, "Mean1"], ", ", rsumm.error.emerg.PA[[1]][4, "Mean1"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Mean1"]), main="Intcpt (emerg), tooth 26"))
  title(sub=paste("Intercept = ", rsumm.error.emerg.PA[[2]][1, "Mean1"], 
        " (", rsumm.error.emerg.PA[[2]][3, "Mean1"], ", ", rsumm.error.emerg.PA[[2]][4, "Mean1"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Mean2"]), main="Intcpt (emerg), tooth 36"))
  title(sub=paste("Intercept = ", rsumm.error.emerg.PA[[2]][1, "Mean2"], 
        " (", rsumm.error.emerg.PA[[2]][3, "Mean2"], ", ", rsumm.error.emerg.PA[[2]][4, "Mean2"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Mean2"]), main="Intcpt (emerg), tooth 46"))
  title(sub=paste("Intercept = ", rsumm.error.emerg.PA[[1]][1, "Mean2"], 
        " (", rsumm.error.emerg.PA[[1]][3, "Mean2"], ", ", rsumm.error.emerg.PA[[1]][4, "Mean2"], ")", sep=""))
dev.off()


###################################################
### chunk number 72: postDens04
###################################################
pdf(paste(figuredir, plname, "_emergScale.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Scale1"]), main="Scale (emerg), tooth 16"))
  title(sub=paste("Scale = ", rsumm.error.emerg.PA[[1]][1, "Scale1"], 
    " (", rsumm.error.emerg.PA[[1]][3, "Scale1"], ", ", rsumm.error.emerg.PA[[1]][4, "Scale1"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Scale1"]), main="Scale (emerg), tooth 26"))
  title(sub=paste("Scale = ", rsumm.error.emerg.PA[[2]][1, "Scale1"], 
    " (", rsumm.error.emerg.PA[[2]][3, "Scale1"], ", ", rsumm.error.emerg.PA[[2]][4, "Scale1"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Scale2"]), main="Scale (emerg), tooth 36"))
  title(sub=paste("Scale = ", rsumm.error.emerg.PA[[2]][1, "Scale2"], 
    " (", rsumm.error.emerg.PA[[2]][3, "Scale2"], ", ", rsumm.error.emerg.PA[[2]][4, "Scale2"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Scale2"]), main="Scale (emerg), tooth 46"))
  title(sub=paste("Scale = ", rsumm.error.emerg.PA[[1]][1, "Scale2"], 
    " (", rsumm.error.emerg.PA[[1]][3, "Scale2"], ", ", rsumm.error.emerg.PA[[1]][4, "Scale2"], ")", sep=""))
dev.off()


###################################################
### chunk number 73: postDens05
###################################################
pdf(paste(figuredir, plname, "_emergCovCor.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Covar"]), main="Covar (emerg), tooth 26 vs. 36"))
  title(sub=paste("Covar = ", rsumm.error.emerg.PA[[2]][1, "Covar"], 
        " (", rsumm.error.emerg.PA[[2]][3, "Covar"], ", ", rsumm.error.emerg.PA[[2]][4, "Covar"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Covar"]), main="Covar (emerg), tooth 16 vs. 46"))
  title(sub=paste("Covar = ", rsumm.error.emerg.PA[[1]][1, "Covar"], 
        " (", rsumm.error.emerg.PA[[1]][3, "Covar"], ", ", rsumm.error.emerg.PA[[1]][4, "Covar"], ")", sep=""))

  eval(call(plfun, mcmc(error.emerg.PA[[2]][, "Correl"]), main="Correl (emerg), tooth 26 vs. 36"))
  title(sub=paste("Correl = ", rsumm.error.emerg.PA[[2]][1, "Correl"], 
        " (", rsumm.error.emerg.PA[[2]][3, "Correl"], ", ", rsumm.error.emerg.PA[[2]][4, "Correl"], ")", sep=""))
  eval(call(plfun, mcmc(error.emerg.PA[[1]][, "Correl"]), main="Correl (emerg), tooth 16 vs. 46"))
  title(sub=paste("Correl = ", rsumm.error.emerg.PA[[1]][1, "Correl"], 
        " (", rsumm.error.emerg.PA[[1]][3, "Correl"], ", ", rsumm.error.emerg.PA[[1]][4, "Correl"], ")", sep=""))
dev.off()


###################################################
### chunk number 74: postDens06
###################################################
pdf(paste(figuredir, plname, "_emergLambda.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
for (i in 1:4){
  eval(call(plfun, mcmc(lambda.emerg.PA[, i]), main=paste("Lambda (emergence), Tooth ", tt[i], sep="")))
  title(sub=paste("lambda = ", rsumm.lambda.emerg.PA[1, i], " (", rsumm.lambda.emerg.PA[3, i], ", ", rsumm.lambda.emerg.PA[4, i], ")", sep=""))
}
dev.off()


###################################################
### chunk number 75: postDens07
###################################################
pdf(paste(figuredir, plname, "_cariesBeta1.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Maxilla"]), main="Maxilla (caries), left"))
  title(sub=paste("Maxilla = ", rsumm.beta.caries.PA[[2]][1, "Maxilla"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Maxilla"], ", ", rsumm.beta.caries.PA[[2]][4, "Maxilla"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Maxilla"]), main="Maxilla (caries), right"))
  title(sub=paste("Maxilla = ", rsumm.beta.caries.PA[[1]][1, "Maxilla"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Maxilla"], ", ", rsumm.beta.caries.PA[[1]][4, "Maxilla"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Girl"]), main="Girl (caries), left"))
  title(sub=paste("Girl = ", rsumm.beta.caries.PA[[2]][1, "Girl"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Girl"], ", ", rsumm.beta.caries.PA[[2]][4, "Girl"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Girl"]), main="Girl (caries), right"))
  title(sub=paste("Girl = ", rsumm.beta.caries.PA[[1]][1, "Girl"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Girl"], ", ", rsumm.beta.caries.PA[[1]][4, "Girl"], ")", sep=""))
dev.off()


###################################################
### chunk number 76: postDens08
###################################################
pdf(paste(figuredir, plname, "_cariesBeta2.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Brush"]), main="Brush (caries), left"))
  title(sub=paste("Brush = ", rsumm.beta.caries.PA[[2]][1, "Brush"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Brush"], ", ", rsumm.beta.caries.PA[[2]][4, "Brush"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Brush"]), main="Brush (caries), right"))
  title(sub=paste("Brush = ", rsumm.beta.caries.PA[[1]][1, "Brush"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Brush"], ", ", rsumm.beta.caries.PA[[1]][4, "Brush"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Seal"]), main="Seal (caries), left"))
  title(sub=paste("Seal = ", rsumm.beta.caries.PA[[2]][1, "Seal"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Seal"], ", ", rsumm.beta.caries.PA[[2]][4, "Seal"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Seal"]), main="Seal (caries), right"))
  title(sub=paste("Seal = ", rsumm.beta.caries.PA[[1]][1, "Seal"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Seal"], ", ", rsumm.beta.caries.PA[[1]][4, "Seal"], ")", sep=""))
dev.off()


###################################################
### chunk number 77: postDens09
###################################################
pdf(paste(figuredir, plname, "_cariesBeta3.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Plaque"]), main="Plaque (caries), left"))
  title(sub=paste("Plaque = ", rsumm.beta.caries.PA[[2]][1, "Plaque"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Plaque"], ", ", rsumm.beta.caries.PA[[2]][4, "Plaque"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Plaque"]), main="Plaque (caries), right"))
  title(sub=paste("Plaque = ", rsumm.beta.caries.PA[[1]][1, "Plaque"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Plaque"], ", ", rsumm.beta.caries.PA[[1]][4, "Plaque"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.PA[[2]][, "Prim5"]), main="Prim5 (caries), left"))
  title(sub=paste("Prim5 = ", rsumm.beta.caries.PA[[2]][1, "Prim5"], 
        " (", rsumm.beta.caries.PA[[2]][3, "Prim5"], ", ", rsumm.beta.caries.PA[[2]][4, "Prim5"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.PA[[1]][, "Prim5"]), main="Prim5 (caries), right"))
  title(sub=paste("Prim5 = ", rsumm.beta.caries.PA[[1]][1, "Prim5"], 
        " (", rsumm.beta.caries.PA[[1]][3, "Prim5"], ", ", rsumm.beta.caries.PA[[1]][4, "Prim5"], ")", sep=""))
dev.off()


###################################################
### chunk number 78: postDens10
###################################################
pdf(paste(figuredir, plname, "_cariesIntcpt.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Mean1"]), main="Intcpt (caries), tooth 16"))
  title(sub=paste("Intercept = ", rsumm.error.caries.PA[[1]][1, "Mean1"], 
        " (", rsumm.error.caries.PA[[1]][3, "Mean1"], ", ", rsumm.error.caries.PA[[1]][4, "Mean1"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Mean1"]), main="Intcpt (caries), tooth 26"))
  title(sub=paste("Intercept = ", rsumm.error.caries.PA[[2]][1, "Mean1"], 
        " (", rsumm.error.caries.PA[[2]][3, "Mean1"], ", ", rsumm.error.caries.PA[[2]][4, "Mean1"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Mean2"]), main="Intcpt (caries), tooth 36"))
  title(sub=paste("Intercept = ", rsumm.error.caries.PA[[2]][1, "Mean2"], 
        " (", rsumm.error.caries.PA[[2]][3, "Mean2"], ", ", rsumm.error.caries.PA[[2]][4, "Mean2"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Mean2"]), main="Intcpt (caries), tooth 46"))
  title(sub=paste("Intercept = ", rsumm.error.caries.PA[[1]][1, "Mean2"], 
        " (", rsumm.error.caries.PA[[1]][3, "Mean2"], ", ", rsumm.error.caries.PA[[1]][4, "Mean2"], ")", sep=""))
dev.off()


###################################################
### chunk number 79: postDens11
###################################################
pdf(paste(figuredir, plname, "_cariesScale.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Scale1"]), main="Scale (caries), tooth 16"))
  title(sub=paste("Scale = ", rsumm.error.caries.PA[[1]][1, "Scale1"], 
    " (", rsumm.error.caries.PA[[1]][3, "Scale1"], ", ", rsumm.error.caries.PA[[1]][4, "Scale1"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Scale1"]), main="Scale (caries), tooth 26"))
  title(sub=paste("Scale = ", rsumm.error.caries.PA[[2]][1, "Scale1"], 
    " (", rsumm.error.caries.PA[[2]][3, "Scale1"], ", ", rsumm.error.caries.PA[[2]][4, "Scale1"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Scale2"]), main="Scale (caries), tooth 36"))
  title(sub=paste("Scale = ", rsumm.error.caries.PA[[2]][1, "Scale2"], 
    " (", rsumm.error.caries.PA[[2]][3, "Scale2"], ", ", rsumm.error.caries.PA[[2]][4, "Scale2"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Scale2"]), main="Scale (caries), tooth 46"))
  title(sub=paste("Scale = ", rsumm.error.caries.PA[[1]][1, "Scale2"], 
    " (", rsumm.error.caries.PA[[1]][3, "Scale2"], ", ", rsumm.error.caries.PA[[1]][4, "Scale2"], ")", sep=""))
dev.off()


###################################################
### chunk number 80: postDens12
###################################################
pdf(paste(figuredir, plname, "_cariesCovCor.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Covar"]), main="Covar (caries), tooth 26 vs. 36"))
  title(sub=paste("Covar = ", rsumm.error.caries.PA[[2]][1, "Covar"], 
        " (", rsumm.error.caries.PA[[2]][3, "Covar"], ", ", rsumm.error.caries.PA[[2]][4, "Covar"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Covar"]), main="Covar (caries), tooth 16 vs. 46"))
  title(sub=paste("Covar = ", rsumm.error.caries.PA[[1]][1, "Covar"], 
        " (", rsumm.error.caries.PA[[1]][3, "Covar"], ", ", rsumm.error.caries.PA[[1]][4, "Covar"], ")", sep=""))

  eval(call(plfun, mcmc(error.caries.PA[[2]][, "Correl"]), main="Correl (caries), tooth 26 vs. 36"))
  title(sub=paste("Correl = ", rsumm.error.caries.PA[[2]][1, "Correl"], 
        " (", rsumm.error.caries.PA[[2]][3, "Correl"], ", ", rsumm.error.caries.PA[[2]][4, "Correl"], ")", sep=""))
  eval(call(plfun, mcmc(error.caries.PA[[1]][, "Correl"]), main="Correl (caries), tooth 16 vs. 46"))
  title(sub=paste("Correl = ", rsumm.error.caries.PA[[1]][1, "Correl"], 
        " (", rsumm.error.caries.PA[[1]][3, "Correl"], ", ", rsumm.error.caries.PA[[1]][4, "Correl"], ")", sep=""))
dev.off()


###################################################
### chunk number 81: postDens13
###################################################
pdf(paste(figuredir, plname, "_cariesLambda.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
for (i in 1:4){
  eval(call(plfun, mcmc(lambda.caries.PA[, i]), main=paste("Lambda (caries), Tooth ", tt[i], sep="")))
  title(sub=paste("lambda = ", rsumm.lambda.caries.PA[1, i], " (", rsumm.lambda.caries.PA[3, i], ", ", rsumm.lambda.caries.PA[4, i], ")", sep=""))
}
dev.off()


###################################################
### chunk number 82: tracePlot01
###################################################
plfun <- "traceplot2"
plname <- "trace"


