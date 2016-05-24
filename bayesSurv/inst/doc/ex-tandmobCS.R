###################################################
### chunk number 1: directories
###################################################
anadir <- "/home/komari/win/work/papers/gsplineTand/anaDMFsymm/"
chaindir.CS <- paste(anadir, "chains/modelCS", sep="")
initdir.CS <- paste(anadir, "chains/end_CS", sep="")
plotdir.CS <- paste(anadir, "summPlots/modelCS/", sep="")
resultdir <- paste(anadir, "results/", sep="")
predCurvesdir <- paste(resultdir, "predCurves/", sep="")
figuredir <- "/home/komari/win/work/papers/gsplineTand/RforCRAN/figuresCS/"


###################################################
### chunk number 2: loadLibData
###################################################
library(bayesSurv)
data(tandmobRoos)


###################################################
### chunk number 3: data01
###################################################
teeth <- c(16, 26, 36, 46)
maxil <- c(1, 1, 0, 0)
dmfs <- list()
for (ti in 1:length(teeth)){
  tt <- teeth[ti]
  remove <- is.na(tandmobRoos[, paste("FBEG.", tt, sep="")])  |
            (tandmobRoos[, paste("TOOTH.", tt, sep="")] == 0)  |
            is.na(tandmobRoos[, paste("T", tt+39, "d", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt+39, "m", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt+39, "f", sep="")])  |
            is.na(tandmobRoos[, paste("T", tt+39, "s", sep="")])  |
            is.na(tandmobRoos[, paste("SEAL.", tt, sep="")])  |
            is.na(tandmobRoos[, "FREQ.BR"])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt, ".1", sep="")])  |
            is.na(tandmobRoos[, paste("PLAQUE.", tt, ".2", sep="")])
  temp <- tandmobRoos[!remove,]  
  temp <- data.frame(Idnr=temp$IDNR,
                     Tooth=rep(tt, dim(temp)[1]), 
                     Ebeg=temp[, paste("EBEG.", tt, sep="")],
                     Eend=temp[, paste("EEND.", tt, sep="")],
                     Fbeg=temp[, paste("FBEG.", tt, sep="")],
                     Fend=temp[, paste("FEND.", tt, sep="")],
                     Maxilla=rep(maxil[ti], dim(temp)[1]),
                     Gender=temp[, "GENDER"],
                     Girl=temp[, "GIRL"], 
                     Brush=temp[, "FREQ.BR"],
                     Plaque.1=temp[, paste("PLAQUE.", tt, ".1", sep="")],
                     Plaque.2=temp[, paste("PLAQUE.", tt, ".2", sep="")],
                     Seal=temp[, paste("SEAL.", tt, sep="")],
                     Prim5d=temp[, paste("T", tt+39, "d", sep="")],
                     Prim5m=temp[, paste("T", tt+39, "m", sep="")],
                     Prim5f=temp[, paste("T", tt+39, "f", sep="")],
                     Prim4d=temp[, paste("T", tt+38, "d", sep="")],
                     Prim4m=temp[, paste("T", tt+38, "m", sep="")],
                     Prim4f=temp[, paste("T", tt+38, "f", sep="")])
  temp$fBrush <- factor(temp$Brush, levels=0:1, labels=c("not.daily", "daily"))
  temp$fPlaque <- ordered(0*(temp$Plaque.1==0 & temp$Plaque.2==0) + 1*(temp$Plaque.1==1) + 2*(temp$Plaque.2==1),
                         levels=0:2, labels=c("none", "pits.fiss", "total"))
  temp$fSeal <- factor(temp$Seal, levels=0:1, labels=c("no", "yes"))
  temp$fPrim5 <- factor(0*(temp$Prim5d==0 & temp$Prim5m==0 & temp$Prim5f==0) + 
                        1*(temp$Prim5d==1) + 2*(temp$Prim5f==1) + 3*(temp$Prim5m==1),
                        levels=0:3, labels=c("sound", "decayed", "filled", "missing"))
  temp$fPrim4 <- factor(0*(temp$Prim4d==0 & temp$Prim4m==0 & temp$Prim4f==0) + 
                        1*(temp$Prim4d==1) + 2*(temp$Prim4m==1) + 3*(temp$Prim4f==1), 
                        levels=0:3, labels=c("sound", "decayed", "missing", "filled"))
  
  temp$Plaque <- 0*(temp$fPlaque=="none") + 1*(temp$fPlaque=="pits.fiss" | temp$fPlaque=="total")  
  temp$Prim5 <- 0*(temp$fPrim5=="sound") + 1*(temp$fPrim5=="decayed" | temp$fPrim5=="filled" | temp$fPrim5=="missing")

  dmfs[[ti]] <- temp
  rm(list="temp") 
}
names(dmfs) <- teeth
n.sample <- sapply(dmfs, nrow)
print(n.sample)


###################################################
### chunk number 4: data02
###################################################
time.0 <- 5


###################################################
### chunk number 5: data03
###################################################
data.CS <- rbind(dmfs$"16", dmfs$"26", dmfs$"36", dmfs$"46")
data.CS <- data.CS[order(data.CS$Idnr), ]
data.CS$Tooth <- factor(data.CS$Tooth)
rownames(data.CS) <- 1:dim(data.CS)[1]

data.CS$Ebeg <- data.CS$Ebeg - time.0
data.CS$Ebeg[data.CS$Ebeg <= 0] <- NA
data.CS$Eend <- data.CS$Eend - time.0
data.CS$Fbeg <- data.CS$Fbeg - time.0
data.CS$Fbeg[data.CS$Fbeg <= 0] <- NA
data.CS$Fend <- data.CS$Fend - time.0

rm(list=c("remove", "ti", "tt", "teeth"))


###################################################
### chunk number 6: data04
###################################################
print(data.CS[1:10,])


###################################################
### chunk number 7: initAFT01
###################################################
emerg.CS.aft <- survreg(Surv(Ebeg, Eend, type="interval2")~Tooth+Girl, dist="loglogistic", data=data.CS)


###################################################
### chunk number 8: initAFT02
###################################################
## To get some impression on the distribution of caries we ad hod impute emergence times as 
## either midpoints of intervals or a~midpoint between 6 and left-censored time or equal to the right-censored time
left <- is.na(data.CS$Ebeg)
interv <- !is.na(data.CS$Ebeg) & !is.na(data.CS$Eend)
right <- is.na(data.CS$Eend)

onset <- data.CS$Ebeg
onset[left] <- 0.5*(0 + data.CS$Eend[left])
onset[right] <- data.CS$Ebeg[right]
onset[interv] <- 0.5*(data.CS$Ebeg[interv] + data.CS$Eend[interv])
vfbeg2 <- data.CS$Fbeg - onset
vfbeg2[vfbeg2 < 0] <-  NA            ### -> both emergence and caries were in one interval, make it left censored now
vfend2 <- data.CS$Fend - onset

caries.CS.aft <- survreg(Surv(vfbeg2, vfend2, type="interval2")~Tooth+Girl+Brush+Plaque+Seal+Prim5,
                         dist="loglogistic", data=data.CS)


###################################################
### chunk number 9: initAFT03
###################################################
summary(emerg.CS.aft)


###################################################
### chunk number 10: initAFT04
###################################################
summary(caries.CS.aft)


###################################################
### chunk number 11: prior01
###################################################
prior.CS.gspl.emerg <- list(
  specification   = 2, 
  K               = 15,  
  izero           = 0,
  neighbor.system = "uniCAR",
  order           = 3,
  equal.lambda    = TRUE,
  c4delta         = 1.5,
  prior.lambda    = "gamma",
  prior.intercept = "normal",
  prior.scale     = "gamma",
  prior.gamma     = "fixed",
  prior.sigma     = "fixed",
  shape.lambda    = 1,
  rate.lambda     = 0.005,
  mean.intercept  = 0,
  var.intercept   = 100,
  shape.scale     = 1,
  rate.scale      = 0.005
)

prior.CS.gspl.caries <- prior.CS.gspl.emerg
prior.CS.b.emerg <- prior.CS.gspl.emerg
prior.CS.b.caries <- prior.CS.gspl.caries


###################################################
### chunk number 12: prior02
###################################################
prior.CS.beta.emerg <- list(mean.prior=rep(0, 4), var.prior=rep(1e2, 4))
prior.CS.beta.caries <- list(mean.prior=rep(0, 8), var.prior=rep(1e2, 8))


###################################################
### chunk number 13: init01
###################################################
init.CS.emerg <- list(
  lambda    = 3000,
  intercept = 0.40,
  scale     = 0.10,
  gamma     = 0,
  sigma     = 0.2,
  beta      = c(-0.01, 0.01, 0.02, -0.02),
  lambda.b    = 1000,
  intercept.b = 0.00,
  scale.b     = 0.05,
  gamma.b     = 0.00,
  sigma.b     = 0.2
)


###################################################
### chunk number 14: init02
###################################################
init.CS.caries <- list(
  lambda    = 3000,
  intercept = 2.16,
  scale     = 0.55,
  gamma     = 0,
  sigma     = 0.2,
  beta      = c(0.03, 0, 0, -0.06, 0.35, -0.24, 0.07, -0.63),
  lambda.b    = 1000,
  intercept.b = 0.00,
  scale.b     = 0.05,
  gamma.b     = 0.00,
  sigma.b     = 0.2
)


###################################################
### chunk number 15: init03
###################################################
init.CS.emerg <- list()
  init.CS.emerg$iter      <- scan(paste(initdir.CS, "/iteration.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$lambda    <- scan(paste(initdir.CS, "/lambda.sim", sep=""), skip=1, nlines=1)
  init.CS.gspline         <- scan(paste(initdir.CS, "/gspline.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$gamma     <- init.CS.gspline[1]
  init.CS.emerg$sigma     <- init.CS.gspline[2]
  init.CS.emerg$intercept <- init.CS.gspline[4]
  init.CS.emerg$scale     <- init.CS.gspline[5]
  init.CS.emerg$beta      <- scan(paste(initdir.CS, "/beta.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$a         <- scan(paste(initdir.CS, "/mlogweight.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$y         <- scan(paste(initdir.CS, "/Y.sim", sep=""), skip=0, nlines=1)
  init.CS.r <- scan(paste(initdir.CS, "/r.sim", sep=""), skip=0, nlines=1)
  init.CS.emerg$r <- vecr2matr(init.CS.r, prior.CS.gspl.emerg$K)

  init.CS.emerg$lambda.b  <- scan(paste(initdir.CS, "/lambda_b.sim", sep=""), skip=1, nlines=1)
  init.CS.gspline.b       <- scan(paste(initdir.CS, "/gspline_b.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$gamma.b     <- init.CS.gspline.b[1]
  init.CS.emerg$sigma.b     <- init.CS.gspline.b[2]
  init.CS.emerg$intercept.b <- init.CS.gspline.b[4]
  init.CS.emerg$scale.b     <- init.CS.gspline.b[5]
  init.CS.emerg$a.b         <- scan(paste(initdir.CS, "/mlogweight_b.sim", sep=""), skip=1, nlines=1)
  init.CS.emerg$b <- scan(paste(initdir.CS, "/b.sim", sep=""), skip=0, nlines=1)
  init.CS.r.b <- scan(paste(initdir.CS, "/r_b.sim", sep=""), skip=0, nlines=1)
  init.CS.emerg$r.b <- vecr2matr(init.CS.r.b, prior.CS.b.emerg$K)

init.CS.caries <- list()
  init.CS.caries$iter      <- scan(paste(initdir.CS, "/iteration.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$lambda    <- scan(paste(initdir.CS, "/lambda_2.sim", sep=""), skip=1, nlines=1)
  init.CS.gspline         <- scan(paste(initdir.CS, "/gspline_2.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$gamma     <- init.CS.gspline[1]
  init.CS.caries$sigma     <- init.CS.gspline[2]
  init.CS.caries$intercept <- init.CS.gspline[4]
  init.CS.caries$scale     <- init.CS.gspline[5]
  init.CS.caries$beta      <- scan(paste(initdir.CS, "/beta_2.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$a         <- scan(paste(initdir.CS, "/mlogweight_2.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$y         <- scan(paste(initdir.CS, "/Y_2.sim", sep=""), skip=0, nlines=1)
  init.CS.r <- scan(paste(initdir.CS, "/r_2.sim", sep=""), skip=0, nlines=1)
  init.CS.caries$r <- vecr2matr(init.CS.r, prior.CS.gspl.caries$K)

  init.CS.caries$lambda.b <- scan(paste(initdir.CS, "/lambda_b2.sim", sep=""), skip=1, nlines=1)
  init.CS.gspline.b       <- scan(paste(initdir.CS, "/gspline_b2.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$gamma.b     <- init.CS.gspline.b[1]
  init.CS.caries$sigma.b     <- init.CS.gspline.b[2]
  init.CS.caries$intercept.b <- init.CS.gspline.b[4]
  init.CS.caries$scale.b     <- init.CS.gspline.b[5]
  init.CS.caries$a.b         <- scan(paste(initdir.CS, "/mlogweight_b2.sim", sep=""), skip=1, nlines=1)
  init.CS.caries$b <- scan(paste(initdir.CS, "/b_2.sim", sep=""), skip=0, nlines=1)
  init.CS.r.b <- scan(paste(initdir.CS, "/r_b2.sim", sep=""), skip=0, nlines=1)
  init.CS.caries$r.b <- vecr2matr(init.CS.r.b, prior.CS.b.caries$K)


###################################################
### chunk number 16: init04
###################################################
str(init.CS.emerg)


###################################################
### chunk number 17: init05
###################################################
str(init.CS.caries)


###################################################
### chunk number 18: sample01
###################################################
nsimul.CS <- list(niter=1100000, nthin=3, nburn=1000000, nwrite=250)


###################################################
### chunk number 19: sample02 eval=FALSE
###################################################
sample.CS <- bayessurvreg3(
    formula=Surv(Ebeg, Eend, type="interval2")~ Tooth + Girl + cluster(Idnr),
    random=~1,
    formula2=Surv(Fbeg, Fend, type="interval2")~ Tooth + Girl + Brush + Plaque + Seal + Prim5 + cluster(Idnr),
    random2=~1,
    onlyX=FALSE,
    dir=chaindir.CS, nsimul=nsimul.CS,
    prior=prior.CS.gspl.emerg, prior2=prior.CS.gspl.caries,
    prior.beta=prior.CS.beta.emerg, prior.beta2=prior.CS.beta.caries,
    prior.b=prior.CS.b.emerg, prior.b2=prior.CS.b.caries,
    init=init.CS.emerg, init2=init.CS.caries,
    store=list(a=TRUE, a2=TRUE, a.b=TRUE, a.b2=TRUE),
    data=data.CS)


###################################################
### chunk number 20: predict01
###################################################
skip <- 0
nwrite <- 5000
pred.grid <- c(seq(0.1, 1.5, by=0.05), seq(1.65, 6, by=0.15))
quants <- c(0.025, 0.5, 0.975)
nquant <- length(quants)


###################################################
### chunk number 21: predict02
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

pred.data <- data.frame(Idnr=1:(4*ncov), Tooth=factor(c(rep(16, ncov), rep(26, ncov), rep(36, ncov), rep(46, ncov))))
pred.data$Tooth26 <- 1*(pred.data$Tooth == 26)
pred.data$Tooth36 <- 1*(pred.data$Tooth == 36)
pred.data$Tooth46 <- 1*(pred.data$Tooth == 46)
pred.data <- cbind(pred.data, rbind(pdata, pdata, pdata, pdata))


###################################################
### chunk number 22: predict03
###################################################
print(pred.data[1:ncov,])


###################################################
### chunk number 23: predict04
###################################################
tooth <- c(16, 26, 36, 46)
start <- list(tooth16=seq(1, ncov, by=2), tooth26=seq(ncov+1, 2*ncov, by=2),
              tooth36=seq(2*ncov+1, 3*ncov, by=2), tooth46=seq(3*ncov+1, 4*ncov, by=2))
end <- list(tooth16=seq(2, ncov, by=2), tooth26=seq(ncov+2, 2*ncov, by=2),
            tooth36=seq(2*ncov+2, 3*ncov, by=2), tooth46=seq(3*ncov+2, 4*ncov, by=2))


###################################################
### chunk number 24: predict05 eval=FALSE
###################################################
pred.CS <- list()
for (k in 1:4){    ## loop over teeth
  cat("Performing TOOTH ", tooth[k], "\n", sep="")
  pred.CS[[k]] <- list()
  
  for (ii in 1:length(start[[k]])){    ## loop over sets of covariates      
      cat("Performing the covariate set number ", ii, "\n", sep="")
      pdata.now <- pred.data[start[[k]][ii]:end[[k]][ii], ]
      nr <- nrow(pdata.now)
      pred.CS[[k]][[ii]] <-
          predictive2(Surv(rep(1, nr), rep(1, nr)) ~ Tooth26+Tooth36+Tooth46+Girl+Brush+Plaque+Seal+Prim5 + cluster(Idnr), random=~1,
                             grid=pred.grid, data=pdata.now,
                             Gspline=list(dim=1, K=15), quantile=quants,
                             skip=skip, by=1, nwrite=nwrite, only.aver=FALSE,
                             predict=list(density=TRUE, Surv=TRUE, hazard=TRUE),
                             dir=chaindir.CS, extens="_2", extens.random="_b2", version=3)
  }
}  
names(pred.CS) <- paste(tooth)


###################################################
### chunk number 25: predict06 eval=FALSE
###################################################
predd.CS <- list()
for (tt in 1:4){    ## loop over teeth
  predd.CS[[tt]] <- list()
  predd.CS[[tt]]$grid <- pred.CS[[tt]][[1]]$grid
  
  predd.CS[[tt]]$Surv <- pred.CS[[tt]][[1]]$Surv
  predd.CS[[tt]]$hazard <- pred.CS[[tt]][[1]]$hazard
  predd.CS[[tt]]$density <- pred.CS[[tt]][[1]]$density

  predd.CS[[tt]]$quant.Surv <- pred.CS[[tt]][[1]]$quant.Surv
  predd.CS[[tt]]$quant.hazard <- pred.CS[[tt]][[1]]$quant.hazard
  predd.CS[[tt]]$quant.density <- pred.CS[[tt]][[1]]$quant.density
  
  for (ii in 2:length(start[[tt]])){
    predd.CS[[tt]]$Surv <- rbind(predd.CS[[tt]]$Surv, pred.CS[[tt]][[ii]]$Surv)
    predd.CS[[tt]]$hazard <- rbind(predd.CS[[tt]]$hazard, pred.CS[[tt]][[ii]]$hazard)    
    predd.CS[[tt]]$density <- rbind(predd.CS[[tt]]$density, pred.CS[[tt]][[ii]]$density)

    templ <- length(predd.CS[[tt]]$quant.Surv)
    for (ij in 1:length(pred.CS[[tt]][[ii]]$quant.Surv)){    
      predd.CS[[tt]]$quant.Surv[[templ + ij]] <- pred.CS[[tt]][[ii]]$quant.Surv[[ij]]
      predd.CS[[tt]]$quant.hazard[[templ + ij]] <- pred.CS[[tt]][[ii]]$quant.hazard[[ij]]
      predd.CS[[tt]]$quant.density[[templ + ij]] <- pred.CS[[tt]][[ii]]$quant.density[[ij]]      
    }
    
  }  
}  
names(predd.CS) <- paste(tooth)
pred.CS <- predd.CS
rm(list="predd.CS")


###################################################
### chunk number 26: predict07 eval=FALSE
###################################################
for (tt in 1:4){
  sink(paste(predCurvesdir, "predDensCaries.", tooth[tt], ".CS.", "dat", sep=""))
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")
  for (j in 1:ncov) cat(pred.CS[[tt]]$density[j,], "\n", sep="   ")
  sink()

  sink(paste(predCurvesdir, "predSurvCaries.", tooth[tt], ".CS.", "dat", sep=""))
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")
  for (j in 1:ncov) cat(pred.CS[[tt]]$Surv[j,], "\n", sep="   ")
  sink()

  sink(paste(predCurvesdir, "predhazardCaries.", tooth[tt], ".CS.", "dat", sep=""))  
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")
  for (j in 1:ncov) cat(pred.CS[[tt]]$hazard[j,], "\n", sep="   ")
  sink()
}


###################################################
### chunk number 27: predict08 eval=FALSE
###################################################
for (tt in 1:4){
  sink(paste(predCurvesdir, "quantDensCaries.", tooth[tt], ".CS.", "dat", sep=""))
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")  
  for (j in 1:ncov){
    for (jj in 1:nquant) cat(pred.CS[[tt]]$quant.density[[j]][jj,], "\n", sep="   ")
  }
  sink()

  sink(paste(predCurvesdir, "quantSurvCaries.", tooth[tt], ".CS.", "dat", sep=""))  
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")  
  for (j in 1:ncov){
    for (jj in 1:nquant) cat(pred.CS[[tt]]$quant.Surv[[j]][jj,], "\n", sep="   ")
  }
  sink()

  sink(paste(predCurvesdir, "quanthazardCaries.", tooth[tt], ".CS.", "dat", sep=""))    
  cat(pred.CS[[tt]]$grid, "\n", sep="   ")  
  for (j in 1:ncov){
    for (jj in 1:nquant) cat(pred.CS[[tt]]$quant.hazard[[j]][jj,], "\n", sep="   ")
  }
  sink()
}


###################################################
### chunk number 28: predict09
###################################################
tnames <- c("16", "26", "36", "46")
qnames <- c("2.5%", "50%", "97.5%")
nquant <- length(qnames)
ncov <- 32


###################################################
### chunk number 29: predict10
###################################################
pred.CS <- list()
for (tt in 1:4){
  pred.CS[[tt]] <- list()
  pred.CS[[tt]]$grid <- scan(paste(predCurvesdir, "predDensCaries.", tnames[tt], ".CS.dat", sep=""), nlines=1)
  ngrid <- length(pred.CS[[tt]]$grid)
  pred.CS[[tt]]$density <- matrix(scan(paste(predCurvesdir, "predDensCaries.", tnames[tt], ".CS.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
  pred.CS[[tt]]$Surv    <- matrix(scan(paste(predCurvesdir, "predSurvCaries.", tnames[tt], ".CS.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
  pred.CS[[tt]]$hazard  <- matrix(scan(paste(predCurvesdir, "predhazardCaries.", tnames[tt], ".CS.dat", sep=""), skip=1),
                                  ncol=ngrid, byrow=TRUE)
}
names(pred.CS) <- tnames


###################################################
### chunk number 30: predict11
###################################################
for (tt in 1:4){
  ngrid <- length(pred.CS[[tt]]$grid)
  pred.CS[[tt]]$quant.density <- list()
  pred.CS[[tt]]$quant.Surv    <- list()
  pred.CS[[tt]]$quant.hazard  <- list()
  for (j in 1:ncov){
    pred.CS[[tt]]$quant.density[[j]] <- matrix(scan(paste(predCurvesdir, "quantDensCaries.", tnames[tt], ".CS.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    pred.CS[[tt]]$quant.Surv[[j]]    <- matrix(scan(paste(predCurvesdir, "quantSurvCaries.", tnames[tt], ".CS.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    pred.CS[[tt]]$quant.hazard[[j]]  <- matrix(scan(paste(predCurvesdir, "quanthazardCaries.", tnames[tt], ".CS.dat", sep=""),
                                                    skip=nquant*(j-1)+1, nlines=nquant), ncol=ngrid, byrow=TRUE)
    rownames(pred.CS[[tt]]$quant.density[[j]]) <- qnames
    rownames(pred.CS[[tt]]$quant.Surv[[j]]) <- qnames
    rownames(pred.CS[[tt]]$quant.hazard[[j]]) <- qnames    
  }
}


###################################################
### chunk number 31: predict12
###################################################
pdfwidth <- 6
pdfheight <- 9
teeth <- c(16, 26, 36, 46)


###################################################
### chunk number 32: predict13
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
### chunk number 33: predict14
###################################################
sets <- list(set1=1:8, set2=9:16, set3=17:24, set4=25:32)
label <- paste(pdata$Gender, ", brush: ", pdata$fBrush, ", seal: ", pdata$fSeal, ", plaq: ", pdata$ffPlaque, ", ",
               pdata$ffPrim5, sep="")

pdf(paste(figuredir, "predSurvCI_CS.pdf", sep=""), width=pdfwidth, height=pdfheight)
for (tt in 1:4){
  for (ss in 1:length(sets)){
    par(mfrow=c(4, 2), bty="n", lwd=1.2)
    for (ii in sets[[ss]]){
      plot(pred.CS[[tt]]$grid, pred.CS[[tt]]$Surv[ii,], type="l", lty=1, xlab="Time since emergence (years)", ylab="Caries free", ylim=c(0,1))
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$quant.Surv[[ii]]["2.5%", ], lty=2)
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$quant.Surv[[ii]]["97.5%", ], lty=2)
      title(main = paste("Tooth ", teeth[tt], sep=""), sub=label[ii])
    }  
  }
}
dev.off()


###################################################
### chunk number 34: predict15
###################################################
pdf(paste(figuredir, "predhazardCI_CS.pdf", sep=""), width=pdfwidth, height=pdfheight)
ylim <- NULL
for (tt in 1:4){
  for (ss in 1:length(sets)){
    par(mfrow=c(4, 2), bty="n", lwd=1.2)
    for (ii in sets[[ss]]){
      plot(pred.CS[[tt]]$grid, pred.CS[[tt]]$quant.hazard[[ii]]["97.5%", ], type="l", lty=2,
           xlab="Time since emergence (years)", ylab="Hazard for caries", ylim=ylim) 
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$hazard[ii,], lty=1)
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$quant.hazard[[ii]]["2.5%", ], lty=2)
      title(main = paste("Tooth ", teeth[tt], sep=""), sub=label[ii])
    }  
  }
}
dev.off()


###################################################
### chunk number 35: predict16
###################################################
col <- c(rep(c("orange", "red"), 4), rep(c("darkblue", "blue"), 4))
lty <- rep(c(2, 2, 4, 4, 1, 1, 3, 3), 2)
glabel <- c("Boy", "Girl")
ncov.gender <- ncov/2

pdf(paste(figuredir, "predSurv_CS.pdf", sep=""), width=pdfheight, height=pdfwidth)
for (tt in 1:4){
  for (gg in 1:2){
    par(mfrow=c(1, 1), bty="n", lwd=1.2)
    plot(pred.CS[[tt]]$grid, pred.CS[[tt]]$Surv[(gg-1)*ncov.gender+1,],
         type="l", lty=lty[1], col=col[1], xlab="Time since emergence (years)", ylab="Caries free", ylim=c(0,1))    
    for (i in 2:ncov.gender){
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$Surv[(gg-1)*ncov.gender+i,], lty=lty[i], col=col[i])
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
### chunk number 36: predict17
###################################################
pdf(paste(figuredir, "predhazard_CS.pdf", sep=""), width=pdfheight, height=pdfwidth)
ylim <- c(0, 0.4)
for (tt in 1:4){
  for (gg in 1:2){
    par(mfrow=c(1, 1), bty="n", lwd=1.2)
    plot(pred.CS[[tt]]$grid, pred.CS[[tt]]$hazard[(gg-1)*ncov.gender+1,],
         type="l", lty=lty[1], col=col[1], xlab="Time since emergence (years)", ylab="Hazard for caries", ylim=ylim)    
    for (i in 2:ncov.gender){
      lines(pred.CS[[tt]]$grid, pred.CS[[tt]]$hazard[(gg-1)*ncov.gender+i,], lty=lty[i], col=col[i])
    }
    title(main = paste("Tooth ", teeth[tt], ", ", glabel[gg], sep=""))
    legend(0, 0.4, legend=c("No plaque, sealed", "No plaque, not sealed", "Plaque, sealed", "Plaque, not sealed"),
           lty=1:4, bty="n", y.intersp=1.2, col="black")
    legend(0, 0.3, legend=c("Daily brush, sound primary", "Daily brush, dmf primary",
                            "Not daily brush, sound primary", "Not daily brush, dmf primary"),
           lty=1, bty="n", y.intersp=1.2, col=c("darkblue", "blue", "orange", "red"))
  }    
}  
dev.off()


###################################################
### chunk number 37: errDens01
###################################################
skip <- 0
nwrite <- 5000


###################################################
### chunk number 38: errDens02 eval=FALSE
###################################################
grid.emerg <- seq(0.35, 0.55, length=100)
dens.emerg <- bayesGspline(chaindir.CS, grid1=grid.emerg, skip=skip, by=1, nwrite=nwrite,
                           extens="", extens.adjust="_b", only.aver=TRUE, version=30)

sink(paste(chaindir.CS, "/densEmerg.dat", sep=""), append=FALSE)
cat(dens.emerg$grid, "\n", sep="  ")
cat(dens.emerg$average, "\n", sep="  ")
sink()


###################################################
### chunk number 39: errDens03 eval=FALSE
###################################################
grid.random.emerg <- seq(-0.5, 0.5, length=100)
dens.random.emerg <- bayesGspline(chaindir.CS, grid1=grid.random.emerg, skip=skip, by=1, nwrite=nwrite, 
                                  extens="_b", extens.adjust="_b", only.aver=TRUE, version=31)

sink(paste(chaindir.CS, "/densRandEmerg.dat", sep=""), append=FALSE)
cat(dens.random.emerg$grid, "\n", sep="  ")
cat(dens.random.emerg$average, "\n", sep="  ")
sink()


###################################################
### chunk number 40: errDens04 eval=FALSE
###################################################
grid.caries <- seq(-1, 4.0, length=100)
dens.caries <- bayesGspline(chaindir.CS, grid1=grid.caries, skip=skip, by=1, nwrite=nwrite, 
                            extens="_2", extens.adjust="_b2", only.aver=TRUE, version=30)

sink(paste(chaindir.CS, "/densCaries.dat", sep=""), append=FALSE)
cat(dens.caries$grid, "\n", sep="  ")
cat(dens.caries$average, "\n", sep="  ")
sink()


###################################################
### chunk number 41: errDens05 eval=FALSE
###################################################
grid.random.caries <- seq(-2.0, 1.5, length=100)
dens.random.caries <- bayesGspline(chaindir.CS, grid1=grid.random.caries, skip=skip, by=1, nwrite=nwrite, 
                                   extens="_b2", extens.adjust="_b2", only.aver=TRUE, version=31)

sink(paste(chaindir.CS, "/densRandCaries.dat", sep=""), append=FALSE)
cat(dens.random.caries$grid, "\n", sep="  ")
cat(dens.random.caries$average, "\n", sep="  ")
sink()


###################################################
### chunk number 42: plotErrDens01
###################################################
pdfwidth <- 9
pdfheight <- 6


###################################################
### chunk number 43: plotErrDens02
###################################################
dens.emerg <- list(grid=scan(paste(chaindir.CS, "/densEmerg.dat", sep=""), nlines=1), 
                   average=scan(paste(chaindir.CS, "/densEmerg.dat", sep=""), skip=1, nlines=1))
dens.random.emerg <- list(grid=scan(paste(chaindir.CS, "/densRandEmerg.dat", sep=""), nlines=1), 
                          average=scan(paste(chaindir.CS, "/densRandEmerg.dat", sep=""), skip=1, nlines=1))
dens.caries <- list(grid=scan(paste(chaindir.CS, "/densCaries.dat", sep=""), nlines=1), 
                    average=scan(paste(chaindir.CS, "/densCaries.dat", sep=""), skip=1, nlines=1))
dens.random.caries <- list(grid=scan(paste(chaindir.CS, "/densRandCaries.dat", sep=""), nlines=1), 
                           average=scan(paste(chaindir.CS, "/densRandCaries.dat", sep=""), skip=1, nlines=1))


###################################################
### chunk number 44: plotErrDens03
###################################################
pdf(paste(figuredir, "err_randomDens_CS.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfcol=c(2,2), bty="n")
plot(dens.emerg$grid, dens.emerg$average, type="l", main="Error, emergence", xlab=expression(zeta), ylab=expression(g(zeta)))
plot(dens.random.emerg$grid, dens.random.emerg$average, type="l", main="Random intcpt, emergence", xlab="d", ylab="g(d)")
plot(dens.caries$grid, dens.caries$average, type="l", main="Error, caries", xlab=expression(epsilon), ylab=expression(g(epsilon)))
plot(dens.random.caries$grid, dens.random.caries$average, type="l", main="Random intcpt, caries", xlab="b", ylab="g(b)")
dev.off()


###################################################
### chunk number 45: chains01
###################################################
nbeta.emerg.CS <- 4
nbeta.caries.CS <- 8


###################################################
### chunk number 46: chains02
###################################################
itersInd <- read.table(paste(chaindir.CS, "/iteration.sim", sep=""), header=TRUE)[,1]
nlines <- length(itersInd)
str(itersInd)


###################################################
### chunk number 47: chains03
###################################################
mixmoment.emerg.CS <- matrix(scan(paste(chaindir.CS, "/mixmoment.sim", sep=""), skip=1, nlines=nlines), ncol=3, byrow=TRUE)
lambda.emerg.CS    <- matrix(scan(paste(chaindir.CS, "/lambda.sim", sep=""), skip=1, nlines=nlines), ncol=1, byrow=TRUE)
logposter.emerg.CS <- matrix(scan(paste(chaindir.CS, "/logposter.sim", sep=""), skip=1, nlines=nlines), ncol=3 , byrow=TRUE)
beta.emerg.CS      <- matrix(scan(paste(chaindir.CS, "/beta.sim", sep=""), skip=1, nlines=nlines), ncol=nbeta.emerg.CS, byrow=TRUE)
mixmoment.b.emerg.CS <- matrix(scan(paste(chaindir.CS, "/mixmoment_b.sim", sep=""), skip=1, nlines=nlines), ncol=3, byrow=TRUE)
lambda.b.emerg.CS    <- matrix(scan(paste(chaindir.CS, "/lambda_b.sim", sep=""), skip=1, nlines=nlines), ncol=1, byrow=TRUE)
logposter.b.emerg.CS <- matrix(scan(paste(chaindir.CS, "/logposter_b.sim", sep=""), skip=1, nlines=nlines), ncol=3 , byrow=TRUE)

colnames(mixmoment.emerg.CS) <- scan(paste(chaindir.CS, "/mixmoment.sim", sep=""), what="character", nlines=1)
colnames(lambda.emerg.CS)    <- scan(paste(chaindir.CS, "/lambda.sim", sep=""), what="character", nlines=1)
colnames(logposter.emerg.CS) <- scan(paste(chaindir.CS, "/logposter.sim", sep=""), what="character", nlines=1)
colnames(beta.emerg.CS)      <- scan(paste(chaindir.CS, "/beta.sim", sep=""), what="character", nlines=1) 
colnames(mixmoment.b.emerg.CS) <- scan(paste(chaindir.CS, "/mixmoment_b.sim", sep=""), what="character", nlines=1)
colnames(lambda.b.emerg.CS)    <- scan(paste(chaindir.CS, "/lambda_b.sim", sep=""), what="character", nlines=1)
colnames(logposter.b.emerg.CS) <- scan(paste(chaindir.CS, "/logposter_b.sim", sep=""), what="character", nlines=1)

scale   <- sqrt(mixmoment.emerg.CS[, "D.1.1"])
scale.b <- sqrt(mixmoment.b.emerg.CS[, "D.1.1"])
intercept <- mixmoment.emerg.CS[,"Mean.1"] + mixmoment.b.emerg.CS[,"Mean.1"]
baseline.emerg.CS <- data.frame(Intcpt.16=intercept, Intcpt.26=intercept+beta.emerg.CS[, "Tooth26"], 
                                Intcpt.36=intercept+beta.emerg.CS[, "Tooth36"], Intcpt.46=intercept+beta.emerg.CS[, "Tooth46"], 
                                Scale=scale, Scale.b=scale.b)
beta.emerg.CS <- as.data.frame(beta.emerg.CS)
logposter.emerg.CS <- as.data.frame(logposter.emerg.CS)
logposter.b.emerg.CS <- as.data.frame(logposter.b.emerg.CS)
lambda.emerg.CS <- data.frame(lambda=lambda.emerg.CS[,1], lambda.b=lambda.b.emerg.CS[,1])
rm(list=c("intercept", "scale", "scale.b", "mixmoment.emerg.CS", "mixmoment.b.emerg.CS", "lambda.b.emerg.CS"))


###################################################
### chunk number 48: chains04
###################################################
str(baseline.emerg.CS)
str(beta.emerg.CS)
str(logposter.emerg.CS)
str(logposter.b.emerg.CS)
str(lambda.emerg.CS)


###################################################
### chunk number 49: chains05
###################################################
mixmoment.caries.CS <- matrix(scan(paste(chaindir.CS, "/mixmoment_2.sim", sep=""), skip=1, nlines=nlines), ncol=3, byrow=TRUE)
lambda.caries.CS    <- matrix(scan(paste(chaindir.CS, "/lambda_2.sim", sep=""), skip=1, nlines=nlines), ncol=1, byrow=TRUE)
logposter.caries.CS <- matrix(scan(paste(chaindir.CS, "/logposter_2.sim", sep=""), skip=1, nlines=nlines), ncol=3 , byrow=TRUE)
beta.caries.CS      <- matrix(scan(paste(chaindir.CS, "/beta_2.sim", sep=""), skip=1, nlines=nlines), ncol=nbeta.caries.CS, byrow=TRUE)
mixmoment.b.caries.CS <- matrix(scan(paste(chaindir.CS, "/mixmoment_b2.sim", sep=""), skip=1, nlines=nlines), ncol=3, byrow=TRUE)
lambda.b.caries.CS    <- matrix(scan(paste(chaindir.CS, "/lambda_b2.sim", sep=""), skip=1, nlines=nlines), ncol=1, byrow=TRUE)
logposter.b.caries.CS <- matrix(scan(paste(chaindir.CS, "/logposter_b2.sim", sep=""), skip=1, nlines=nlines), ncol=3 , byrow=TRUE)

colnames(mixmoment.caries.CS) <- scan(paste(chaindir.CS, "/mixmoment_2.sim", sep=""), what="character", nlines=1)
colnames(lambda.caries.CS)    <- scan(paste(chaindir.CS, "/lambda_2.sim", sep=""), what="character", nlines=1)
colnames(logposter.caries.CS) <- scan(paste(chaindir.CS, "/logposter_2.sim", sep=""), what="character", nlines=1)
colnames(beta.caries.CS)      <- scan(paste(chaindir.CS, "/beta_2.sim", sep=""), what="character", nlines=1) 
colnames(mixmoment.b.caries.CS) <- scan(paste(chaindir.CS, "/mixmoment_b2.sim", sep=""), what="character", nlines=1)
colnames(lambda.b.caries.CS)    <- scan(paste(chaindir.CS, "/lambda_b2.sim", sep=""), what="character", nlines=1)
colnames(logposter.b.caries.CS) <- scan(paste(chaindir.CS, "/logposter_b2.sim", sep=""), what="character", nlines=1)

scale   <- sqrt(mixmoment.caries.CS[, "D.1.1"])
scale.b <- sqrt(mixmoment.b.caries.CS[, "D.1.1"])
intercept <- mixmoment.caries.CS[,"Mean.1"] + mixmoment.b.caries.CS[,"Mean.1"]
baseline.caries.CS <- data.frame(Intcpt.16=intercept, Intcpt.26=intercept+beta.caries.CS[, "Tooth26"], 
                                 Intcpt.36=intercept+beta.caries.CS[, "Tooth36"], Intcpt.46=intercept+beta.caries.CS[, "Tooth46"], 
                                 Scale=scale, Scale.b=scale.b)
beta.caries.CS <- as.data.frame(beta.caries.CS)
logposter.caries.CS <- as.data.frame(logposter.caries.CS)
logposter.b.caries.CS <- as.data.frame(logposter.b.caries.CS)
lambda.caries.CS <- data.frame(lambda=lambda.caries.CS[,1], lambda.b=lambda.b.caries.CS[,1])
rm(list=c("intercept", "scale", "scale.b", "mixmoment.caries.CS", "mixmoment.b.caries.CS", "lambda.b.caries.CS"))


###################################################
### chunk number 50: chains06
###################################################
str(baseline.caries.CS)
str(beta.caries.CS)
str(logposter.caries.CS)
str(logposter.b.caries.CS)
str(lambda.caries.CS)


###################################################
### chunk number 51: summary01
###################################################
summ.baseline.emerg.CS <- give.summary(baseline.emerg.CS)
summ.baseline.caries.CS <- give.summary(baseline.caries.CS)
summ.beta.emerg.CS <- give.summary(beta.emerg.CS)
summ.beta.caries.CS <- give.summary(beta.caries.CS)
summ.lambda.emerg.CS <- give.summary(lambda.emerg.CS)
summ.lambda.caries.CS <- give.summary(lambda.caries.CS)


###################################################
### chunk number 52: summary02
###################################################
rsumm.baseline.emerg.CS <- round(summ.baseline.emerg.CS, dig=4)
rsumm.baseline.caries.CS <- round(summ.baseline.caries.CS, dig=4)
rsumm.beta.emerg.CS <- round(summ.beta.emerg.CS, dig=4)
rsumm.beta.caries.CS <- round(summ.beta.caries.CS, dig=4)
rsumm.lambda.emerg.CS <- round(summ.lambda.emerg.CS, dig=4)
rsumm.lambda.caries.CS <- round(summ.lambda.caries.CS, dig=4)


###################################################
### chunk number 53: summary03
###################################################
print(rsumm.beta.emerg.CS)
print(rsumm.baseline.emerg.CS)
print(rsumm.lambda.emerg.CS)


###################################################
### chunk number 54: summary04
###################################################
print(rsumm.beta.caries.CS)
print(rsumm.baseline.caries.CS)
print(rsumm.lambda.caries.CS)


###################################################
### chunk number 55: postDens01
###################################################
plfun <- "densplot2"
plname <- "dens"

pdfwidth <- 9
pdfheight <- 6
tt <- c(16, 26, 36, 46)


###################################################
### chunk number 56: postDens02
###################################################
pdf(paste(figuredir, plname, "_emergBeta.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.emerg.CS[, "Tooth26"]), main="Tooth26 (emerg)"))
  title(sub=paste("Tooth26 = ", rsumm.beta.emerg.CS[1, "Tooth26"], 
        " (", rsumm.beta.emerg.CS[3, "Tooth26"], ", ", rsumm.beta.emerg.CS[4, "Tooth26"], ")", sep=""))
  eval(call(plfun, mcmc(beta.emerg.CS[, "Tooth36"]), main="Tooth36 (emerg)"))
  title(sub=paste("Tooth36 = ", rsumm.beta.emerg.CS[1, "Tooth36"], 
        " (", rsumm.beta.emerg.CS[3, "Tooth36"], ", ", rsumm.beta.emerg.CS[4, "Tooth36"], ")", sep=""))
  eval(call(plfun, mcmc(beta.emerg.CS[, "Tooth46"]), main="Tooth46 (emerg)"))
  title(sub=paste("Tooth46 = ", rsumm.beta.emerg.CS[1, "Tooth46"], 
        " (", rsumm.beta.emerg.CS[3, "Tooth46"], ", ", rsumm.beta.emerg.CS[4, "Tooth46"], ")", sep=""))

  eval(call(plfun, mcmc(beta.emerg.CS[, "Girl"]), main="Girl (emerg)"))
  title(sub=paste("Girl = ", rsumm.beta.emerg.CS[1, "Girl"], 
        " (", rsumm.beta.emerg.CS[3, "Girl"], ", ", rsumm.beta.emerg.CS[4, "Girl"], ")", sep=""))
dev.off()


###################################################
### chunk number 57: postDens03
###################################################
pdf(paste(figuredir, plname, "_emergIntcpt.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(baseline.emerg.CS[, "Intcpt.16"]), main="Intcpt (emerg), tooth 16"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Intcpt.16"], 
        " (", rsumm.baseline.emerg.CS[3, "Intcpt.16"], ", ", rsumm.baseline.emerg.CS[4, "Intcpt.16"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.emerg.CS[, "Intcpt.26"]), main="Intcpt (emerg), tooth 26"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Intcpt.26"], 
        " (", rsumm.baseline.emerg.CS[3, "Intcpt.26"], ", ", rsumm.baseline.emerg.CS[4, "Intcpt.26"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.emerg.CS[, "Intcpt.36"]), main="Intcpt (emerg), tooth 36"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Intcpt.36"], 
        " (", rsumm.baseline.emerg.CS[3, "Intcpt.36"], ", ", rsumm.baseline.emerg.CS[4, "Intcpt.36"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.emerg.CS[, "Intcpt.46"]), main="Intcpt (emerg), tooth 46"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Intcpt.46"], 
        " (", rsumm.baseline.emerg.CS[3, "Intcpt.46"], ", ", rsumm.baseline.emerg.CS[4, "Intcpt.46"], ")", sep=""))
dev.off()


###################################################
### chunk number 58: postDens04
###################################################
pdf(paste(figuredir, plname, "_emergScaleLambda.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(baseline.emerg.CS[, "Scale"]), main="Scale (emerg)"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Scale"], 
        " (", rsumm.baseline.emerg.CS[3, "Scale"], ", ", rsumm.baseline.emerg.CS[4, "Scale"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.emerg.CS[, "Scale.b"]), main="Scale of b (emerg)"))
  title(sub=paste("Intercept = ", rsumm.baseline.emerg.CS[1, "Scale.b"], 
        " (", rsumm.baseline.emerg.CS[3, "Scale.b"], ", ", rsumm.baseline.emerg.CS[4, "Scale.b"], ")", sep=""))

  eval(call(plfun, mcmc(lambda.emerg.CS[, "lambda"]), main="Lambda (emerg)"))
  title(sub=paste("lambda = ", rsumm.lambda.emerg.CS[1, "lambda"],
        " (", rsumm.lambda.emerg.CS[3, "lambda"], ", ", rsumm.lambda.emerg.CS[4, "lambda"], ")", sep=""))

  eval(call(plfun, mcmc(lambda.emerg.CS[, "lambda.b"]), main="Lambda of b (emerg)"))
  title(sub=paste("lambda.b = ", rsumm.lambda.emerg.CS[1, "lambda.b"],
        " (", rsumm.lambda.emerg.CS[3, "lambda.b"], ", ", rsumm.lambda.emerg.CS[4, "lambda.b"], ")", sep="")) 
dev.off()


###################################################
### chunk number 59: postDens05
###################################################
pdf(paste(figuredir, plname, "_cariesBeta1.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.caries.CS[, "Tooth26"]), main="Tooth26 (caries)"))
  title(sub=paste("Tooth26 = ", rsumm.beta.caries.CS[1, "Tooth26"], 
        " (", rsumm.beta.caries.CS[3, "Tooth26"], ", ", rsumm.beta.caries.CS[4, "Tooth26"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.CS[, "Tooth36"]), main="Tooth36 (caries)"))
  title(sub=paste("Tooth36 = ", rsumm.beta.caries.CS[1, "Tooth36"], 
        " (", rsumm.beta.caries.CS[3, "Tooth36"], ", ", rsumm.beta.caries.CS[4, "Tooth36"], ")", sep=""))
  eval(call(plfun, mcmc(beta.caries.CS[, "Tooth46"]), main="Tooth46 (caries)"))
  title(sub=paste("Tooth46 = ", rsumm.beta.caries.CS[1, "Tooth46"], 
        " (", rsumm.beta.caries.CS[3, "Tooth46"], ", ", rsumm.beta.caries.CS[4, "Tooth46"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.CS[, "Girl"]), main="Girl (caries)"))
  title(sub=paste("Girl = ", rsumm.beta.caries.CS[1, "Girl"], 
        " (", rsumm.beta.caries.CS[3, "Girl"], ", ", rsumm.beta.caries.CS[4, "Girl"], ")", sep=""))
dev.off()


###################################################
### chunk number 60: postDens06
###################################################
pdf(paste(figuredir, plname, "_cariesBeta2.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(beta.caries.CS[, "Brush"]), main="Brush (caries)"))
  title(sub=paste("Brush = ", rsumm.beta.caries.CS[1, "Brush"], 
        " (", rsumm.beta.caries.CS[3, "Brush"], ", ", rsumm.beta.caries.CS[4, "Brush"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.CS[, "Seal"]), main="Seal (caries)"))
  title(sub=paste("Seal = ", rsumm.beta.caries.CS[1, "Seal"], 
        " (", rsumm.beta.caries.CS[3, "Seal"], ", ", rsumm.beta.caries.CS[4, "Seal"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.CS[, "Plaque"]), main="Plaque (caries)"))
  title(sub=paste("Plaque = ", rsumm.beta.caries.CS[1, "Plaque"], 
        " (", rsumm.beta.caries.CS[3, "Plaque"], ", ", rsumm.beta.caries.CS[4, "Plaque"], ")", sep=""))

  eval(call(plfun, mcmc(beta.caries.CS[, "Prim5"]), main="Prim5 (caries)"))
  title(sub=paste("Prim5 = ", rsumm.beta.caries.CS[1, "Prim5"], 
        " (", rsumm.beta.caries.CS[3, "Prim5"], ", ", rsumm.beta.caries.CS[4, "Prim5"], ")", sep=""))
dev.off()


###################################################
### chunk number 61: postDens07
###################################################
pdf(paste(figuredir, plname, "_cariesIntcpt.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(baseline.caries.CS[, "Intcpt.16"]), main="Intcpt (caries), tooth 16"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Intcpt.16"], 
        " (", rsumm.baseline.caries.CS[3, "Intcpt.16"], ", ", rsumm.baseline.caries.CS[4, "Intcpt.16"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.caries.CS[, "Intcpt.26"]), main="Intcpt (caries), tooth 26"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Intcpt.26"], 
        " (", rsumm.baseline.caries.CS[3, "Intcpt.26"], ", ", rsumm.baseline.caries.CS[4, "Intcpt.26"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.caries.CS[, "Intcpt.36"]), main="Intcpt (caries), tooth 36"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Intcpt.36"], 
        " (", rsumm.baseline.caries.CS[3, "Intcpt.36"], ", ", rsumm.baseline.caries.CS[4, "Intcpt.36"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.caries.CS[, "Intcpt.46"]), main="Intcpt (caries), tooth 46"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Intcpt.46"], 
        " (", rsumm.baseline.caries.CS[3, "Intcpt.46"], ", ", rsumm.baseline.caries.CS[4, "Intcpt.46"], ")", sep=""))
dev.off()


###################################################
### chunk number 62: postDens08
###################################################
pdf(paste(figuredir, plname, "_cariesScaleLambda.pdf", sep=""), width=pdfwidth, height=pdfheight)
par(mfrow=c(2,2), bty="n")
  eval(call(plfun, mcmc(baseline.caries.CS[, "Scale"]), main="Scale (caries)"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Scale"], 
        " (", rsumm.baseline.caries.CS[3, "Scale"], ", ", rsumm.baseline.caries.CS[4, "Scale"], ")", sep=""))

  eval(call(plfun, mcmc(baseline.caries.CS[, "Scale.b"]), main="Scale of b (caries)"))
  title(sub=paste("Intercept = ", rsumm.baseline.caries.CS[1, "Scale.b"], 
        " (", rsumm.baseline.caries.CS[3, "Scale.b"], ", ", rsumm.baseline.caries.CS[4, "Scale.b"], ")", sep=""))

  eval(call(plfun, mcmc(lambda.caries.CS[, "lambda"]), main="Lambda (caries)"))
  title(sub=paste("lambda = ", rsumm.lambda.caries.CS[1, "lambda"],
        " (", rsumm.lambda.caries.CS[3, "lambda"], ", ", rsumm.lambda.caries.CS[4, "lambda"], ")", sep=""))

  eval(call(plfun, mcmc(lambda.caries.CS[, "lambda.b"]), main="Lambda of b (caries)"))
  title(sub=paste("lambda.b = ", rsumm.lambda.caries.CS[1, "lambda.b"],
        " (", rsumm.lambda.caries.CS[3, "lambda.b"], ", ", rsumm.lambda.caries.CS[4, "lambda.b"], ")", sep="")) 
dev.off()


###################################################
### chunk number 63: tracePlot01
###################################################
plfun <- "traceplot2"
plname <- "trace"


