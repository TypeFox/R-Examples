## Simulating ED values
"simFct" <- function(noSim, edVal = c(10, 20, 50), type = c("non-parametric", "parametric"), 
response = c("bin", "con"), fct = LL.2(), coefVec, method = c("sp", "p", "np"), 
doseVec, nVec, pVec, rVec, resVar, pfct = fct, reference = NULL, span = NA, 
minmax = "response", lower = NULL, upper = NULL, seedVal = 200810201)
{
    method <- match.arg(method)
    response <- match.arg(response)
    type <- match.arg(type)

    set.seed(seedVal)

    lenData <- length(doseVec)  # replace lenpv?

    ## Parametric simulations
    ## Drawing random dose-response curves
    if (type == "parametric")
    {
        ## Model fit to simulate from
        if (response == "bin")
        {
            simMat <- rdrm(noSim, fct, coefVec, doseVec, yerror = "rbinom", ypar = nVec, onlyY = TRUE)
        } else {
            simMat <- rdrm(noSim, fct, coefVec, doseVec, ypar = c(0, sqrt(resVar)), onlyY = TRUE)
        }       
        ## drop = FALSE also in rdrm???
        print(simMat$y[1, ])        
    }

    ## Non-parametric simulations
    if (type == "non-parametric")
    {
        lenpv <- length(pVec)
    
        simMat <- matrix(NA, noSim, lenpv)
        if (response == "bin")
        {
            for (i in 1:noSim)
            {
                simMat[i, ] <- rbinom(lenpv, nVec, pVec)
            }
        } else {
            for (i in 1:noSim)
            {
                simMat[i, ] <- rnorm(lenpv, pVec, sqrt(resVar))
            }
        }
        simMat <- list(y = simMat)
        print(simMat$y[1, ])
    }
    

    lenev <- length(edVal)
    aicVec <- rep(NA, noSim)
    edMat <- array(NA, c(lenev, 3, noSim))
    mixVec <- rep(NA, noSim)
    spanVec <- rep(span, noSim)

    ## Generalized cross-validation criterion
    gcvFct <- function(doseVec, y)  # define function outside the i loop
    {
        gcvVec <- rep(NA, 20)
        for (j in 1:20)
        {
            tempLoess <- try(loess(y ~ doseVec, degree = 1, span = j/20), silent = TRUE)
            if (inherits(tempLoess, "try-error"))
            {
                gcvVec[j] <- NA
            } else {
                gcvVec[j] <- sum(residuals(tempLoess)^2) / (lenData - tempLoess$trace.hat)^2
            }
        }
        ((1:20)/20)[which.min(gcvVec)]
    }

    for (i in 1:noSim) 
    {
        ## Converting to proportions
        if (response == "bin")
        {
            y <- simMat$y[i, ] / nVec
        } else {
            y <- simMat$y[i, ]
        }
        
        ## Obtaining model-robust fit 
        if (method == "sp")
        {      
            parModel <- drm(y ~ doseVec, fct = pfct)
            
            if (is.na(span)) 
            {
                spanVec[i] <- gcvFct(doseVec, y)
            }
#            print(spanVec[i])
            loessModel <- loess(y ~ doseVec, degree = 1, span = spanVec[i])  # span = span

            tempModel <- mrdrm(parModel, loessModel)
            if (inherits(tempModel, "try-error"))
            {
                tempModel <- list(edMat = NA, mixing = NA, aic = NA)
            }  else {
                 aicVec[i] <- tempModel$gof[3]
                 mixVec[i] <- tempModel$lambda                 
                 edMat[, , i] <- ED(tempModel, edVal, interval = "approximate", minmax = minmax, 
                 lower = lower, upper = upper, display = FALSE)[, c(1:3)]
            }
        }
        if (method == "p")
        {
            if (response == "con")
            {
                tempModel <- try(drm(formula = y ~ doseVec, fct = pfct), silent = TRUE)
            } else {
                tempModel <- try(drm(formula = y ~ doseVec, weights = nVec, fct = pfct, type = "binomial"), 
                silent = TRUE)
            }
            if (inherits(tempModel, "try-error"))
            {
                edMat[, , i] <- NA 
                mixVec[i] <- NA
            } else {
                tempED <- ED(tempModel, edVal, display = FALSE, ci = "delta")[, c(1, 3, 4)]
                if (inherits(tempModel, "try-error"))
                {
                    edMat[, , i] <- NA
                    mixVec[i] <- NA
                } else {
                    edMat[, , i] <- tempED
                    mixVec[i] <- 0
                    aicVec[i] <- AIC(tempModel)
                }
            }
        }     
    }
    list(edArray = edMat, mixVec = mixVec, edVal = edVal, aicVec = aicVec, spanVec = spanVec)
}

## Calculating coverage percentage
coverFct <- function(mfit, simres, edVec = NULL)
{
    edVal <- simres$edVal
    if (is.null(edVec)) 
    {
        edVec <- ED(mfit, edVal, display = FALSE)[, 1]
    }

#    notNA <- sum(!is.na(simres$edArray[1, 2,]))
    lenem <- length(edVal)
    cpVec <- rep(NA, lenem)
    cplVec <- rep(NA, lenem)
    cpuVec <- rep(NA, lenem)        
    mvVec <- rep(NA, lenem)    
    mwVec <- rep(NA, lenem)
    notNA <- rep(NA, lenem)    
    names(cpVec) <- edVal
    for (i in 1:lenem)
    {
        notNA[i] <- sum((!is.na(simres$edArray[i, 2, ])) & (!is.na(simres$edArray[i, 3, ])))
        cplVec[i] <- sum(is.na(simres$edArray[i, 2, ]) & (simres$edArray[i, 3, ] > edVec[i]), na.rm = TRUE)
        cpuVec[i] <- sum(is.na(simres$edArray[i, 3, ]) & (simres$edArray[i, 2, ] < edVec[i]), na.rm = TRUE)        
        cpVec[i] <- sum((simres$edArray[i, 2, ] < edVec[i]) & (simres$edArray[i, 3, ] > edVec[i]), na.rm = TRUE) / notNA[i]
        mvVec[i] <- mean(simres$edArray[i, 1, ], na.rm = TRUE)
        mwVec[i] <- mean(simres$edArray[i, 3, ] - simres$edArray[i, 2, ], na.rm = TRUE)
    }
    list(coverage = cpVec, covLow = cplVec, covUp = cpuVec, true = edVec, mean = mvVec, width = mwVec, 
    notNAs = notNA, NAs = length(simres$edArray[1, 2,]) - notNA, mixingAverage = mean(simres$mixVec, na.rm = TRUE))
}


#misspec <- function(lambda, delta = 0.5)
#{
#    function(x)
#    {
#        Lfct <- function(x, mu, tau) {1/(1+exp(-((x-mu)/tau)))}
#        (1-lambda)*Lfct(x, 0.5, 0.1) + lambda*(delta*Lfct(x, 0.25, 0.05) + (1-delta)*Lfct(x, 0.75, 0.05))
#    }
#}
#
#msFct <- misspec(0.5)

if (FALSE)
{

## Simulations based on the design and probabilities in the dataset 'deguelin' 

deguelin.m1 <- drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")

## Semi-parametric
simres<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = deguelin$dose, nVec = deguelin$n)  # , rVec = deguelin$r)
coverFct(deguelin.m1, simres)

true.sr.sp1b<-simFct2(100, edVal = c(50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = deguelin$dose, nVec = deguelin$n)  # , rVec = deguelin$r)
coverFct(deguelin.m1, true.sr.sp1)

simres4<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = c(49, deguelin$n))  # , rVec = c(1, deguelin$r))
coverFct(deguelin.m1, simres4)

simres7<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7)) 
coverFct(deguelin.m1, simres7)

simres7b<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), seedVal=200802211) 
coverFct(deguelin.m1, simres7b)

simres7c<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), seedVal=200804011) 
coverFct(deguelin.m1, simres7c)

simres7d<-simFct2(10, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("sp"), 
doseVec = c(2.5, deguelin$dose), nVec = rep(20, 7), seedVal=200804012) 
coverFct(deguelin.m1, simres7d)


## Parametric
simres2<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("p"), 
doseVec = deguelin$dose, nVec = deguelin$n)  # , rVec = deguelin$r)

simres3<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = c(49, deguelin$n))  # , rVec = c(1, deguelin$r))
coverFct(deguelin.m1, simres3)

simres5<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7)) 
coverFct(deguelin.m1, simres5)

simres6<-simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7)) 
coverFct(deguelin.m1, simres6)

## Non-parametric
np.simres1 <- simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7), seedVal = 200802191) 
coverFct(deguelin.m1, np.simres1)

np.simres2 <- simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), seedVal = 200802192) 
coverFct(deguelin.m1, np.simres2)

np.simres3 <- simFct2(1000, edVal = c(10,20,50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(deguelin.m1), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(50, 7), seedVal = 200802193) 
coverFct(deguelin.m1, np.simres3)



## Under misspecification

msFct <- misspec(0.5)

evFct <- function(edVal, maxx = 1)
{
    lenev <- length(edVal)
    edVec <- rep(NA, lenev)
    for (i in 1:lenev)
    {
        edVec[i] <- uniroot(function(x){msFct(x/maxx)-edVal[i]/100}, c(1, 99))$root 
    }
    edVec
}
edVec <- evFct(c(10,20,50), maxx = 52)


## Parametric
mis.sr.p1 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802194)
coverFct(deguelin.m1, mis.sr.p1, edVec)

mis.sr.p2 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802195)
coverFct(deguelin.m1, mis.sr.p2, edVec)

mis.sr.p3 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), 
doseVec = c(1, deguelin$dose), nVec = rep(50, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802196)
coverFct(deguelin.m1, mis.sr.p3, edVec)


## Semi-parametric
mis.sr.sp1 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802197)
coverFct(deguelin.m1, mis.sr.sp1, edVec)

mis.sr.sp1b <- simFct2(100, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802197)
coverFct(deguelin.m1, mis.sr.sp1b, edVec)

mis.sr.sp2 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802198)
coverFct(deguelin.m1, mis.sr.sp2, edVec)

mis.sr.sp3 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(50, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802199)
coverFct(deguelin.m1, mis.sr.sp3, edVec)

mis.sr.sp3b <- simFct2(100, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = c(1, deguelin$dose), nVec = rep(50, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802199)
coverFct(deguelin.m1, mis.sr.sp3b, edVec)


## Non-parametric
mis.sr.np1 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(10, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802201)
coverFct(deguelin.m1, mis.sr.np1, edVec)

mis.sr.np2 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(20, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802202)
coverFct(deguelin.m1, mis.sr.np2, edVec)

mis.sr.np3 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), 
doseVec = c(1, deguelin$dose), nVec = rep(50, 7), pVec =  msFct(c(1, deguelin$dose)/52), seedVal = 200802203)
coverFct(deguelin.m1, mis.sr.np3, edVec)



## Using 14 dose levels

dose14 <- seq(1, 50, length.out = 14)

## Semi-parametric
mis.14.sp1 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = dose14, nVec = rep(10, 14), pVec =  msFct(dose14/52), seedVal = 200802212)
coverFct(deguelin.m1, mis.14.sp1, edVec)

mis.14.sp2 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = dose14, nVec = rep(20, 14), pVec =  msFct(dose14/52), seedVal = 200802213)
coverFct(deguelin.m1, mis.14.sp2, edVec)

mis.14.sp3 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = dose14, nVec = rep(50, 14), pVec =  msFct(dose14/52), seedVal = 200802214)
coverFct(deguelin.m1, mis.14.sp3, edVec)

## Using 56 dose levels

dose56 <- seq(1, 50, length.out = 56)

## Semi-parametric
mis.56.sp1 <- simFct2(1000, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = dose56, nVec = rep(10, 56), pVec =  msFct(dose56/52), seedVal = 200802281)
coverFct(deguelin.m1, mis.56.sp1, evFct(c(10,20,50), maxx = 52))

mis.56.sp2 <- simFct2(100, edVal = c(10,20,50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = dose56, nVec = rep(20, 56), pVec =  msFct(dose56/52), seedVal = 200802281)
coverFct(deguelin.m1, mis.56.sp2, evFct(c(10,20,50), maxx = 52))




## bin.mat
bin.mat.m1 <- drm(matured/total~conc, weights=total, data = bin.mat[c(1,4,7,10,13),], fct=LL.2())

## Parametric
bin.mat.true.p1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("p"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, seedVal=200802271)
coverFct(bin.mat.m1, bin.mat.true.p1)

bin.mat.true.p2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("p"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), seedVal=200802272)
coverFct(bin.mat.m1, bin.mat.true.p2)

bin.mat.true.p3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("p"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), seedVal=200802273)
coverFct(bin.mat.m1, bin.mat.true.p3)


## Semi-parametric
true.sp1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("sp"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, seedVal=200802221)
coverFct(bin.mat.m1, true.sp1)

#true.aic.sp1<-true.sp1<-simFct2(10, edVal = c(10, 20, 50), type = c("parametric"), 
#response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("sp"), 
#doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, seedVal=200802221,aic=TRUE)
#true.aic.sp1$aic
#
#true.aic.p1<-true.sp1<-simFct2(10, edVal = c(10, 20, 50), type = c("parametric"), 
#response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("p"), 
#doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, seedVal=200802221,aic=TRUE)
#true.aic.p1$aic
## the AIC values are not compatible between the parametric and the semi-parametric models for binomial data!


true.sp2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("sp"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), seedVal=200802222)
coverFct(bin.mat.m1, true.sp2)

true.sp3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), fct = LL.2(), coefVec = coef(bin.mat.m1), method = c("sp"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), seedVal=200802223)
coverFct(bin.mat.m1, true.sp3)

## Non-parametric
bin.mat.true.np1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), method = c("np"), doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, 
fct = LL.2(), coefVec = coef(bin.mat.m1), 
seedVal=200802261)
coverFct(bin.mat.m1, bin.mat.true.np1)

bin.mat.true.np2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), method = c("np"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), 
fct = LL.2(), coefVec = coef(bin.mat.m1),
seedVal=200802262)
coverFct(bin.mat.m1, bin.mat.true.np2)

bin.mat.true.np3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("bin"), method = c("np"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), 
fct = LL.2(), coefVec = coef(bin.mat.m1),
seedVal=200802263)
coverFct(bin.mat.m1, bin.mat.true.np3)


## Under misspecification

evFct <- function(edVal, maxx = 1)
{
    lenev <- length(edVal)
    edVec <- rep(NA, lenev)
    for (i in 1:lenev)
    {
        edVec[i] <- uniroot(function(x){msFct(x/maxx)-edVal[i]/100}, c(1e-6, maxx - 1e-6))$root 
    }
    edVec
}
edVec <- evFct(c(10,20,50))


## Parametric
bin.mat.mis.p1<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, 
pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802251)
coverFct(bin.mat.m1, bin.mat.mis.p1, edVec)

bin.mat.mis.p2<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802252)
coverFct(bin.mat.m1, bin.mat.mis.p2, edVec)

bin.mat.mis.p3<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("p"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802253)
coverFct(bin.mat.m1, bin.mat.mis.p3, edVec)

## Semi-parametric
bin.mat.mis.sp1<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, 
pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802254)
coverFct(bin.mat.m1, bin.mat.mis.sp1, edVec)

bin.mat.mis.sp2<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802255)
coverFct(bin.mat.m1, bin.mat.mis.sp2, edVec)

bin.mat.mis.sp3<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("sp"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802256)
coverFct(bin.mat.m1, bin.mat.mis.sp3, edVec)

## Non-parametric
bin.mat.mis.np1<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = bin.mat[c(1,4,7,10,13),]$total, 
pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802257)
coverFct(bin.mat.m1, bin.mat.mis.np1, edVec)

bin.mat.mis.np2<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(20, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802258)
coverFct(bin.mat.m1, bin.mat.mis.np2, edVec)

bin.mat.mis.np3<-simFct2(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("bin"), method = c("np"), 
doseVec = bin.mat[c(1,4,7,10,13),]$conc, nVec = rep(50, 5), pVec = 1 - msFct(bin.mat[c(1,4,7,10,13),]$conc), 
seedVal=200802259)
coverFct(bin.mat.m1, bin.mat.mis.np3, edVec)




## ryegrass
ryegrass.m1 <- drm(rootl~conc, data = ryegrass, fct=LL.4())

## Parametric
ryegrass.true.p1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("p"), 
doseVec = unique(ryegrass$conc), resVar = summary(ryegrass.m1)$resVar, seedVal=200803031)
coverFct(ryegrass.m1, ryegrass.true.p1)

ryegrass.true.p2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("p"), 
doseVec = rep(unique(ryegrass$conc), rep(2, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal=200803032)
coverFct(ryegrass.m1, ryegrass.true.p2)

ryegrass.true.p3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("p"), 
doseVec = rep(unique(ryegrass$conc), rep(3, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal=200803033)
coverFct(ryegrass.m1, ryegrass.true.p3)

## Semi-parametric

## Tester
ryegrass.tester <- simFct(50, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("sp"), 
doseVec = unique(ryegrass$conc), resVar = summary(ryegrass.m1)$resVar, seedVal=200810202)

ryegrass.true.sp1<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("sp"), 
doseVec = unique(ryegrass$conc), resVar = summary(ryegrass.m1)$resVar, seedVal = 200810211)
coverFct(ryegrass.m1, ryegrass.true.sp1)



dVec <- c(0, 0.12, 0.235, 0.47, unique(ryegrass$conc)[-1], 60, 120, 240, 480)
## No. replicates: 1
ryegrass.true.sp1new<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  dVec, resVar = 0.25, seedVal = 200810213)
coverFct(ryegrass.m1, ryegrass.true.sp1new, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp.1<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec = dVec, resVar = 0.25, seedVal = 200810213, span = 0.75)
coverFct(NULL, true.sp.1, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp.1a<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec = dVec, resVar = 0.0025, seedVal = 200810213, span = 0.75)
coverFct(NULL, true.sp.1a, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))


## More help by fixing lower and upper limits
true.sp.1b <- simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec = dVec, resVar = 0.0025, seedVal = 200810213, lower = 0.5, upper = 8, span = 0.75)
coverFct(NULL, true.sp.1ab, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp.1ab <- simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec = dVec, resVar = 0.0025, seedVal = 200810213, lower = 0.05, upper = 0.8, span = 0.75)
coverFct(NULL, true.sp.1ab, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.p.1<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("p"), 
doseVec = dVec, resVar = 0.0025, seedVal = 200810191)
coverFct(NULL, true.p.1, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## No. replicates: 3
ryegrass.true.sp2new<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(3, 14)), resVar = 0.25, seedVal = 200810213)
coverFct(ryegrass.m1, ryegrass.true.sp2new, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp.3<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(3, 14)), resVar = 0.25, seedVal = 200811081, span = 0.75)
coverFct(NULL, true.sp.3, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## More help by fixing lower and upper limits
true.sp.3b <- simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec = rep(dVec, rep(3, 14)), resVar = 0.25, seedVal = 200811081, lower = 0.5, upper = 8, span = 0.75)
coverFct(NULL, true.sp.3b, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## Response between 0 and 1
true.sp.3a<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(3, 14)), resVar = 0.0025, seedVal = 200811081, span = 0.75)
coverFct(NULL, true.sp.3a, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp.3ab <- simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec = rep(dVec, rep(3, 14)), resVar = 0.0025, seedVal = 200811081, lower = 0.05, upper = 0.8, span = 0.75)
coverFct(NULL, true.sp.3ab, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.p.3<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("p"), 
doseVec =  rep(dVec, rep(3, 14)), resVar = 0.0025, seedVal = 200810192)
coverFct(ryegrass.m1, true.p.3, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## No. replicates: 5
ryegrass.true.sp3new<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.25, seedVal = 200810214)
coverFct(ryegrass.m1, ryegrass.true.sp3new, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp5<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.25, seedVal = 200811081, span = 0.75)
coverFct(NULL, true.sp5, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## More help by fixing lower and upper limits
true.sp5b<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.5, 8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.25, seedVal = 200811081, lower = 0.5, upper = 8, span = 0.75)
coverFct(NULL, true.sp5b, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

## Response between 0 and 1
true.sp5a<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.0025, seedVal = 200811081, span = 0.75)
coverFct(NULL, true.sp5a, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.sp5ab<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("sp"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.0025, seedVal = 200811081, lower = 0.05, upper = 0.8, span = 0.75)
coverFct(NULL, true.sp5ab, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))

true.p.5<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = c(3, 0.05, 0.8, 60), method = c("p"), 
doseVec =  rep(dVec, rep(5, 14)), resVar = 0.0025, seedVal = 200810193)
coverFct(ryegrass.m1, true.p.5, c(60*(10/90)^(1/3), 60*(20/80)^(1/3), 60))



#misFct <- function(a, b, c, d, x0, x)
#{
#    indX <- (x < x0)
#    x1 <- x[indX]
#    x2 <- x[!indX]
#    
#    c(a + b * x1, c + d * x2)   
#}
# misFct(8, -0.02, 4, -0.025, 30, c(0:13)*10)

#misFct2 <- function(a, b, c, d, e, x0, x)
#{
#    indX1 <- (x < x0)
#    indX2 <- (x >= x0) & (x < x0 + 0.025)
#    indX3 <- (x >= x0 + 0.025)        
#    x1 <- x[indX1]
#    x2 <- x[indX2]
#    x3 <- x[indX3]
#    
#    c(a + 0 * x1, b + c * x2, d + e * x3)   
#}

## Misspecified model
#dVec2 <- c(0.06, 0.12, 0.235, 0.47, unique(ryegrass$conc)[-1], 60, 120, 240, 480)
#fVec <- misFct2(8, 4.5+3.5/0.025*0.1, -3.5/0.025, 1+3.5/499.9*500, -3.5/499.9, 0.075, dVec2)
##coverFct(NULL, mis.sp.1, c(0.080715, 0.0864285, 71.5))
#fVec <- misFct2(8, 5+3/0.025*0.1, -3/0.025, 1+4/499.9*500, -4/499.9, 0.075, dVec2)
#
### No. replicates:1
#mis.sp.1 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
#response = c("con"), method = c("sp"),  doseVec = rep(dVec2, rep(1, 14)), 
#pVec = rep(fVec, rep(1, 14)),  seedVal = 200810222, pfct = LL.4(), resVar = 0.5)
#
#coverFct(NULL, mis.sp.1, c(0.0816666, 0.0883333, 125))
#
#
#mis.p.1 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
#response = c("con"), method = c("p"),  doseVec = rep(dVec2, rep(1, 14)), 
#pVec = rep(fVec, rep(1, 14)),  seedVal = 200810222, pfct = LL.4(), resVar = 0.5)
#
#coverFct(NULL, mis.p.1, c(0.0816666, 0.0883333, 125))


dVec3<-seq(0, 500, length = 14)
horFct<-function(x){(8+0.1*x)/(1+(x/125)^3.5)}
fVec <- horFct(dVec3)
dVec4<-seq(0, 500, length = 50)
fVec2 <- horFct(dVec4)
dVec5<-seq(0, 500, length = 100)
fVec3 <- horFct(dVec5)


## No. replicates:1
mis.sp.1ab <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(1, 14)), 
pVec = rep(fVec/10, rep(1, 14)),  seedVal = 200710221, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose")
coverFct(NULL, mis.sp.1ab, c(158.96, 169.39, 211.35))

## More help by fixing lower and upper limits
mis.sp.1aab <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(1, 14)), 
pVec = rep(fVec/10, rep(1, 14)),  seedVal = 200710221, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose",
lower = 0, upper = 0.8)
coverFct(NULL, mis.sp.1aab, c(158.96, 169.39, 211.35))

mis.p.1 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("p"),  doseVec = rep(dVec3, rep(1, 14)), 
pVec = rep(fVec, rep(1, 14)),  seedVal = 200710222, pfct = LL.4(), resVar = 1)
coverFct(NULL, mis.p.1, c(158.96, 169.39, 211.35))

## No. replicates:1 - 100 doses
mis.sp.1.100a <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec5, rep(1, 100)), 
pVec = rep(fVec3/10, rep(1, 100)), seedVal = 200811032, pfct = LL.4(), resVar = 0.01, span = 0.2, minmax = "dose", 
lower = 0, upper = 0.8)
coverFct(NULL, mis.sp.1.100a, c(158.96, 169.39, 211.35))

mis.p.1 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("p"),  doseVec = rep(dVec3, rep(1, 14)), 
pVec = rep(fVec, rep(1, 14)),  seedVal = 200710222, pfct = LL.4(), resVar = 1)
coverFct(NULL, mis.p.1, c(158.96, 169.39, 211.35))


## No. replicates:3
mis.sp.3a <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(3, 14)), 
pVec = rep(fVec/10, rep(3, 14)), seedVal = 200810216, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose")
coverFct(NULL, mis.sp.3a, c(158.96, 169.39, 211.35))

## More help by fixing lower and upper limits
mis.sp.3ab <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(3, 14)), 
pVec = rep(fVec/10, rep(3, 14)),  seedVal = 200810216, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose",
lower = 0, upper = 0.8)
coverFct(NULL, mis.sp.3ab, c(158.96, 169.39, 211.35))

mis.p.3 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("p"),  doseVec = rep(dVec3, rep(3, 14)), 
pVec = rep(fVec/10, rep(3, 14)),  seedVal = 200810227, pfct = LL.4(), resVar = 0.01)
coverFct(NULL, mis.p.3, c(158.96, 169.39, 211.35))

mis.sp.3.100 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec5, rep(3, 100)), 
pVec = rep(fVec3, rep(3, 100)), seedVal = 200811053, pfct = LL.4(), resVar = 1, span = 0.2, minmax = "dose", 
lower = 0, upper = 8)
coverFct(NULL, mis.sp.3.100, c(158.96, 169.39, 211.35))



## No. replicates:5
mis.sp.5a <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(5, 14)), 
pVec = rep(fVec/10, rep(5, 14)), seedVal = 200810228, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose")
coverFct(NULL, mis.sp.5, c(158.96, 169.39, 211.35))

## More help by fixing lower and upper limits
mis.sp.5ab <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("sp"),  doseVec = rep(dVec3, rep(5, 14)), 
pVec = rep(fVec/10, rep(5, 14)),  seedVal = 200810228, pfct = LL.4(), resVar = 0.01, span = 0.35, minmax = "dose",
lower = 0, upper = 0.8)
coverFct(NULL, mis.sp.5b, c(158.96, 169.39, 211.35))

mis.p.5 <- simFct(1000, edVal = c(10, 20, 50), type = c("non-parametric"), 
response = c("con"), method = c("p"),  doseVec = rep(dVec3, rep(5, 14)), 
pVec = rep(fVec/10, rep(5, 14)),  seedVal = 200810229, pfct = LL.4(), resVar = 0.01)
coverFct(NULL, mis.p.5, c(158.96, 169.39, 211.35))


## Figure 2
horFct2<-function(x){horFct(x)/10}
ll4Fct<-function(x){(0.5+(8-0.5)/(1+(x/60)^3))/10}
curve(ll4Fct, xlim=c(0, 500), ylim=c(0, 1.4), xlab="Dose", ylab="Response", lwd=2)
curve(horFct2, xlim=c(0, 500), add=TRUE, lty=2, lwd=2)






ryegrass.true.sp2<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("sp"), 
doseVec = rep(unique(ryegrass$conc), rep(2, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal = 200810212)
coverFct(ryegrass.m1, ryegrass.true.sp2)

ryegrass.true.sp3<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("sp"), 
doseVec = rep(unique(ryegrass$conc), rep(3, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal = 200810206)
coverFct(ryegrass.m1, ryegrass.true.sp3)

ryegrass.true.sp5<-simFct(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("sp"), 
doseVec = rep(unique(ryegrass$conc), rep(5, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal = 200810207)
coverFct(ryegrass.m1, ryegrass.true.sp5)


## Non-parametric
ryegrass.true.np1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("np"), 
doseVec = unique(ryegrass$conc), resVar = summary(ryegrass.m1)$resVar, seedVal=200803037)
coverFct(ryegrass.m1, ryegrass.true.np1)

ryegrass.true.np2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("np"), 
doseVec = rep(unique(ryegrass$conc), rep(2, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal=200803038)
coverFct(ryegrass.m1, ryegrass.true.np2)

ryegrass.true.np3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(ryegrass.m1), method = c("np"), 
doseVec = rep(unique(ryegrass$conc), rep(3, 7)), resVar = summary(ryegrass.m1)$resVar, seedVal=200803039)
coverFct(ryegrass.m1, ryegrass.true.np3)

## Under misspecification
lettuce.m1 <- drm(weight ~ conc, data = lettuce, fct = LL.3())
lettuce.m2 <- drm(weight ~ conc, data = lettuce, fct = BC.4())

## Parametric
lettuce.mis.p1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("p"), 
doseVec = unique(lettuce$conc), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803041)
coverFct(lettuce.m2, lettuce.mis.p1)

lettuce.mis.p2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("p"), 
doseVec = rep(unique(lettuce$conc), rep(2, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803042)
coverFct(lettuce.m2, lettuce.mis.p2)

lettuce.mis.p3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("p"), 
doseVec = rep(unique(lettuce$conc), rep(3, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803043)
coverFct(lettuce.m2, lettuce.mis.p3)

## Semi-parametric
lettuce.mis.sp1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("sp"), 
doseVec = unique(lettuce$conc), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), reference = 0, seedVal=200803044)
coverFct(lettuce.m2, lettuce.mis.sp1)

lettuce.mis.sp2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("sp"), 
doseVec = rep(unique(lettuce$conc), rep(2, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), 
reference = 0, seedVal=200803045)
coverFct(lettuce.m2, lettuce.mis.sp2)

lettuce.mis.sp3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("sp"), 
doseVec = rep(unique(lettuce$conc), rep(3, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), 
reference = 0, seedVal=200803046)
coverFct(lettuce.m2, lettuce.mis.sp3)

## Non-parametric
lettuce.mis.np1<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("np"), 
doseVec = unique(lettuce$conc), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803047)
coverFct(lettuce.m2, lettuce.mis.np1)

lettuce.mis.np2<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("np"), 
doseVec = rep(unique(lettuce$conc), rep(2, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803048)
coverFct(lettuce.m2, lettuce.mis.np2)

lettuce.mis.np3<-simFct2(1000, edVal = c(10, 20, 50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("np"), 
doseVec = rep(unique(lettuce$conc), rep(3, 7)), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), seedVal=200803049)
coverFct(lettuce.m2, lettuce.mis.np3)

## Comparison of parametric and semi-parametric models

## Under the true model
exp.a.m1<-drm(y~x, data=exp.a, fct=LL.4())

true.aic.sp1<-true.sp1<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(exp.a.m1), method = c("sp"), 
doseVec = unique(exp.a$x), seedVal=200805052, resVar = summary(exp.a.m1)$resVar, aic = TRUE)
true.aic.sp1$aic

true.aic.p1<-true.p1<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("con"), fct = LL.4(), coefVec = coef(exp.a.m1), method = c("p"), 
doseVec = unique(exp.a$x), seedVal=200805052, resVar = summary(exp.a.m1)$resVar, aic = TRUE)
true.aic.p1$aic

## Under misspecification
mis.aic.sp1<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("sp"), 
doseVec = seq(min(uniVec<-unique(lettuce$conc)), max(uniVec), length.out=14), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), reference = 0, 
seedVal=200805053, aic = TRUE)
mis.aic.sp1$aic

mis.aic.p1<-simFct2(1000, edVal = c(50), type = c("parametric"), 
response = c("con"), fct = BC.4(), coefVec = coef(lettuce.m2), method = c("p"), 
doseVec = seq(min(uniVec<-unique(lettuce$conc)), max(uniVec), length.out=14), resVar = summary(lettuce.m2)$resVar, pfct = LL.3(), reference = 0, 
seedVal=200805053, aic = TRUE)
mis.aic.p1$aic
plot(mis.aic.p1$aic-(mis.aic.sp1$aic+2))


## Testing area
mr.test <- mrdrm(bin.mat[c(1,4,7,10,13),]$conc, c(12, 5, 4, 2, 0)/bin.mat[c(1,4,7,10,13),]$total, 
bin.mat[c(1,4,7,10,13),]$total, 
type = "bin", fct = LL.2(), respLev = c(10, 20, 50), reference = NULL, level = 0.95, robust = FALSE, 
mixVec = seq(0, 1, by = 0.05), logex = TRUE, bwLower = 0)

mr.test1 <- mrdrm(lettuce$conc, lettuce$weight, 
type = "con", fct = LL.3(), respLev = c(10, 20, 50), reference = 0, level = 0.95, robust = FALSE, 
mixVec = seq(0, 1, by = 0.05), bwLower = 0, logex = TRUE)
mr.test1$edMat
mr.test1$mixing
plotmr(mr.test1)


mr.test2 <- mrdrm(lettuce$conc, lettuce$weight, 
type = "con", fct = LL.3(), respLev = c(10, 20, 50), reference = 0, level = 0.95, robust = TRUE, 
mixVec = seq(0, 1, by = 0.05), bwLower = 0, logex = TRUE)
mr.test2$edMat
mr.test2$mixing
plotmr(mr.test2)

mr.test3 <- mrdrm(lettuce$conc, lettuce$weight, 
type = "con", fct = BC.4(), respLev = c(10, 20, 50), reference = 0, level = 0.95, robust = FALSE, 
mixVec = seq(0, 1, by = 0.05), bwLower = 0, logex = TRUE)
mr.test3$edMat
mr.test3$mixing
plotmr(mr.test3)



simres<-simFct1(1000, 50)
edMat <- ED(deguelin.m1, c(10, 20, 50), reference="control")[, 1]

## Cannot be estimated
#sum(simres$edArray[2,,1] > edMat[1])
#sum(simres$edArray[3,,1] < edMat[1])
#sum(simres$edArray[2,,2] > edMat[2])
#sum(simres$edArray[3,,2] < edMat[2])

sum(simres$edArray[2,,1] > edMat[3])
sum(simres$edArray[3,,1] < edMat[3])

mean(simres$edArray[1,,1])
mean(apply(simres$edArray[2:3,,1], 2, diff))

## Simulations based on the design and probabilities modified from the dataset 'deguelin' 
simres2<-simFct1(1000, design="moddeguelin", type="parametric", seedVal=20080103)

sum(simres2$edArray[2,,3] > 9.8228)
sum(simres2$edArray[3,,3] < 9.8228)

sum(simres2$edArray[2,,2] > 4.7186)
sum(simres2$edArray[3,,2] < 4.7186)

sum(simres2$edArray[2,,1] > 3.0729)
sum(simres2$edArray[3,,1] < 3.0729)

mean(simres2$edArray[1,,1])
mean(simres2$edArray[1,,2])
mean(simres2$edArray[1,,3])

mean(apply(simres2$edArray[2:3,,3], 2, diff))
mean(apply(simres2$edArray[2:3,,2], 2, diff))
mean(apply(simres2$edArray[2:3,,1], 2, diff))


        tempModel1 <- SP.mrr(
        formula = deguelin$r/deguelin$n ~ log(deguelin$dose),
        cases=deguelin$n,
        NP.control=NP.control.lr.wls(
        response.type="bin",
        optim="bandwidth.grid",
        weight.function="gaus",
        bandwidth.type="fixed.width",
        bandwidth=seq(0.3,1.0,by=0.05)
        ),
        P.control=P.control.glm.binomial.ml(
        link="logit"
        ),
        SP.control=SP.control.mrr(
        response.type="bin",
        optim="mixing.grid",
        mixing=seq(0,1,by=0.05)
        ))

        tempModel2 <- SP.mrr(
        formula = c(1, deguelin$r)/c(49, deguelin$n) ~ log(c(1, deguelin$dose)),
        cases=c(49, deguelin$n),
        NP.control=NP.control.lr.wls(
        response.type="bin",
        optim="bandwidth.grid",
        weight.function="gaus",
        bandwidth.type="fixed.width",
        bandwidth=seq(0.3,1.0,by=0.05)
        ),
        P.control=P.control.glm.binomial.ml(
        link="logit"
        ),
        SP.control=SP.control.mrr(
        response.type="bin",
        optim="mixing.grid",
        mixing=seq(0,1,by=0.05)
        ))
        
        m1<-drm(c(1, deguelin$r)/c(49, deguelin$n) ~ c(1, deguelin$dose), 
        weights = c(49, deguelin$n), fct=LL.2(), type="binomial")
        
        plot(c(1,deguelin$dose), predict(m1)[,1], type="l")
        lines(c(1,deguelin$dose), predict(m1)[,1]+predict(m1)[,2]*1.96, lty=2)
        lines(c(1,deguelin$dose), predict(m1)[,1]-predict(m1)[,2]*1.96, lty=2)
        abline(h=0.5, lty=3)
        
        m2<-glm(c(1, deguelin$r)/c(48, deguelin$n) ~ log(c(1, deguelin$dose)), 
        weights = c(48, deguelin$n), family=binomial)        
}
