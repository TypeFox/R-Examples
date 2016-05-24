require("DoseFinding")

## ## commented out for time reasons

## resp <- c(1.23, 1.31, 1.32, 1.36, 1.38)
## dose <- c(0, 1.25, 2.5, 5, 10)
## sdev <- c(0.015, 0.014, 0.015, 0.016, 0.015)
## V <- diag(sdev^2)
## mods <- Mods(emax=c(2.65, 12.5), linear=NULL, linInt = c(1, 1, 1, 1),
##              logistic=c(29, 9.55), quadratic = -0.0075,
##              doses=dose)
## mmfit <- MCPMod(dose, resp, S=V, type="general", models=mods, Delta=0.12)
## efit <- mmfit$mods$emax
## plot(efit, plotData = "meansCI", CI=TRUE)
## plot(efit, plotData = "meansCI", CI=FALSE)
## ## plot(efit, plotData = "raw") # should throw an error
## plot(efit, plotData = "means", CI = TRUE)
## plot(efit, plotData = "means", CI = FALSE)
## plot(efit, plotData = "none", CI =TRUE)
## plot(efit, plotData = "none", CI =FALSE)

## plot(mmfit, plotData = "meansCI", CI=TRUE)
## plot(mmfit, plotData = "meansCI", CI=FALSE)
## ## plot(mmfit, plotData = "raw") # should throw an error
## plot(mmfit, plotData = "means", CI = TRUE)
## plot(mmfit, plotData = "means", CI = FALSE)
## plot(mmfit, plotData = "none", CI =TRUE)
## plot(mmfit, plotData = "none", CI =FALSE)

## data(IBScovars)
## models <- Mods(emax = c(0.5, 1), betaMod=c(1,1), linear = NULL, doses=c(0,4))
## mmfit <- MCPMod(dose, resp, data=IBScovars, models=models, Delta=0.12)
## efit <- mmfit$mods$emax
## plot(efit, plotData = "meansCI", CI=TRUE)
## plot(efit, plotData = "meansCI", CI=FALSE)
## plot(efit, plotData = "raw", CI=FALSE)
## plot(efit, plotData = "raw", CI=TRUE)
## plot(efit, plotData = "means", CI = TRUE)
## plot(efit, plotData = "means", CI = FALSE)
## plot(efit, plotData = "none", CI =TRUE)
## plot(efit, plotData = "none", CI =FALSE)

## plot(mmfit, plotData = "meansCI", CI=TRUE)
## plot(mmfit, plotData = "meansCI", CI=FALSE)
## plot(mmfit, plotData = "raw", CI=TRUE)
## plot(mmfit, plotData = "raw", CI=FALSE)
## plot(mmfit, plotData = "means", CI = TRUE)
## plot(mmfit, plotData = "means", CI = FALSE)
## plot(mmfit, plotData = "none", CI =TRUE)
## plot(mmfit, plotData = "none", CI =FALSE)

## data(IBScovars)
## models <- Mods(emax = c(0.5, 1), betaMod=c(1,1), linInt = c(1, 1, 1, 1),
##                linear = NULL, doses=0:4)
## anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
## drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
## vCov <- vcov(anovaMod)[2:5,2:5]
## dose <- sort(unique(IBScovars$dose))[-1]
## mmfit <- MCPMod(dose, drFit, S=vCov, type = "general", models=models, Delta=0.12, placAdj=TRUE)
## efit <- mmfit$mods$emax
## plot(efit, plotData = "meansCI", CI=TRUE)
## plot(efit, plotData = "meansCI", CI=FALSE)
## ## plot(efit, plotData = "raw", CI=FALSE) # should throw an error
## ## plot(efit, plotData = "raw", CI=TRUE) # should throw an error
## plot(efit, plotData = "means", CI = TRUE)
## plot(efit, plotData = "means", CI = FALSE)
## plot(efit, plotData = "none", CI =TRUE)
## plot(efit, plotData = "none", CI =FALSE)

## plot(mmfit, plotData = "meansCI", CI=TRUE)
## plot(mmfit, plotData = "meansCI", CI=FALSE)
## ## plot(mmfit, plotData = "raw", CI=TRUE) # should throw an error
## ## plot(mmfit, plotData = "raw", CI=FALSE) # should throw an error
## plot(mmfit, plotData = "means", CI = TRUE)
## plot(mmfit, plotData = "means", CI = FALSE)
## plot(mmfit, plotData = "none", CI =TRUE)
## plot(mmfit, plotData = "none", CI =FALSE)

## ## neurodeg example (in 0.9-6 not all means were visible)
## doses <- c(0,1,3,10,30)
## muH  <- c(-5.099, -4.581, -3.22, -2.879, -3.52)     # estimated slope
## covH <- structure(c(0.149, 0.009, 0.009, 0.009, 0.009, 0.009, 0.149, 
##                     0.009, 0.009, 0.009, 0.009, 0.009, 0.149, 0.009,
##                     0.009, 0.009, 0.009, 0.009, 0.149, 0.009, 0.009,
##                     0.009, 0.009, 0.009, 0.149), .Dim = c(5L, 5L))
## fit <- fitMod(doses, muH, S=covH, model="emax", type = "general")
## plot(fit, plotData="meansCI", CI=TRUE)

