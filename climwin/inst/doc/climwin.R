## ------------------------------------------------------------------------
library(climwin)

## ---- eval = FALSE-------------------------------------------------------
#  
#  MassWin <- climatewin(xvar = list(Temp = MassClimate$Temp),
#                        cdate = MassClimate$Date,
#                        bdate = Mass$Date,
#                        baseline = lm(Mass ~ 1, data = Mass),
#                        cinterval = "day",
#                        furthest = 100, closest = 0,
#                        type = "fixed", cutoff.day = 20, cutoff.month = 05,
#                        stat = "mean",
#                        func = "lin")
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  head(MassWin[[1]]$Dataset)
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  MassWin[[1]]$BestModel
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  Call:
#  lm(formula = Yvar ~ climate, data = modeldat)
#  
#  Coefficients:
#  (Intercept)      climate
#      163.544       -4.481
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  head(MassWin[[1]]$BestModelData)
#  

## ---- fig.width = 4, fig.height = 4--------------------------------------

plotdelta(dataset = MassOutput)


## ---- fig.width = 4, fig.height = 4--------------------------------------

plotweights(dataset = MassOutput)


## ---- fig.width = 4, fig.height = 4--------------------------------------

plotbetas(dataset = MassOutput)


## ---- eval = FALSE-------------------------------------------------------
#  
#  MassRand <- randwin(repeats = 5,
#                      xvar = list(Temp = MassClimate$Temp),
#                      cdate = MassClimate$Date,
#                      bdate = Mass$Date,
#                      baseline = lm(Mass ~ 1, data = Mass),
#                      cinterval = "day",
#                      furthest = 100, closest = 0,
#                      type = "fixed", cutoff.day = 20, cutoff.month = 05,
#                      stat = "mean",
#                      func = "lin")

## ---- fig.width = 4, fig.height = 4--------------------------------------

plothist(dataset = MassOutput, datasetrand = MassRand)


## ---- fig.width = 4, fig.height = 4--------------------------------------

plotwin(dataset = MassOutput)


## ------------------------------------------------------------------------

MassSingle <- singlewin(xvar = list(Temp = MassClimate$Temp),
                        cdate = MassClimate$Date,
                        bdate = Mass$Date,
                        baseline = lm(Mass ~ 1, data = Mass),
                        cinterval = "day",
                        furthest = 72, closest = 15,
                        type = "fixed", cutoff.day = 20, cutoff.month = 05,
                        stat = "mean",
                        func = "lin")

## ---- fig.width = 4, fig.height = 4--------------------------------------

plotbest(dataset = MassOutput,
         bestmodel = MassSingle$BestModel, 
         bestmodeldata = MassSingle$BestModelData)


## ---- fig.width = 10, fig.height = 7.5-----------------------------------

plotall(dataset = MassOutput,
        datasetrand = MassRand,
        bestmodel = MassSingle$BestModel, 
        bestmodeldata = MassSingle$BestModelData)


## ---- eval = FALSE-------------------------------------------------------
#  
#  MassWin2 <- climatewin(xvar = list(Temp = MassClimate$Temp),
#                         cdate = MassClimate$Date,
#                         bdate = Mass$Date,
#                         baseline = lm(Mass ~ 1, data = Mass),
#                         cinterval = "day",
#                         furthest = 100, closest = 0,
#                         type = "fixed", cutoff.day = 20, cutoff.month = 05,
#                         stat = c("max", "mean"),
#                         func = c("lin", "quad"))
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  MassWin2$combos
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  MassWin2[[3]]$BestModel
#  

## ---- eval = FALSE-------------------------------------------------------
#  Call:
#  lm(formula = Yvar ~ climate + I(climate^2), data = modeldat)
#  
#  Coefficients:
#   (Intercept)       climate  I(climate^2)
#     139.39170      -1.33767       0.03332

