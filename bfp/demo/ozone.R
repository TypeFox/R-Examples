## load packages and the ozone data set
library(Hmisc)
library(bfp)
data(ozone)

## look at the ozone data
describe (ozone)
names (ozone)

## strong correlations between two variables:
symnum (cor (subset (na.omit (ozone), select = - weekday)))

## so we drop tempElMonte
ozoneFull <- na.omit (subset (ozone, select = - tempElMonte))

## add a workday variable
workday <- ozoneFull$weekday
levels (workday)
levels (workday) <- c (rep ("workday", 5), rep ("weekend", 2))
ozoneFull$workday <- workday

## and a day of the year variable
library (doBy)
ozoneFull <- orderBy (~ month + day, ozoneFull)
ozoneFull$dayOfYear <- seq_len (nf <- nrow (ozoneFull))

## determine test and training data set
set.seed (321)                          
whichTest <- sample (nf, 30, replace = FALSE)
ozoneTest <- ozoneFull [whichTest, ]
ozoneTraining <- ozoneFull[- whichTest, ]

ozoneTest <- orderBy (~ month + day, ozoneTest)
testDates <- with (ozoneTest, paste (month, day, sep = "/"))
testDates <- matrix (testDates, nrow = 3, ncol = 10, byrow = TRUE)

## run sampler
set.seed (432)
system.time(ozoneModels <- BayesMfp (hourAverageMax ~
                                     bfp (dayOfYear) +
                                     bfp (pressure500Height) +
                                     bfp (windSpeed) +
                                     bfp (humidity) +
                                     bfp (tempSandburg) +
                                     bfp (inversionBaseHeight) +
                                     bfp (pressureGradientDaggett) +
                                     bfp (inversionBaseTemp) +
                                     bfp (visibility),
                                     data = ozoneTraining,
                                     nModels = 3000,
                                     method = "sampling",
                                     chainlength = 2e+5))

## get a summary of the results
sumOzoneModels <- summary(ozoneModels)
str(sumOzoneModels)

## the top 50 models:
sumOzoneModels$dataframe[1:50,]

## now we add the workday variable
set.seed (567)
system.time(ozoneModelsWork <- BayesMfp(hourAverageMax ~
                                        bfp (dayOfYear) +
                                        bfp (pressure500Height) +
                                        bfp (windSpeed) +
                                        bfp (humidity) +
                                        bfp (tempSandburg) +
                                        bfp (inversionBaseHeight) +
                                        bfp (pressureGradientDaggett) +
                                        bfp (inversionBaseTemp) +
                                        bfp (visibility) +
                                        uc (workday),
                                        data = ozoneTraining,
                                        nModels = 3000,
                                        method = "sampling",
                                        chainlength = 2e+5))

## what is the MAP (maximum a posteriori) model?
mapModel <- ozoneModels[1]
summary (mapModel)

## get nice plots of the curve estimates conditional on the MAP model and
## the posterior expected g
myLty <- c (1, 2, 2, 3, 3)
myLwd <- 2

post <- getPosteriorParms (mapModel)
plotCurveEstimate (mapModel, "dayOfYear", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[0]^2~paste("(", x[0], ",",~alpha, ",", ~p[0],")")),
                   xlab = expression (z[0]),
                   post = post,
                   legendPos = "topright")
plotCurveEstimate (mapModel, "humidity",  lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[6]^1~paste("(", x[6], ",",~alpha, ",", ~p[6],")")),
                   xlab = expression (z[6]),
                   post = post)
plotCurveEstimate (mapModel, "tempSandburg", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[7]^1~paste("(", x[7], ",",~alpha, ",", ~p[7],")")),
                   xlab = expression (z[7]),
                   post = post)
plotCurveEstimate (mapModel, "inversionBaseHeight", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[8]^1~paste("(", x[8], ",",~alpha, ",", ~p[8],")")),
                   xlab = expression (z[8]),
                   post = post)
plotCurveEstimate (mapModel, "pressureGradientDaggett", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[9]^1~paste("(", x[9], ",",~alpha, ",", ~p[9],")")),
                   xlab = expression (z[9]),
                   post = post)
plotCurveEstimate (mapModel, "visibility", lty = myLty, col = 1, lwd = myLwd,
                   ylab = expression (f[10]^2~paste("(", x[10], ",",~alpha, ",", ~p[10],")")),
                   xlab = expression (z[10]),
                   post = post,
                   legendPos = "topright")

## get point predictions for the test data set
mapPredictions <- predict (mapModel, newdata = ozoneTest)

## the corresponding root mean squared error is:
mapRMSE <- sqrt (mean ((mapPredictions - ozoneTest$hourAverageMax)^2))



## now do Bayesian model averaging over all saved models
set.seed (456)
system.time(BmaModelSamples <- BmaSamples (ozoneModels, gridSize = 0))
sumBma <- summary (BmaModelSamples)
print (sumBma, table = FALSE)

## histograms of intercept, regression variance and shrinkage parameter samples:
hist (BmaModelSamples$fixed, nclass = 100)
hist (BmaModelSamples$sigma2, nclass = 100)
hist(BmaModelSamples$shrinkage, nclass=100)

## again nice plots of curve estimates
plotCurveEstimate (BmaModelSamples, "dayOfYear")
plotCurveEstimate (BmaModelSamples, "humidity")
plotCurveEstimate (BmaModelSamples, "tempSandburg")
plotCurveEstimate (BmaModelSamples, "inversionBaseHeight")
plotCurveEstimate (BmaModelSamples, "pressureGradientDaggett")
plotCurveEstimate (BmaModelSamples, "visibility")

## get corresponding point predictions and RMSE
bmaPredictions <- bmaPredict (ozoneModels, newdata = ozoneTest)
bmaRMSE <- sqrt (mean ((bmaPredictions - ozoneTest$hourAverageMax)^2))
bmaRMSE



## try a corresponding analysis with the mfp package:
library (mfp)
ozoneMfp <- mfp (hourAverageMax ~
                 fp (dayOfYear) +
                 fp (pressure500Height) +
                 fp (windSpeed) +
                 fp (humidity) +
                 fp (tempSandburg) +
                 fp (inversionBaseHeight) +
                 fp (pressureGradientDaggett) +
                 fp (inversionBaseTemp) +
                 fp (visibility),
                 data = ozoneTraining)
ozoneMfp
str (ozoneMfp)

sumMfp <- summary (ozoneMfp)
str (sumMfp)

## get predictions from the found model:
ozoneMfpLm <- lm(ozoneMfp$formula, data = ozoneTraining)
mfpPredictions <- predict (ozoneMfpLm, newdata = ozoneTest)
mfpRMSE <- sqrt (mean ((mfpPredictions - ozoneTest$hourAverageMax)^2))


## and finally another model proposed by Casella & Moreno (for the full dataset admittedly)

casella <- lm (hourAverageMax ~ humidity + I (windSpeed^2) + I (tempSandburg^2) + I (pressureGradientDaggett^2) +
               factor (month):visibility + pressure500Height:tempSandburg + windSpeed:visibility +
               humidity:inversionBaseHeight,
               data = ozoneTraining)
summary (casella)

casellaPredictions <- predict (casella, newdata = ozoneTest)
casellaRMSE <- sqrt (mean ((casellaPredictions - ozoneTest$hourAverageMax)^2))


## make a graph comparing the prediction performance of the different approaches
predMat <- cbind(casellaPredictions, mfpPredictions, bmaPredictions, mapPredictions)

myPch <- c (3, 4, 19, 19)
myCol <- c (rep("black", 2), "grey", "black")

table (ozoneTest$hourAverageMax)
x <- jitter (ozoneTest$hourAverageMax)

matplot (x,
         predMat,
         pch = myPch,
         col = myCol,
         type = "p",
         xlim = (lims <- c (0, max (ozoneTest$hourAverageMax)) * 1.1),
         ylim = lims,
         xlab = "actual 1-hour-average ozone level [ppm]",
         ylab = "predicted value")
rug (x)
abline (0, 1)
legend ("bottomright", pch = myPch, col = myCol, legend = c ("Casella", "mfp", "BMA", "MAP"), bty = "n")




