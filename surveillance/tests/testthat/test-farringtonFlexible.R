data("salmonella.agona")
# sts object
lala <- paste(salmonella.agona$start[1],salmonella.agona$start[2],"1",sep=" ")
firstMonday <- as.POSIXlt(lala, format = "%Y %W %u")
salm.ts <- salmonella.agona$observed
dates <- as.Date(firstMonday) + 7 * 0:(length(salm.ts) - 1)
start=c(salmonella.agona$start[1],salmonella.agona$start[2])
salm <- new("sts",epoch = as.numeric(dates), start = start, freq = 52,
observed = salm.ts, epochAsDate = TRUE)

################################################################################
context("weights function")
################################################################################
test_that("gamma = 1 if everything below the threshold",{
  s <- rep(0,10)
  weightsThreshold <- 0
  weights <- algo.farrington.assign.weights(s,weightsThreshold)
  expect_equal(weights,rep(1,10))   
})


test_that(" A case that was checked by hand",{
  s <- rep(2,10)
  s[1:5] <- 0
  weightsThreshold <- 0
  weights <- algo.farrington.assign.weights(s,weightsThreshold)
  expect_equal(weights[1:5],rep(1.6,5))
  expect_equal(weights[6:10],rep(0.4,5))
})
################################################################################
# END OF WEIGHTS FUNCTION TESTS
################################################################################

################################################################################
context("residuals function")
################################################################################

test_that(" residuals should be zero",{
  x <- rpois(10,1)
  y <- exp(x)
  model <- glm(y~x,family = quasipoisson(link="log"))
  phi <- max(summary(model)$dispersion,1)
  s <- anscombe.residuals(model,phi)
  expect_equal(as.numeric(s),rep(0,10))
})

test_that(" residuals should not be zero",{
  x <- rpois(1000,1)
  y <- exp(x)+runif(1)
  model <- glm(y~x,family = quasipoisson(link="log"))
  phi <- max(summary(model)$dispersion,1)
  s <- anscombe.residuals(model,phi)
  expect_true(mean(s)>0)
})
################################################################################
# END OF RESIDUALS FUNCTION TESTS
################################################################################

################################################################################
context("formula function")
################################################################################
test_that("We get the right formula",{
  expect_equal(formulaGLM(populationOffset=FALSE,timeBool=TRUE,factorsBool=FALSE),"response ~ 1+wtime")
  expect_equal(formulaGLM(populationOffset=FALSE,timeBool=FALSE,factorsBool=FALSE),"response ~ 1")
  expect_equal(formulaGLM(populationOffset=TRUE,timeBool=TRUE,factorsBool=FALSE),"response ~ 1+wtime+offset(log(population))")
  expect_equal(formulaGLM(populationOffset=TRUE,timeBool=TRUE,factorsBool=TRUE),"response ~ 1+wtime+offset(log(population))+seasgroups")

})
################################################################################
# END OF FORMULA FUNCTION TESTS
################################################################################

################################################################################
context("reference time points function")
################################################################################
test_that("We get the expected timepoints with weekly data",{
  # Case with weekly data with dates
  dayToConsider <- as.Date("2013-06-06")
  b <- 3
  freq <- 52
  epochAsDate <- TRUE
  epochStr <- "week"
  lala <- algo.farrington.referencetimepoints(dayToConsider,b=b,freq=freq,epochAsDate,epochStr)
  # Do we get the same day as dayToConsider?
  expect_equal(as.numeric(format(lala, "%w")),rep(4,4))
  # Actually for this example I know the dates one should get
expect_equal(sort(lala),sort(c(as.Date("2010-06-03"),as.Date("2013-06-06"),as.Date("2012-06-07"),as.Date("2011-06-09"))))
})

test_that("We get the expected timepoints with monthly data",{
  dayToConsider <- 48
  b <- 3
  freq <- 12
  epochAsDate <- FALSE
  epochStr <- "month"
  lala <- algo.farrington.referencetimepoints(dayToConsider,b=b,freq=freq,epochAsDate,epochStr)
  expect_equal(lala,c(48,36,24,12))
})
test_that("one gets a warning if too many years back",{
  dayToConsider <- 48
  b <- 3
  freq <- 12
  epochAsDate <- FALSE
  epochStr <- "month"
  expect_that(algo.farrington.referencetimepoints(dayToConsider,b=8,freq=freq,epochAsDate,epochStr), gives_warning("Some reference"))
  
  # apply code
   control1 <-  list(range=250,noPeriods=10,populationOffset=FALSE,
                     fitFun="algo.farrington.fitGLM.flexible",
                     b=10,w=3,weightsThreshold=2.58,
                     pastWeeksNotIncluded=26,
                     pThresholdTrend=1,trend=TRUE,
                     thresholdMethod="muan",alpha=0.05,glmWarnings=FALSE)
  expect_error(farringtonFlexible(salm,control=control1),"Some reference")
})
################################################################################
# END OF REFERENCE TIME POINTS FUNCTION TESTS
################################################################################

################################################################################
context("fit glm function")
################################################################################

# Case with convergence
control<-  list(range=250,noPeriods=10,populationOffset=TRUE,
                fitFun="algo.farrington.fitGLM.flexible",
                b=40,w=3,weightsThreshold=2.58,
                pastWeeksNotIncluded=26,
                pThresholdTrend=1,trend=TRUE,
                thresholdMethod="muan",alpha=0.05,glmWarnings=FALSE)
 response=salm@observed[1:120]     
dataGLM <- data.frame(response=response,wtime=1:120,
	                  population=runif(120)*100,
                      seasgroups=as.factor(rep(1:12,10)))
     
arguments <- list(dataGLM=dataGLM,
                   timeTrend=TRUE,
                   populationOffset=TRUE,
                   factorsBool=TRUE,reweight=TRUE,
                   weightsThreshold=0.5,glmWarnings=control$glmWarnings,
				   control=control)
model <- do.call(algo.farrington.fitGLM.flexible, args=arguments)

test_that("The fit glm function gives the right class of output?",{
  expect_equal(class(model),c("glm","lm"))
})

test_that("The fit glm function gives as many coefficients as expected",{
  expect_equal(dim(summary(model)$coefficients)[1],length(levels(dataGLM$seasgroups))-1+1+1)
})

test_that("wtime, response, phi and weights were added to the model",{
  expect_true(is.null(model$phi)==FALSE)
  expect_true(is.null(model$wtime)==FALSE)
  expect_true(is.null(model$response)==FALSE)
  expect_true(is.null(model$population)==FALSE)
  expect_true(is.null(model$weights)==FALSE)
})

test_that("reweighting was done",{
  expect_true(sum(model$weights!=rep(1,length(model$weights)))==length(model$weights))
})

test_that("there are no weights if very high threshold",{
  arguments$reweight <- TRUE
  arguments$weightsThreshold <- 100000
  model <- do.call(algo.farrington.fitGLM.flexible, args=arguments)
  expect_true(sum(model$weights==rep(1,length(model$weights)))==length(model$weights))
})

test_that("there is not a too small overdispersion",{
  expect_true(model$phi>=1)
})

################################################################################
# END OF FIT GLM FUNCTION TESTS
################################################################################

################################################################################
context("block function")
################################################################################

referenceTimePoints <- c(as.Date("2010-06-03"),as.Date("2013-06-06"),as.Date("2012-06-07"),as.Date("2011-06-09"))
firstDay <- as.Date("1990-06-07")
vectorOfDates <- dates <- as.Date(firstDay) + 7 * 0:1300
freq <- 52
dayToConsider <- as.Date("2013-06-06")
b <- 3
w <- 3
epochAsDate <- TRUE

# p=1
p <- 1
lala <- blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,
epochAsDate)
test_that("the reference window has the right length",{
  expect_equal(length(vectorOfDates[is.na(lala)==FALSE&lala==p]),w+1+b*(2*w+1))
  
  # p>1
  p <- 8
  lala <- blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,
  epochAsDate)
  # reference windows
  expect_equal(length(vectorOfDates[is.na(lala)==FALSE&lala==p]),w+1+b*(2*w+1))
})

lili <- as.factor(lala[is.na(lala)==FALSE])

test_that("there are as many levels as expected",{ 
  expect_equal(length(levels(lili)),p)

})
p <- 8
lala <- blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,
               epochAsDate)
lili <- as.factor(lala[is.na(lala)==FALSE])
lolo <- lili[lili!=p]
test_that("periods of roughly the same length each year",{ 
    expect_equal(as.numeric(abs(diff(table(lolo))[1:(p-2)])<=b),rep(1,(p-2)))
})
################################################################################
# END OF BLOCKS FUNCTION TESTS
################################################################################

################################################################################
context("Farrington threshold function")
################################################################################

predFit <- 5
predSeFit <- 0.2
wtime <- 380
skewness.transform <- "2/88"
alpha <- 0.05
y <- 8
method <- "delta"
phi <- 1

test_that("the function recognizes wrong exponents",{
  expect_that(algo.farrington.threshold.farrington(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method),throws_error("proper exponent"))
})

test_that("some results we know are found",{
  skewness.transform <- "none"
  lala <- algo.farrington.threshold.farrington(predFit,predSeFit,phi,
                                                 skewness.transform,
												  alpha,y,method)												  
  # Should always be ok
  lala <- as.numeric(lala)
  expect_true(lala[3]<=1&lala[1]>=0)	
  expect_true(lala[2]>lala[1])
  expect_true(lala[1]>=0)			
  
  # Here we know the results								  
  expect_equal(abs(as.numeric(lala)-c(1.3073128, 8.6926872, 0.0907246, 0.8124165))<rep(1e-7,4),rep(TRUE,4))
  
  skewness.transform <- "1/2"
  lala <- algo.farrington.threshold.farrington(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method)												  
  											  
  expect_equal(abs(as.numeric(lala)-c( 1.9891097, 9.3744842, 0.0000000, 0.6857951))<rep(1e-7,4),rep(TRUE,4))
  
  skewness.transform <- "2/3"
  lala <- algo.farrington.threshold.farrington(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method)												  
  											  
  expect_equal(abs(as.numeric(lala)-c( 1.808448e+00,  9.115482e+00, 1.596176e-112,  7.289546e-01))<rep(1e-6,4),rep(TRUE,4))

})
################################################################################
# END OF THRESHOLD FUNCTION FARRINGTON TESTS
################################################################################


################################################################################
context("Noufaily threshold function")
################################################################################
											  
predFit <- log(5)
predSeFit <- log(2)
wtime <- 380
skewness.transform <- "none"
alpha <- 0.05
y <- 11
phi <- 1.5
method <- "muan"
lala <- algo.farrington.threshold.noufaily(predFit,predSeFit,phi,
                                                 skewness.transform,
												  alpha,y,method)
test_that("some results we know are found",{
  # Should always be ok
  lala <- as.numeric(lala)
  expect_true(lala[3]<=1&lala[1]>=0)	
  expect_true(lala[2]>lala[1])
  expect_true(lala[1]>=0)	
  
  # Here we calculated some examples
  expect_equal(abs(as.numeric(lala)-c(7.0000000, 26.0000000,  0.8597797,  0.3850080))<rep(1e-6,4),rep(TRUE,4))
  phi <- 1.0
  method <- "muan"
  lala <- algo.farrington.threshold.noufaily(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method)
  expect_equal(abs(as.numeric(lala)-c(8.0000000, 24.0000000,  0.9093099 , 0.4193982))<rep(1e-6,4),rep(TRUE,4))
  
  phi <- 1.5
  method <- "nbPlugin"
  lala <- algo.farrington.threshold.noufaily(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method)
  expect_equal(abs(as.numeric(lala)-c(1.00000000, 11.00000000,  0.03763657,  1.00000000))<rep(1e-6,4),rep(TRUE,4))
  												  
  phi <- 1.0
  method <- "nbPlugin"
  lala <- algo.farrington.threshold.noufaily(predFit,predSeFit,phi,
                                                   skewness.transform,
  												  alpha,y,method)										
  expect_equal(abs(as.numeric(lala)-c( 1.00000000, 10.00000000,  0.01369527,  1.11918153))<rep(1e-6,4),rep(TRUE,4))
})
################################################################################
# END OF THRESHOLD FUNCTION NOUFAILY TESTS
################################################################################

################################################################################
context("data GLM function")
################################################################################
b <- 3
freq <- 52
dayToConsider <- as.Date("2013-05-30")
epochAsDate <- TRUE
epochStr <- "week"
firstDay <- as.Date("1990-06-07")
vectorOfDates <- dates <- as.Date(firstDay) + 7 * 0:1300
w <- 3
noPeriods <- 10
observed <- rnorm(1301)+runif(1301)+30
population <- rnorm(1301)+10
verbose <- FALSE
pastWeeksNotIncluded <- w
k <- 1200

lala <- algo.farrington.data.glm(dayToConsider, b, freq, 
                                     epochAsDate,epochStr,
									 vectorOfDates,w,noPeriods,
									 observed,population,
									 verbose,pastWeeksNotIncluded,k)#
test_that("the output is a data.frame",{
  expect_true(class(lala)=="data.frame")
})

test_that("the data frame contains all variables",{
  expect_equal(names(lala)==c( "response", "wtime","population","seasgroups","vectorOfDates"),rep(TRUE,5))
})

test_that("the time variable is ok with diff 1",{
  expect_equal(diff(lala$wtime)==rep(1,length(lala$wtime)-1),rep(TRUE,length(lala$wtime)-1))
})

test_that("the factor variable has the right number of levels",{
  expect_true(length(levels(lala$seasgroups))==noPeriods)
})

observed[1150] <- NA
lala <- algo.farrington.data.glm(dayToConsider, b, freq, 
                                     epochAsDate,epochStr,
									 vectorOfDates,w,noPeriods,
									 observed,population,
									 verbose,pastWeeksNotIncluded,k)

test_that("the data frame has the right dimensions",{
  expect_equal(dim(lala),c(156,5))
})
################################################################################
# END OF DATA GLM FUNCTION TESTS
################################################################################

################################################################################
context("GLM function")
################################################################################

dataGLM <- lala
timeTrend <- TRUE
populationOffset <- TRUE
factorsBool <- TRUE
reweight <- TRUE
weightsThreshold <- 1
pThresholdTrend <- 1
b <- 3
noPeriods <- 10
typePred <- "link"
fitFun <- "algo.farrington.fitGLM.flexible"
glmWarnings <- FALSE
epochAsDate <- TRUE
dayToConsider <- as.Date("2013-05-30")
diffDates <- 7
populationNow <- 10

test_that("the output has the needed variables",{
  finalModel <- algo.farrington.glm(dataGLM,timeTrend,populationOffset,factorsBool,
                                  reweight,weightsThreshold,pThresholdTrend,b,
  								noPeriods,typePred,fitFun,glmWarnings,epochAsDate,
  								dayToConsider,diffDates,populationNow,verbose=FALSE)
  expect_equal(names(finalModel)==c("pred","doTrend","coeffTime","phi"),rep(TRUE,4))
})

test_that("no time trend in no time trend",{
  pThresholdTrend <- 1
  b <- 2
  finalModel <- algo.farrington.glm(dataGLM,timeTrend,populationOffset,factorsBool,
                                  reweight,weightsThreshold,pThresholdTrend,b,
  								noPeriods,typePred,fitFun,glmWarnings,epochAsDate,
  								dayToConsider,diffDates,populationNow,verbose=FALSE)
  expect_equal(finalModel$doTrend,FALSE)
})
################################################################################
# END OF GLM FUNCTION TESTS
################################################################################

################################################################################
context("Alarms")
################################################################################
data("salmonella.agona")
# sts object
lala <- paste(salmonella.agona$start[1],salmonella.agona$start[2],"1",sep=" ")
firstMonday <- as.POSIXlt(lala, format = "%Y %W %u")
salm.ts <- salmonella.agona$observed
dates <- as.Date(firstMonday) + 7 * 0:(length(salm.ts) - 1)
start=c(salmonella.agona$start[1],salmonella.agona$start[2])
salm <- new("sts",epoch = as.numeric(dates), start = start, freq = 52,
            observed = salm.ts, epochAsDate = TRUE)
test <- farringtonFlexible(salm,control=list(thresholdMethod="nbPlugin",alpha=0.1))
test_that("there are only alarms when expected",{
  # No alarm when observed is 0
  expect_true(sum(test@alarm[test@observed==0])==0)
  # No alarm when the observed counts are UNDER the threshold
  expect_true(sum(observed(test)>upperbound(test),na.rm=TRUE)==sum(test@alarm==TRUE))
})
################################################################################
# RESIDUALS FUNCTION
################################################################################

################################################################################
context("no convergence")
################################################################################
timeSeries <- rep(0,698)
timeSeries[696] <- 1

algoControl <- list(noPeriods=10,alpha = 0.01,verbose = F,
                    b=5,w=4,weightsThreshold=2.58,pastWeeksNotIncluded=26,
                    pThresholdTrend=1,thresholdMethod='nbPlugin',limit54 = c(4,5),
                    range = (length(timeSeries) - 1):length(timeSeries), glmWarnings = FALSE)
seriesSTSObject <- new('sts', observed = timeSeries,
                       epoch = as.numeric(seq(as.Date('2001-01-01'),length.out=length(timeSeries), by='1 week')),
                       epochAsDate = TRUE)
test_that("The code does not produce any error",{
# It is ok if the code does not produce any error
 expect_that(farringtonFlexible(seriesSTSObject, control = algoControl), gives_warning())
})

################################################################################
context("NA")
################################################################################
timeSeries <- observed <- rnorm(698)*10+runif(698)*100+30

algoControl <- list(noPeriods=10,alpha = 0.01,verbose = F,
                    b=5,w=4,weightsThreshold=2.58,pastWeeksNotIncluded=w,
                    pThresholdTrend=1,thresholdMethod='nbPlugin',limit54 = c(4,5),
                    range = (length(timeSeries) - 1):length(timeSeries), glmWarnings = FALSE)
seriesSTSObject <- new('sts', observed = timeSeries,
                       epoch = as.numeric(seq(as.Date('2001-01-01'),length.out=length(timeSeries), by='1 week')),
                       epochAsDate = TRUE)
test_that("The code does not produce any error",{
  farringtonFlexible(seriesSTSObject, control = algoControl)

  results1 <- farringtonFlexible(seriesSTSObject, control = algoControl)
  expect_that(results1, is_a("sts"))
  seriesSTSObject@observed[680:690] <- NA
  results2 <- farringtonFlexible(seriesSTSObject, control = algoControl)
  expect_that(results2, is_a("sts"))
})

