library("testthat")
library("surveillance")
##################################################################
context("Checking the provided reporting triangle")
# Control slot for the proposed algorithm with D=10 correction
rangeTest <- 410:412
alpha <- 0.05

controlDelay <-  list(range = rangeTest, b = 4, w = 3,
                      pastAberrations = TRUE, mc.munu=10, mc.y=10,
                      verbose = FALSE,populationOffset=FALSE,
                      alpha = alpha, trend = TRUE,
                      limit54=c(0,50),
                      noPeriods = 10, pastWeeksNotIncluded = 26,
                      delay=TRUE)
test_that("The absence of reporting triangle throws an error",{
  data("salmNewport")
  expect_error(bodaDelay(salmNewport, controlDelay),"You have to")
})
test_that("The function spots uncorrect reporting triangles",{
  data('salmAllOnset')
  stsFake <- salmAllOnset
  stsFake@control$reportingTriangle$n <- head(stsFake@control$reportingTriangle$n,n=10)
  expect_error(bodaDelay(stsFake, controlDelay),"The reporting triangle number")
  stsFake <- salmAllOnset
  stsFake@control$reportingTriangle$n[1,] <- stsFake@control$reportingTriangle$n[1,]/2
  expect_error(bodaDelay(stsFake, controlDelay),"The reporting triangle is wrong")
})
##################################################################
context("Data glm function")

# Parameters
epochAsDate <- TRUE
epochStr <- "week"
freq <- 52
b <- controlDelay$b
w <- controlDelay$w
populationOffset <- controlDelay$populationOffset
noPeriods <- controlDelay$noPeriods
verbose <- controlDelay$verbose
reportingTriangle <- salmAllOnset@control$reportingTriangle
timeTrend <- controlDelay$trend
alpha <- controlDelay$alpha
populationOffset <- controlDelay$populationOffset
factorsBool <- controlDelay$factorsBool
pastAberrations <- controlDelay$pastAberrations
glmWarnings <- controlDelay$glmWarnings
delay <- controlDelay$delay
k <- controlDelay$k
verbose <- controlDelay$verbose
pastWeeksNotIncluded <- controlDelay$pastWeeksNotIncluded
mc.munu <- controlDelay$mc.munu
mc.y <- controlDelay$mc.y
vectorOfDates <- as.Date(salmAllOnset@epoch, origin="1970-01-01")
dayToConsider <- vectorOfDates[rangeTest[1]]
observed <- salmAllOnset@observed
population <- salmAllOnset@populationFrac
dataGLM <- bodaDelay.data.glm(dayToConsider=dayToConsider, 
                              b=b, freq=freq, 
                              epochAsDate=epochAsDate,
                              epochStr=epochStr,
                              vectorOfDates=vectorOfDates,w=w,
                              noPeriods=noPeriods,
                              observed=observed,population=population,
                              verbose=verbose,
                              pastWeeksNotIncluded=pastWeeksNotIncluded,
                              reportingTriangle=reportingTriangle, 
                              delay=delay) 
delay <- FALSE
dataGLMNoDelay <- bodaDelay.data.glm(dayToConsider=dayToConsider, 
                                  b=b, freq=freq, 
                                  epochAsDate=epochAsDate,
                                  epochStr=epochStr,
                                  vectorOfDates=vectorOfDates,w=w,
                                  noPeriods=noPeriods,
                                  observed=observed,population=population,
                                  verbose=verbose,
                                  pastWeeksNotIncluded=pastWeeksNotIncluded,
                                  reportingTriangle=reportingTriangle, 
                                  delay=delay)
test_that("the output is a data.frame",{
  expect_true(class(dataGLM)=="data.frame")
  expect_true(class(dataGLMNoDelay)=="data.frame")
})

test_that("the data frame contains all variables",{
  expect_equal(names(dataGLM)==c( "response", "wtime","population","seasgroups","vectorOfDates","delay"),rep(TRUE,6))
  expect_equal(names(dataGLMNoDelay)==c( "response", "wtime","population","seasgroups","vectorOfDates"),rep(TRUE,5)) 
  })

test_that("the variables have the right class",{
  expect_equal(class(dataGLM$response),"numeric")
  expect_equal(class(dataGLM$wtime),"numeric")
  expect_equal(class(dataGLM$population),"numeric")
  expect_equal(class(dataGLM$seasgroups),"factor")
  expect_equal(class(dataGLM$vectorOfDates),"Date")
  expect_equal(class(dataGLM$delay),"numeric")
 
  expect_equal(class(dataGLMNoDelay$response),"numeric")
  expect_equal(class(dataGLMNoDelay$wtime),"numeric")
  expect_equal(class(dataGLMNoDelay$population),"numeric")
  expect_equal(class(dataGLMNoDelay$seasgroups),"factor")
  expect_equal(class(dataGLMNoDelay$vectorOfDates),"Date")
})

test_that("the time variable is ok with diff 1",{
  delayWtime <- as.numeric(levels(as.factor(dataGLM$wtime)))
  expect_equal(diff(delayWtime)==rep(1,length(delayWtime)-1),rep(TRUE,length(delayWtime)-1))
  
  expect_equal(diff(dataGLMNoDelay$wtime)==rep(1,length(dataGLMNoDelay$wtime)-1),rep(TRUE,length(dataGLMNoDelay$wtime)-1))
})

test_that("the factor variable has the right number of levels",{
  expect_true(length(levels(dataGLM$seasgroups))==noPeriods)
  expect_true(length(levels(dataGLMNoDelay$seasgroups))==noPeriods)
})



##################################################################
context("Fit glm function")
# if (interactive() && require("INLA")) { # do not test INLA-related code on CRAN
#   ## CAVE: _R_CHECK_TIMINGS_ as queried by surveillance.options("allExamples")
#   ## is no reliable condition to skip the test on CRAN (see
#   ## https://stat.ethz.ch/pipermail/r-devel/2012-September/064812.html
#   ## ), and especially seems not to be set by the daily Windows checks.
#   argumentsGLM <- list(dataGLM=dataGLM,reportingTriangle=reportingTriangle,
#                        timeTrend=timeTrend,alpha=alpha,
#                        populationOffset=populationOffset,
#                        factorsBool=TRUE,pastAberrations=FALSE,
#                        glmWarnings=glmWarnings,
#                        verbose=verbose,delay=delay,k=k,control=controlDelay,
#                        inferenceMethod="INLA")
#   
#   model <- do.call(bodaDelay.fitGLM, args=argumentsGLM)
#   test_that("The fit glm function gives the right class of output?",{
#     expect_equal(class(model),"inla")
#   })
# }


argumentsGLM <- list(dataGLM=dataGLM,reportingTriangle=reportingTriangle,
                     timeTrend=timeTrend,alpha=alpha,
                     populationOffset=populationOffset,
                     factorsBool=TRUE,pastAberrations=FALSE,
                     glmWarnings=glmWarnings,
                     verbose=verbose,delay=delay,k=k,control=controlDelay,
                     inferenceMethod="asym")

model <- do.call(bodaDelay.fitGLM, args=argumentsGLM)
test_that("The fit glm function gives the right class of output?",{
  expect_equal(class(model)==c("negbin", "glm", "lm" ),rep(TRUE,3))
})
################################################################################
context("formula function")
################################################################################
test_that("We get the right formula",{
  expect_equal(formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE),"response ~ 1+wtime")
  expect_equal(formulaGLMDelay(timeBool=FALSE,factorsBool=FALSE),"response ~ 1")
  expect_equal(formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE),"response ~ 1+wtime")
  expect_equal(formulaGLMDelay(timeBool=TRUE,factorsBool=TRUE),"response ~ 1+wtime+as.factor(seasgroups)")
  expect_equal(formulaGLMDelay(timeBool=TRUE,factorsBool=TRUE,delay=TRUE),"response ~ 1+wtime+as.factor(seasgroups)+as.factor(delay)")
  expect_equal(formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE,outbreak=TRUE),"response ~ 1+wtime+f(outbreakOrNot,model='linear', prec.linear = 1)")  
})
