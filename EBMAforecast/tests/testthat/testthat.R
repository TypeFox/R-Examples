
data(calibrationSample)
data(testSample)


	context("Test if predictions between 0 and 1")
test_that("error for predcalibration greater 1",{
##test 1 for error predcalibration not between 0 and 1
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setPredCalibration(this.ForecastData)<-matrix(1.001,ncol=3,nrow=696) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})

test_that("error for predcalibration smaller 0",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setPredCalibration(this.ForecastData)<-matrix(-0.001,ncol=3,nrow=696) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})

test_that("error if predtest greater 1",{
##test 2 for error if predtest not between 0 and 1
#reset forecastdata
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setPredTest(this.ForecastData)<-matrix(1.001,ncol=3,nrow=348) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})
#
test_that("error if predtest smaller 0",{
##test 2 for error if predtest not between 0 and 1
#reset forecastdata
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setPredTest(this.ForecastData)<-matrix(-0.001,ncol=3,nrow=348) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})


context("Outcome set with values either 0 or 1 test")
##test 3 for error if outcomeCalibration not 0 or 1
#reset forecastdata
test_that("error if outcomeCalibration greater 1",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeCalibration(this.ForecastData)<-c(rep(1,600),rep(1.5,96)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})

test_that("error if outcomeCalibration is 0.5",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeCalibration(this.ForecastData)<-c(rep(1,600),rep(0.5,96)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})


test_that("error if outcomeCalibration smaller 0",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeCalibration(this.ForecastData)<-c(rep(1,600),rep(-0.00015,96)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})
##test 4 for error if outcomeTest not 0 or 1
#reset forecastdata
test_that("error if outcomeTest is larger 1",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeTest(this.ForecastData)<-c(rep(1,300),rep(1.5,48)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})

test_that("error if outcomeTest is 0.5",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeTest(this.ForecastData)<-c(rep(1,300),rep(0.5,48)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})
#

test_that("error if outcomeTest smaller 0",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
setOutcomeTest(this.ForecastData)<-c(rep(1,300),rep(-0.00015,48)) 
expect_that(as(this.ForecastData,"ForecastDataLogit"), throws_error())
})

context("Vector size test")
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))

###test 5 for error if length of vectors not the same
test_that("error if length of vectors are not the same",{
expect_that(setPredCalibration(this.ForecastData)<-c(rep(1,240)), throws_error()) ### too short
expect_that(setPredTest(this.ForecastData)<-c(rep(1,240)), throws_error()) ### too short
expect_that(setOutcomeCalibration(this.ForecastData)<-c(rep(1,240)), throws_error())### too short
#expect_that(setOutcomeTest(this.ForecastData)<-c(rep(1,240)), throws_error())### too short
expect_that(setPredCalibration(this.ForecastData)<-c(rep(1,940)), throws_error()) ### too long
expect_that(setPredTest(this.ForecastData)<-c(rep(1,940)), throws_error()) ### too long
expect_that(setOutcomeCalibration(this.ForecastData)<-c(rep(1,940)), throws_error())### too long
#expect_that(setOutcomeTest(this.ForecastData)<-c(rep(1,900)), throws_error())### too long
})

#### test 6 for error if columns in predCalibration and predTest differ
test_that("error if number of columns in predCalibration and predTest differ",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))

expect_that(setPredCalibration(this.ForecastData)<-matrix(1,ncol=4,nrow=696), throws_error())
expect_that(setPredCalibration(this.ForecastData)<-matrix(1,ncol=2,nrow=696), throws_error())
expect_that(setPredTest(this.ForecastData)<-matrix(1,ncol=4,nrow=348), throws_error())
expect_that(setPredTest(this.ForecastData)<-matrix(1,ncol=2,nrow=348), throws_error())
})

### test 7  check that results for calibration set and test set are the same as in paper after ensemble
context("Results Check for logit")
test_that("results are the same as presented in paper (calibration period)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", maxIter=25000, exp=3)
test_mat<-round(check1@modelWeights,2)
check_against<-(c(0.85,0.15,0.00))
expect_that(test_mat[[1]], equals(check_against[1]))
expect_that(test_mat[[2]], equals(check_against[2]))
expect_that(test_mat[[3]], equals(check_against[3]))
})


# # context("get tests")
# this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
# test_that("getPredCalibration gives PredCalibration",{
	# expect_that(getPredCalibration(this.ForecastData),equals(this.ForecastData@predCalibration))
# })


# test_that("getOutcomeCalibration gives OutcomeCalibration",{
	# expect_that(getOutcomeCalibration(this.ForecastData),equals(this.ForecastData@outcomeCalibration))
# })

# test_that("getPredTest gives predTest",{
	# expect_that(getPredTest(this.ForecastData),equals(this.ForecastData@predTest))
# })

# test_that("getOutcomeTest gives OutcomeTest",{
	# expect_that(getOutcomeTest(this.ForecastData),equals(this.ForecastData@outcomeTest))
# })

# test_that("getModelNames gives ModelNames",{
	# expect_that(getModelNames(this.ForecastData),equals(this.ForecastData@modelNames))
# })


context("set tests")
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))
test_that("setPredCalibration works",{
	setPredCalibration(this.ForecastData)<-matrix(1,ncol=3,nrow=696) 
	expect_that(this.ForecastData@predCalibration, equals(array(1,dim=c(696,3,1))))
})

test_that("setOutcomeCalibration works",{
	setOutcomeCalibration(this.ForecastData)<-rep(1,696) 
	expect_that(this.ForecastData@outcomeCalibration, equals(rep(1,696)))
})

test_that("setPredTest works",{
	setPredTest(this.ForecastData)<-matrix(1,ncol=3,nrow=348) 
	expect_that(this.ForecastData@predTest,  equals(array(1,dim=c(348,3,1))))
})

test_that("setOutcomeTest works",{
	setOutcomeTest(this.ForecastData)<-rep(1,348) 
	expect_that(this.ForecastData@outcomeTest, equals(rep(1,348)))
})

test_that("setModelNames works",{
	names<-c("Frank","Aaron","David")
	setModelNames(this.ForecastData)<-names
	expect_that(this.ForecastData@modelNames, equals(names))
})


#context("NA test")
##### test 8 check that NA's are not taken
#test_that("error if NA's are fed into ForecastData (predCalibration)",{
#expect_that(setPredCalibration(this.ForecastData)<-matrix(NA,ncol=3,nrow=696), throws_error())
#})
test_that("error if NA's are fed into ForecastData (outcomeCalibration)",{
expect_that(setOutcomeCalibration(this.ForecastData)<-c(rep(NA,696)), throws_error())
})

#test_that("error if NA's are fed into ForecastData (predTest)",{
#expect_that(setPredTest(this.ForecastData)<-matrix(NA,ncol=3,nrow=348), throws_error())
#})

test_that("error if NA's are fed into ForecastData (outcomeTest)",{
expect_that(setOutcomeTest(this.ForecastData)<-c(rep(NA,348)), throws_error())
})



context("test that makeForecastData takes arrays,matrix,and data.frame objects")
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
test_that("makeForecastData takes DF",{
calibrationSample.df<-as.data.frame(calibrationSample)
testSample.df<-as.data.frame(testSample)
this.ForecastData.df <- makeForecastData(.predCalibration=calibrationSample.df[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample.df[,"Insurgency"],.predTest=testSample.df[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample.df[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
expect_that(this.ForecastData,equals(this.ForecastData.df))
})

test_that("makeForecastData takes matrix",{
calibrationSample.m<-as.matrix(calibrationSample)
testSample.m<-as.matrix(testSample)
this.ForecastData.m <- makeForecastData(.predCalibration=calibrationSample.m[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample.m[,"Insurgency"],.predTest=testSample.m[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample.m[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
expect_that(this.ForecastData,equals(this.ForecastData.m))
})

test_that("makeForecastData takes array",{
calibrationSample.a<-as.array(calibrationSample)
testSample.a<-as.array(testSample)
this.ForecastData.a <- makeForecastData(.predCalibration=calibrationSample.a[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample.a[,"Insurgency"],.predTest=testSample.a[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample.a[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
expect_that(this.ForecastData,equals(this.ForecastData.a))
})



context("test for functionality of options in logit EBMA (logit)")
test_that("tolerance changes if option is used (logit)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.00141, maxIter=25000, exp=3)		
expect_that(check1@tol,equals(0.00141))	
})

test_that("maximum iteration changes if option is used (logit)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.0000000001, maxIter=15, exp=3)		
expect_that(check1@maxIter,equals(15))
})

test_that("exponent changes if option is used (logit)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.0001, maxIter=25000, exp=15)		
expect_that(check1@exp,equals(15))
})

test_that("model parameters are turned of, all parameters are 0,1 (logit)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, useModelParams=FALSE)
parameters<-matrix(c(0,1,0,1,0,1),ncol=3)		
for(i in 1:2){
	for(j in 1:3){
		expect_that(matrix(check1@modelParams,ncol=3)[i,j], equals(parameters[i,j]))
			}
}
})


context("test for functionality of options in logit EBMA by checking if results are different")
test_that("tolerance changes if option is used (logit - results)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.000141, maxIter=25000, exp=3)
check2<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.1, maxIter=25000, exp=3)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	
})

test_that("maximum iteration changes if option is used (logit - results)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=3, exp=3)
check2<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, exp=3)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	
})

test_that("exponent changes if option is used (logit - results)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, exp=1)
check2<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, exp=25)		
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	
})

test_that("model parameters are turned of, all parameters are 0,1 (logit - results)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, useModelParams=TRUE)
check2<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.001, maxIter=25000, useModelParams=FALSE)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	
})

test_that("If different initial weights are used, iterations should increase (logit)",{
  this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
  check1<-calibrateEnsemble(this.ForecastData, model="logit", maxIter=25000, exp=1)
  check2<-calibrateEnsemble(this.ForecastData, model="logit", maxIter=25000, exp=1,W=c(0.01,0.19,0.8))
  expect_false((check1@iter==check2@iter)[[1]])  
})

test_that("model option = normal changes results (logit - results)",{
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"], .modelNames=c("LMER", "SAE", "GLM"))
check1<-calibrateEnsemble(this.ForecastData, model="logit", tol=0.01, maxIter=25000, exp=3,useModelPara=FALSE)
check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.01, maxIter=25000, exp=3,useModelPara=FALSE)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])
})


context("test for functionality of options in normal EBMA")
#create data frame
set.seed(123)
predictions<-matrix(NA, nrow=400, ncol=4)
predictions[,1]<-rnorm(400,mean=2.6,sd=5)
predictions[,2]<-rnorm(400,mean=6,sd=10)
predictions[,3]<-rnorm(400,mean=0.4,sd=8)
predictions[,4]<-rnorm(400,mean=-2,sd=15)
true<-rep(NA,400)
true<-rnorm(400,mean=2.2,sd=2)

test.pred<-matrix(NA, nrow=40, ncol=4)
test.pred[,1]<-rnorm(40,mean=2.3,sd=7)
test.pred[,2]<-rnorm(40,mean=3.3,sd=12)
test.pred[,3]<-rnorm(40,mean=1.3,sd=11)
test.pred[,4]<-rnorm(40,mean=2.2,sd=18)
test.true<-rnorm(40,mean=2.2,sd=2)


test_that("tolerance changes if option is used (normal)",{
this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.000141, maxIter=25000, exp=3)
expect_that(check1@tol,equals(0.000141))	
})

test_that("maximum iteration changes if option is used (normal)",{
this.ForecastData <-makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check111<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.0000000001, maxIter=25, exp=3)		
expect_that(check111@maxIter,equals(25))
})

test_that("exponent changes if option is used (normal)",{
this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))	
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.0001, maxIter=25000, exp=15)		
expect_that(check1@exp,equals(15))
})

test_that("If different initial weights are used, iterations should increase (normal)",{
    this.ForecastData <-makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
    check1<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=1)
    check2<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=1,W=c(0.8,0.1,0.05,0.05))
    expect_false((check1@iter==check2@iter)[[1]])  
})

test_that("model parameters are turned of, all parameters are 0,1 (normal)",{
this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.0001, maxIter=25000, useModelParams=FALSE)
parameters<-matrix(c(0,1,0,1,0,1,0,1),ncol=4)	
for(i in 1:2){
	for(j in 1:4){
		expect_that(matrix(check1@modelParams,ncol=4)[i,j], equals(parameters[i,j]))
			}
}
})


### same test check if results change

context("test for functionality of options in normal EBMA, look for different results (normal - results)")
#create data frame
set.seed(123)
predictions<-matrix(NA, nrow=400, ncol=4)
predictions[,1]<-rnorm(400,mean=2.6,sd=5)
predictions[,2]<-rnorm(400,mean=6,sd=10)
predictions[,3]<-rnorm(400,mean=0.4,sd=8)
predictions[,4]<-rnorm(400,mean=-2,sd=15)
true<-rep(NA,400)
true<-rnorm(400,mean=2.2,sd=2)

test.pred<-matrix(NA, nrow=40, ncol=4)
test.pred[,1]<-rnorm(40,mean=2.3,sd=7)
test.pred[,2]<-rnorm(40,mean=3.3,sd=12)
test.pred[,3]<-rnorm(40,mean=1.3,sd=11)
test.pred[,4]<-rnorm(40,mean=2.2,sd=18)
test.true<-rnorm(40,mean=2.2,sd=2)

test_that("tolerance changes if option is used (normal - results)",{
this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.000141, maxIter=25000, exp=3)
check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=1, maxIter=25000, exp=3)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	
})

test_that("maximum iteration changes if option is used (normal - results)",{
this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))	
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.01, maxIter=1, exp=3)
check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.01, maxIter=25000, exp=3)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])	

})

#test_that("exponent changes if option is used (normal - results)",{
#	this.ForecastData <- #makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m#1", "m2", "m3","m4"))
#check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.000001, maxIter=25000, exp=1)
#check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.000001, maxIter=25000, exp=70)
#expect_false((check1@modelWeights==check2@modelWeights)[[1]])
#})

test_that("model parameters are turned of, all parameters are 0,1 (normal - results)",{
	this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.00001, maxIter=25000, exp=3,useModelPara=FALSE)
check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.00001, maxIter=25000, exp=3,useModelPara=TRUE)
expect_false((check1@modelWeights==check2@modelWeights)[[1]])
})


test_that("predType changes prediction (normal - results)",{
	this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
check1<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.00001, maxIter=25000, exp=3,useModelPara=FALSE,predType="posteriorMedian")
check2<-calibrateEnsemble(this.ForecastData, model="normal", tol=0.00001, maxIter=25000, exp=3,useModelPara=FALSE,predType="posteriorMean")
expect_true((check1@modelWeights==check2@modelWeights)[[1]])
expect_false((check1@predTest[,1,1]==check2@predTest[,1,1])[[1]])
})

test_that("model option = logit changes results (normal - results)",{
		this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
expect_error(calibrateEnsemble(this.ForecastData, model="logit", tol=0.01, maxIter=25000, exp=3,useModelPara=FALSE))
})

context("test that results are same as in Raftery package")
#create data frame
data(presidentialForecast)
tyn=15
a=1
train.years=14
dates <- rep(NA, tyn)
   for (i in 1:tyn){
     dates[i] <- paste("2011", "01", 10+i, "01", sep="")
    }

   pred.date <- dates[tyn]
full.forecasts<-presidentialForecast[,c(1:6)]
full.observed<-presidentialForecast[,7]
library(ensembleBMA)
test_that("same result as in Raftery",{
   my.E.data <- ensembleData(forecasts=(full.forecasts)^(1/a), dates=dates, observations=full.observed,
                             initializationTime=1, forecastHour=1) #Make a dataset of the appropriate format for the ensembleBMA package
   fit.eBMA <- ensembleBMAnormal(my.E.data, trainingDays=train.years, dates=pred.date, minCRPS=FALSE,
                              control=controlBMAnormal(biasCorrection="none",tol=0.000000001))
my.data<-makeForecastData(.predCalibration=full.forecasts[c(1:14),],.outcomeCalibration=full.observed[c(1:14)],.predTest=full.forecasts[15,],.outcomeTest=full.observed[15], c("Campbell", "Lewis-Beck","EWT2C2","Fair","Hibbs","Abramowitz"))
check1<-calibrateEnsemble(my.data, model="normal", maxIter=25000,useModelPara=FALSE,tol=0.000000001)
round(check1@modelWeights,4)                            
## this needs to be fixed
round(fit.eBMA$weights,4)
check2<-as.numeric(round(as.matrix(check1@modelWeights)[1:5,],3))
expect_that(as.numeric(round(as.matrix(fit.eBMA$weights)[1:5,],3)),equals(check2))
})



context("test that results are same as in Raftery package with missing obs")
data(presidentialForecast)
tyn=15
a=1
train.years=14

dates <- rep(NA, tyn)
   for (i in 1:tyn){
     dates[i] <- paste("2011", "01", 10+i, "01", sep="")
    }

   pred.date <- dates[tyn]
full.forecasts<-presidentialForecast[,c(1:6)]
full.observed<-presidentialForecast[,7]
full.forecasts[1,6]<-NA
full.forecasts[3,2]<-NA
full.forecasts[2,2]<-NA
full.forecasts[7,2]<-NA
full.forecasts[6,1]<-NA
full.forecasts[14,2]<-NA
full.forecasts[7,6]<-NA
library(ensembleBMA)
test_that("same result as in Raftery",{
   my.E.data <- ensembleData(forecasts=(full.forecasts)^(1/1), dates=dates, observations=full.observed,
                             initializationTime=1, forecastHour=1) #Make a dataset of the appropriate format for the ensembleBMA package
   fit.eBMA <- ensembleBMAnormal(my.E.data, trainingDays=train.years, dates=pred.date, minCRPS=FALSE,
                              control=controlBMAnormal(biasCorrection="none",tol=0.00000001))
my.data<-makeForecastData(.predCalibration=full.forecasts[c(1:14),],.outcomeCalibration=full.observed[c(1:14)],.predTest=full.forecasts[15,],.outcomeTest=full.observed[15], c("Campbell", "Lewis-Beck","EWT2C2","Fair","Hibbs","Abramowitz"))
check13<-calibrateEnsemble(my.data, model="normal", maxIter=25000,useModelPara=FALSE,tol=0.00000001)
## this needs to be fixed
check2<-as.numeric(round(as.matrix(check13@modelWeights),3))
expect_that(as.numeric(round(as.matrix(fit.eBMA$weights),3)),equals(check2))
})

context("Make sure demos actually run")
test_that("logit EBMA model",{
  demo(EBMAforecast)})
  
test_that("logit EBMA model, with missing obs",{
  data(calibrationSample)
  data(testSample)
  missing = sample(1:nrow(calibrationSample), 300,replace=FALSE)
  calibrationSample[missing[1:100],c("LMER")] = NA
  calibrationSample[missing[101:200],c("SAE")] = NA
  calibrationSample[missing[201:300],c("GLM")] = NA


this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
                                                                 .outcomeCalibration=calibrationSample[,"Insurgency"],
                                                                .predTest=testSample[,c("LMER", "SAE", "GLM")],
                                                                 .outcomeTest=testSample[,"Insurgency"],
                                                                 .modelNames=c("LMER", "SAE", "GLM"))
  
this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.0001, maxIter=25000, exp=3)
})



test_that("logit EBMA model",{
  demo(presForecast)})

test_that("normal EBMA model, with missing obs",{
  data(presidentialForecast)
  missing = sample(1:nrow(presidentialForecast), 12,replace=FALSE)
  full.forecasts<-presidentialForecast[,c(1:6)]
  full.forecasts[missing[c(1,2)],1] = NA
  full.forecasts[missing[c(3,4)],1] = NA
  full.forecasts[missing[c(5,6)],1] = NA
  full.forecasts[missing[c(7,8)],1] = NA
  full.forecasts[missing[c(9,10)],1] = NA
  full.forecasts[missing[c(11,12)],1] = NA
  full.observed<-presidentialForecast[,7]
  
  this.ForecastData<-makeForecastData(.predCalibration=full.forecasts[c(1:14),],.outcomeCalibration=full.observed[c(1:14)],.predTest=full.forecasts[15,], .outcomeTest=full.observed[15], .modelNames=c("Campbell", "Lewis-Beck","EWT2C2","Fair","Hibbs","Abramowitz"))
  thisEnsemble<-calibrateEnsemble(this.ForecastData, model="normal", useModelParams=FALSE, tol=0.000000001)
  })


context("Test function to create predictions without reestimating model")
set.seed(123)
predictions<-matrix(NA, nrow=400, ncol=4)
predictions[,1]<-rnorm(400,mean=2.6,sd=5)
predictions[,2]<-rnorm(400,mean=6,sd=10)
predictions[,3]<-rnorm(400,mean=0.4,sd=8)
predictions[,4]<-rnorm(400,mean=-2,sd=15)
true<-rep(NA,400)
true<-rnorm(400,mean=2.2,sd=2)

test.pred<-matrix(NA, nrow=40, ncol=4)
test.pred[,1]<-rnorm(40,mean=2.3,sd=7)
test.pred[,2]<-rnorm(40,mean=3.3,sd=12)
test.pred[,3]<-rnorm(40,mean=1.3,sd=11)
test.pred[,4]<-rnorm(40,mean=2.2,sd=18)
test.true<-rnorm(40,mean=2.2,sd=2)


new.pred = matrix(rnorm(80,5,7),ncol=4)

test_that("EBMApredict for normal EBMA model,model parameters TRUE",{
      this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
      check1<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=3)
      newPred = EBMApredict(check1,new.pred)
      test_pred2 = EBMApredict(check1,test.pred)
      expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
})

test_that("EBMApredict for normal EBMA model,model parameters TRUE, pedictionType mean",{
  this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
  check1<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=1,useModelParams=TRUE, predType="posteriorMean")
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,test.pred)
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
  
})

test_that("EBMApredict for normal EBMA model,model parameters FALSE",{
  this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
  check1<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=1,useModelParams=FALSE)
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,test.pred)
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
})
  

test_that("EBMApredict for normal EBMA model,model parameters FALSE, pedictionType mean",{
  this.ForecastData <- makeForecastData(.predCalibration=predictions,.outcomeCalibration=true,.predTest=test.pred,.outcomeTest=test.true, .modelNames=c("m1", "m2", "m3","m4"))
  check1<-calibrateEnsemble(this.ForecastData, model="normal", maxIter=25000, exp=1,useModelParams=FALSE, predType="posteriorMean")
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,test.pred)
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
})


data(calibrationSample)
data(testSample)
new.pred = matrix(runif(60,0.02,0.98),ncol=3)
this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],.outcomeCalibration=calibrationSample[,"Insurgency"],.predTest=testSample[,c("LMER", "SAE", "GLM")],.outcomeTest=testSample[,"Insurgency"],.modelNames=c("LMER", "SAE", "GLM"))

test_that("EBMApredict for logit EBMA model,model parameters FALSE",{
  check1 <- calibrateEnsemble(this.ForecastData,model="logit",exp=3,useModelParams=FALSE)
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,this.ForecastData@predTest,Outcome=this.ForecastData@outcomeTest)
  test_pred2 = EBMApredict(check1,this.ForecastData@predTest)
  
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
  
})



test_that("EBMApredict for logit EBMA model,model parameters TRUE",{
  check1 <- calibrateEnsemble(this.ForecastData,model="logit",exp=3,useModelParams=TRUE)
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,this.ForecastData@predTest)
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
})


test_that("EBMApredict for logit EBMA model,model parameters TRUE",{
  check1 <- calibrateEnsemble(this.ForecastData,model="logit",exp=3,useModelParams=TRUE)
  newPred = EBMApredict(check1,new.pred)
  test_pred2 = EBMApredict(check1,this.ForecastData@predTest,Outcome=this.ForecastData@outcomeTest)
  expect_that(as.numeric(test_pred2@predTest[,1,]),equals(as.numeric(check1@predTest[,1,])))
})
