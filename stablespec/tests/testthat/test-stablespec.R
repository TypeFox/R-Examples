context("Expected input arguments and output in stableSpec and plotStability")

test_that("Incorrect/missing argument data yields errors", {

#------argument data
  expect_error(stableSpec(theData=NULL, nSubset=25,
                       iteration=20, nPop=100,
                       mutRate=0.075, crossRate=0.85,
                       longitudinal=FALSE, numTime=1, seed=NULL,
                       co="covariance", consMatrix=matrix(1:4, 2, 2),
                       threshold=0.6, toPlot=FALSE),
               "Data cannot be missing")

  expect_error(stableSpec(theData=1:10, nSubset=25,
                          iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Data should be either a data frame or a matrix of numerical values.")

  expect_error(stableSpec(theData=c("a", "b"), nSubset=25,
                          iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Data should be either a data frame or a matrix of numerical values.")

  expect_error(stableSpec(theData=data.frame(letter=letters[1:3], number=1:3),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Data should be either a data frame or a matrix of numerical values.")
})

test_that("Incorrect/missing numerical or logical or character arguments yields errors", {
#-------numerical/logical/character arguments

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset="notNumerical", iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument nSubset should be positive numeric, e.g., 10.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration="notNumerical", nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument iteration or NSGA-II generations should be positive numeric, e.g., 20.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop="notNumerical",
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument nPop should be positive numeric, e.g., 50.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate="notNumerical", crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument mutRate should be positive numeric, e.g., 0.075.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate="notNumerical",
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument crossRate should be positive numeric, e.g., 0.85.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime="notNumerical", seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument numTime should be positive numeric, e.g., 1 if cross-sectional data.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold="notNumerical", toPlot=FALSE),
               "Argument threshold should be positive numeric, e.g., 0.6.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=1:3, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument longitudinal should be either logical TRUE or FALSE.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=NULL, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument longitudinal cannot be missing.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed="notNumerical",
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument seed should be numeric vector.")



  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co=20, consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument co should be a vector of characters, e.g., either covariance or correlation.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="wrongvector", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Argument co should be either covariance or correlation matrix.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=1:3,
                          threshold=0.6, toPlot=FALSE),
               "The constraints should be formed in a matrix.")


  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot="notLogical"),
               "Argument toPlot should be either logical TRUE or FALSE.")

})


#----possible conflicting arguments
test_that("Conflicting arguments yield errors.", {

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=TRUE, numTime=1, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Longitudinal data should have more than one time slices, e.g., numTime = 2, with two time slices.")

  expect_error(stableSpec(theData=data.frame(number1=1:3, number2=4:6, number3=7:9, number4=10:12),
                          nSubset=25, iteration=20, nPop=100,
                          mutRate=0.075, crossRate=0.85,
                          longitudinal=FALSE, numTime=2, seed=NULL,
                          co="covariance", consMatrix=matrix(1:4, 2, 2),
                          threshold=0.6, toPlot=FALSE),
               "Cross-sectional data should have only one time slice, e.g., numTime = 1")
})




test_that("Incorrect/missing argument in plotStability yields errors", {

  skip_on_cran()
  library(graph)


  result <- stableSpec(theData=adhd, nSubset=1, iteration=5, nPop=10,
                       mutRate=0.075, crossRate=0.85, longitudinal=FALSE, numTime=1,
                       seed=NULL, co="covariance",
                       consMatrix=matrix(c(2, 1, 3, 1, 4, 1, 5, 1, 6, 1), 5, 2, byrow=TRUE),
                       threshold=0.1, toPlot=FALSE)

  expect_error(plotStability(listOfFronts=1:3, threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument listOfFronts should be a list.")

  expect_error(plotStability(listOfFronts=NULL, threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument listOfFronts cannot be missing.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold="a",
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument threshold should be positive numeric, e.g., 0.6.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=1:3,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument stableCausal should be a list.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=NULL,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument stableCausal cannot be missing.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=1:3,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument stableCausal_l1 should be a list.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=NULL,
                             stableEdge=result$edgeStab,longitudinal=FALSE),
               "Argument stableCausal_l1 cannot be missing.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=1:3,longitudinal=FALSE),
               "Argument stableEdge should be a list.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=NULL,longitudinal=FALSE),
               "Argument stableEdge cannot be missing")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,
                             longitudinal=1:3),
               "Argument longitudinal should be either logical TRUE or FALSE.")

  expect_error(plotStability(listOfFronts=result$listOfFronts,
                             threshold=0.6,
                             stableCausal=result$causalStab,
                             stableCausal_l1=result$causalStab_l1,
                             stableEdge=result$edgeStab,
                             longitudinal=NULL),
               "Argument longitudinal cannot be missing.")


  # outputs

  expect_true(is.list(result$listOfFronts))
  expect_true(is.list(result$causalStab))
  expect_true(is.list(result$causalStab_l1))
  expect_true(is.list(result$edgeStab))
  expect_true(is.vector(result$allSeed))
  expect_true(is.numeric(result$allSeed))
  expect_true(is.matrix(result$relCausalPath))
  expect_true(is.matrix(result$relCausalPath_l1))
  expect_true(is.matrix(result$relEdge))
  expect_is(result$graph, "graphNEL")


})






