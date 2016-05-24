context("Test print.summary.ltmleMSM") 

test_that("transformOutcome NOTE is printed by print.summary.ltmleMSM", {

  data(sampleDataForLtmleMSM)
  Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
  Lnodes <- c("CD4_1", "CD4_2")
  Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))

  data <- sampleDataForLtmleMSM$data
  data[,Ynodes] <- matrix(runif(length(Ynodes)*nrow(data)), nrow(data))*10
  data <- as.data.frame(apply(data, 2, function(x) {
    x[is.na(x)] <- max(x, na.rm=TRUE)
    x
    }))

  result <- ltmleMSM(data, Anodes=Anodes, Lnodes=Lnodes,
             Ynodes=Ynodes, survivalOutcome=FALSE,
                     regimes=sampleDataForLtmleMSM$regimes, 
                     summary.measures=sampleDataForLtmleMSM$summary.measures,
                     final.Ynodes=Ynodes, 
                     working.msm="Y ~ time + I(pmax(time - switch.time, 0))", 
                     estimate.time=FALSE, IC.variance.only=TRUE)

  expect_that(print(result), prints_text("NOTE: The MSM is modeling the transformed outcome"))
  expect_that(print(summary(result)), prints_text("NOTE: The MSM is modeling the transformed outcome"))    

})