library(datafsm)
context("Main evolve_model function")

test_that("evolve_model() returns correct type of object", {
        cdata <- data.frame(period = 1:5, outcome = c(1,2,1,1,1),
                            my.decision1 = c(1,0,1,1,1), other.decision1 = c(0,0,0,1,1))
        result <- evolve_model(cdata, cv=FALSE)
        expect_is(result, "ga_fsm")
})

test_that("evolve_model() returns warnings and errors", {
  cdata <- as.matrix(data.frame(period = 1:5, outcome = c(1,2,1,1,1),
                      my.decision1 = c(1,0,1,1,1), other.decision1 = c(0,0,0,1,1)))
  expect_warning(evolve_model(cdata, cv=FALSE), "did not supply a data.frame")
  
  cdata <- data.frame(period = 1:5, outcome = c(NA,2,1,1,1),
                      my.decision1 = c(1,0,1,1,1), other.decision1 = c(0,0,0,1,1))
  expect_error(evolve_model(cdata, cv=FALSE), "missing")
  
  cdata <- data.frame(period = 1:5, outcome = c(1,2,1,1,1),
                      my.decision1 = c(1,0,1,1,1), other.decision1 = c(0,0,0,1,1),
                      joe.decision1 = c(0,0,0,1,1), jack.decision1 = c(0,0,0,1,1) )
  expect_error(evolve_model(cdata, cv=FALSE), "predictor", all=FALSE)
  
  cdata <- data.frame(period = 1:5, outcome = c(1,1,1,1,1),
                      my.decision1 = c(1,0,1,3,1), other.decision1 = c(0,0,0,1,1))
  expect_error(evolve_model(cdata, cv=FALSE), "unique", all=FALSE)
  
  cdata <- data.frame(period = 1:5,
                      my.decision1 = c(1,0,1,3,1), other.decision1 = c(0,0,0,1,1))
  expect_error(evolve_model(cdata, cv=FALSE), "predictor", all=FALSE)
  
})

