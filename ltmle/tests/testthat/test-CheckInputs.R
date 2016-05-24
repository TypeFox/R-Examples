context("Test CheckInputs")

test_that("Lnodes, Anodes, Cnodes, Ynodes include all columns after the first one", {
  n <- 10 
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = NULL, 
     Ynodes = "Y", 
     abar = c(1,1),
     survivalOutcome = FALSE,
     estimate.time = FALSE),
     throws_error("All nodes after the first of A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes"))
  })

test_that("In formulas RHS variables are parents of the current node", {
  n <- 10 
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )  

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~W", "A2~Y+L2+A1+W"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       throws_error("Some nodes in gform\\[2\\] are not parents of A2"))

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    Qform = c(L2="Q.kplus1~A1+W", Y="Q.kplus1~L2+A2+Y"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       throws_error("Some nodes in Qform\\[2\\] are not parents of Y"))

  #No baseline covars before first A node
  dat2 <- dat[, -1]
  expect_that(ltmle(dat2, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~A1", "A2~1"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),
       throws_error("Some nodes in gform\\[1\\] are not parents of A1")) 

  expect_that(ltmle(dat2, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~1", "A2~1"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       is_a("ltmle"))            
})

test_that("Non-binary outcome can be in [0, 1]", {
  n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, IC.variance.only=TRUE, estimate.time = FALSE), is_a("ltmle"))
})

test_that("survivalOutcome=TRUE requires binary outcomes.", {
  n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, survivalOutcome=TRUE, IC.variance.only=TRUE, estimate.time = FALSE), 
    throws_error("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA"))
})

test_that("binaryOutcome flag set correctly", {
  n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, survivalOutcome=FALSE, IC.variance.only=TRUE,  estimate.time = FALSE)$binaryOutcome, is_false())
  expect_that(ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=1, survivalOutcome=FALSE)$binaryOutcome, is_true())  
})


test_that("Y outside of [0, 1] is scaled", {
  n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = runif(n)*10)

  expect_that(l <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, IC.variance.only=TRUE,  estimate.time = FALSE), 
    shows_message("Some Ynodes are not in \\[0, 1\\], and Yrange was NULL, so all Y nodes are being\\ntransformed to \\(Y-min.of.all.Ys\\)/range.of.all.Ys"))  
  expect_that(l$transformOutcome, is_true())
  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, Yrange=c(0, 3), IC.variance.only=TRUE, estimate.time = FALSE), 
    shows_message("Some Ynodes are not in \\[Yrange\\[1\\], Yrange\\[2\\]\\], Y values are truncated"))
  expect_that(ltmle(transform(data, Y=Y/10), Anodes="A", Ynodes="Y", abar=1, IC.variance.only=TRUE, estimate.time = FALSE)$transformOutcome, is_false())
})

test_that("survivalOutcome required if outcome is binary", {
  n <- 10
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y1 = rbinom(n, 1, .5),
                     Y2 = rbinom(n, 1, .5))

  expect_that(ltmle(data, Anodes="A", Ynodes=c("Y1", "Y2"), abar=1, estimate.time = FALSE),
    throws_error("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified."))
})
