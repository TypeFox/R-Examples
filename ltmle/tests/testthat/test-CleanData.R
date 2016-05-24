context("Test CleanData")

TestCleanData <- function(seed, set.a) {
  set.seed(seed)
  n <- 1000
  ua <- rep(TRUE, n)   #ua = uncensored and alive
  L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
  W <- rnorm(n)
  C1 <- BinaryToCensoring(is.uncensored=rexpit(2 + W))
  ua <- ua & C1 == "uncensored"
  L1[ua] <- rnorm(n)[ua] + W[ua]
  A1[ua] <- rexpit(L1[ua])
  A1[ua & L1 < -2] <- 1
  A1[ua & L1 >  2] <- rbinom(n, size=1, prob=0.1)[ua & L1 >  2]
  Y1[ua] <- rexpit((W + L1 - A1)[ua])
  ua <- ua & !Y1
  C2.binary[ua] <- rexpit((1 + 0.7 * L1 - A1)[ua])
  C2 <- BinaryToCensoring(is.uncensored=C2.binary)
  ua <- ua & C2 == "uncensored"
  L2[ua] <- (0.5 * L1 - 0.9 * A1 + rnorm(n))[ua]
  A2[ua] <- rexpit((0.5 * L1 + 0.8 * L2)[ua]) | A1[ua]
  if (set.a) {
    A2[!ua] <- A1[!ua]  #set A2 to nonsense value - will be ignored
  }
  Y2[ua] <- rexpit((0.7 * L1 + L2 - 0.8 * A1 - A2)[ua])
  Y2[Y1 == 1] <- 1  # if a patient dies at time 1, record death at time 2 as well
  data <- data.frame(W, C1, L1, A1, Y1, C2, L2, A2, Y2)
  result <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
                   Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(0, 1), 
                   deterministic.g.function=NULL, survivalOutcome=TRUE)
  return(result$estimates)
}

test_that("Anodes changing after death/censoring does not change result", {
  prev.seed <- .Random.seed
  seed <- round(runif(1) * 1000)
  expect_message(r1 <- TestCleanData(seed, TRUE), regexp="Your data did not conform and has been adjusted") 
  r2 <- TestCleanData(seed, FALSE)
  .Random.seed <<- prev.seed
  expect_equal(r1, r2) 
})


