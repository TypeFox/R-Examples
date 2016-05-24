# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ----------------------------------------------------------------------------#
# Test whether empirical and and theoretical rejection probability coincide   #                                    
# --------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

context("Emprical vs theoretical rejection probability for the class selBias")

test_that("the empirical and theoretical frequency coincide for the class selbias", {
  # Number of patients
  N <- 4
  # Number of generated sequences
  r <- 10000  
  # = 1 - confidence level for simulations
  alpha <- 0.01
  # type-I-error probability of the t.test
  toe <- runif(1, min = 0.1, max = 0.15)  
  # variance of the endpoint
  sigma <- runif(1, min = 0.5, max = 5)
  # magnitude of selection bias
  eta <- runif(1, max = 4)
  # type of selection bias
  type <- sample(c("CS", "DS"), 1)
  endp <- normEndp(mu = c(runif(1, max = 3), runif(1, max = 3)), sigma = c(sigma, sigma))
  selBiasEx <- selBias(type = type, eta = eta, method = "exact", alpha = toe)
  selBiasSim <- selBias(type = type, eta = eta, method = "sim", alpha = toe)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Test for the average rejection probability of selection bias for cr
  
  # true rejection probability
  rejProbEx <- summary(assess(getAllSeq(crPar(N)), selBiasEx, endp = endp))[1,1]
  # simualted rejection probability
  rejProbSim <- summary(assess(genSeq(crPar(N), r), selBiasSim, endp = endp))[1,1]
  tol <- qnorm(1-alpha/2)*sqrt(rejProbSim*(1-rejProbSim)/r)  

  expect_true(abs(rejProbEx - rejProbSim) < tol)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 2. Test for the average rejection probability of selection bias for a given sequence
  
  # true rejection probability of a sequence
  RS <- getAllSeq(crPar(N))
  # sample sequence
  RRS <- sample(nrow(RS@M), 1)
  # true rejection probability of a sequence
  rejProbEx <- assess(RS, selBiasEx, endp = endp)@D[RRS, 3]
  # simulated rejection probability of a sequence
  simRS <- genSeq(crPar(N), r)
  simRS@M <- matrix(rep(RS@M[RRS, ], r), nrow = r, byrow = T)
  rejProbSim <- summary(assess(simRS, selBiasSim, endp = endp))[1, 1]
  tol <- qnorm(1-alpha/2)*sqrt(rejProbSim*(1-rejProbSim)/r)  
  
  expect_true(abs(rejProbEx - rejProbSim) < tol)

  }
)
  