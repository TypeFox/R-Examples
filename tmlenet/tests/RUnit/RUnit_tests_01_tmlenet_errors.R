# ---------------------------------------------------------------------------------
# TEST SET 1
# SIMPLE NETWORK TMLE
# ---------------------------------------------------------------------------------

test.tmleneterrrors <- function() {
  # Set x% of community to A=1 (returns probability P(A=1))
  f.A_x <- function(data, x, ...) rbinom(n = nrow(data), size = 1, prob = x[1])
  # Deterministically set every A=0
  f.A_0 <- function(data, ...) f.A_x(data, 0, ...)
  # Deterministically set every A=1
  f.A_1 <- function(data, ...) f.A_x(data, 1, ...)

  # incorrect intervention f_gstar1:
  f.A_wrong <- function(x, ...) 1

  #***************************************************************************************
  # EXAMPLE WITH SIMULATED DATA FOR 6 FRIENDS AND 3 W's (OLD SIMULATION 3)
  #***************************************************************************************
  data(df_netKmax6) # Load the network data
  Kmax <- 6 # Max number of friends in the network

  #----------------------------------------------------------------------------------
  # Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
  # Intercept based TMLE
  #----------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)

  def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
            def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
            def.sA(netA = A[[0:Kmax]])

  # Test for non-existing predictors in hform.g0/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + netW3_sum + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for non-existing predictors in hform.g0/hform.gstar:
  checkException(
    res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ netW3_sum",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for non-existing outcomes in hform.g0/hform.gstar:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "sum.netW3 ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for different outcomes in hform.g0/hform.gstar (non-existing for hform.gstar):
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netW3 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Test for existing but different outcomes in hform.g0/hform.gstar:
  checkException(
       res_K6_1 <- tmlenet(data = df_netKmax6,
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "sum.netAW2 ~ sum.netW3 + nF",
                    Qform = "Y ~ sum.netW3 + sum.netAW2",

                    Anode = "A", Ynode = "Y",
                    Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                    f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Throws exception since netW3_sum, sum_1mAW2_nets from Qform don't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ netW3_sum + sum_1mAW2_nets",
                  hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",

                  Anode = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))
  )
  # Throws exception since Ynode arg is omitted, but blah from Qform LHS doesn't exist:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Qform = " blah ~ sum.netW3 + sum.netAW2",
                  hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                  hform.gstar = "netA ~ sum.netW3",
                  Anode = "A",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_0, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))

  # Throws exception when f_gstar1 function doesn't have argument "data":
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Anode = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = f.A_wrong, sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))
  # Throws an exception when f_gstar1 is a vector of 1<length<n:
  checkException(
     res_K6_1 <- tmlenet(data = df_netKmax6,
                  Anode = "A", Ynode = "Y",
                  Kmax = Kmax, IDnode = "IDs", NETIDnode = "Net_str",
                  f_gstar1 = rep(1L,50), sW = def_sW, sA = def_sA, optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1)))




}






