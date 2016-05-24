# ---------------------------------------------------------------------------------
# TEST SET 1
# SIMPLE NETWORK TMLE
# ---------------------------------------------------------------------------------

test.examples <- function() {
  # Stochastically sets (100*x)% of community to A=1
  # Returns a function that will sample A with probability x:=P(A=1))
  make_f.gstar <- function(x, ...) {
    eval(x)
    f.A_x <- function(data, ...){
      rbinom(n = nrow(data), size = 1, prob = x[1])
    }
    return(f.A_x)
  }
  # Deterministic f_gstar setting every A=0:
  f.A_0 <- make_f.gstar(x = 0)
  # Deterministic f_gstar setting every A=1:
  f.A_1 <- make_f.gstar(x = 1)
  # Stochastic f_gstar that sets A=1 with probability 0.2:
  f.A_.2 <- make_f.gstar(x = 0.2)

  #***************************************************************************************
  # EXAMPLE WITH SIMULATED DATA FOR 6 FRIENDS AND 3 W's (OLD SIMULATION 3)
  #***************************************************************************************
  tmlenet:::checkpkgs(pkgs = c("stringr"))
  require(stringr)

  # Max number of friends in the network:
  Kmax <- 6
  # Load simulation function:
  source("./datgen_nets/sim3_datgen_k6.R")  # to load from inside run-time test dir
  # source("../datgen_nets/sim3_datgen_k6.R") # to load from current file dir
  # Simulate network data:
  # set.seed(543)
  n <- 1000
  df_netKmax6 <- gendata_pop(nC = 1, n_arr = n, k_arr = Kmax, EC_arr = EC, f.g_list = "f.A", f.g_args_list = list(NULL), rndseed = 543)
  # save(df_netKmax6, file = "./df_netKmax6.rda")

  print(head(df_netKmax6))
  #   IDs W1 W2 W3 A Y nFriends                  Net_str
  # 1  I1  3  0  1 1 1        1                     I537
  # 2  I2  3  1  0 0 1        5    I6 I58 I595 I641 I654
  # 3  I3  5  0  0 0 0        5 I163 I637 I650 I722 I783
  # 4  I4  2  1  0 1 1        2                 I49 I995
  # 5  I5  3  1  1 0 1        3           I358 I369 I762
  # 6  I6  2  0  1 1 1        2                I682 I917
  print(class(df_netKmax6$A)) # [1] "integer"
  print(class(df_netKmax6$nFriends)) # [1] "numeric"
  print(table(df_netKmax6$W1))
   #  0   1   2   3   4   5
   # 50 170 302 242 180  56
  print(c(mean(df_netKmax6$W1), mean(df_netKmax6$W2), mean(df_netKmax6$W3)))
  # [1] 2.500 0.553 0.589
  print(mean(df_netKmax6$A)) # [1] 0.198
  print(mean(df_netKmax6$Y)) # [1] 0.435
  print(mean(df_netKmax6$nFriends)) # [1] 3.307

  #----------------------------------------------------------------------------------
  # Example 1. Mean population outcome under deterministic intervention A=0 with 6 friends
  # Intercept based TMLE
  #----------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)
  def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
            def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
            def.sA(netA = A[[0:Kmax]])

  #***************************************************************************************
  # # correct version(s):
  #***************************************************************************************
  # No Ynode:
 res_K6_1a <- tmlenet(data = df_netKmax6, Kmax = Kmax, Anode = "A", f_gstar1 = f.A_0, sW = def_sW, sA = def_sA,
                      Qform = "Y ~ sum.netW3 + sum.netAW2",
                      hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                      hform.gstar = "netA ~ sum.netW3",
                      IDnode = "IDs", NETIDnode = "Net_str",
                      optPars = list(n_MCsims = 1))

  tmle_idx <- rownames(res_K6_1a$EY_gstar1$estimates)%in%"tmle"
  h_iptw_idx <- rownames(res_K6_1a$EY_gstar1$estimates)%in%"h_iptw"
  gcomp_idx <- rownames(res_K6_1a$EY_gstar1$estimates)%in%"gcomp"

  # Test estimates:
  checkTrue(abs(res_K6_1a$EY_gstar1$estimates[tmle_idx] - 0.5051903) < 10^(-06))
  checkTrue(abs(res_K6_1a$EY_gstar1$estimates[h_iptw_idx] - 0.5065960) < 10^(-06))
  checkTrue(abs(res_K6_1a$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_1a$EY_gstar1$vars[tmle_idx] - 0.0009268804) < 10^(-06))
  checkTrue(abs(res_K6_1a$EY_gstar1$vars[h_iptw_idx] - 0.0021023317) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_1a$EY_gstar1$CIs[tmle_idx][1] - 0.4455197) < 10^(-06)) &
              (abs(res_K6_1a$EY_gstar1$CIs[tmle_idx][2] - 0.5648608) < 10^(-06)))
  checkTrue((abs(res_K6_1a$EY_gstar1$CIs[h_iptw_idx][1] - 0.4167293) < 10^(-06)) &
              (abs(res_K6_1a$EY_gstar1$CIs[h_iptw_idx][2] - 0.5964627) < 10^(-06)))
  # res_K6_1a$EY_gstar1$other.vars

  #----------------------------------------------------------------------------------
  # Example 2. Same as above but for covariate-based TMLE
  #----------------------------------------------------------------------------------
  res_K6_2 <- tmlenet(data = df_netKmax6, Kmax = Kmax, Anode = "A", Ynode = "Y", f_gstar1 = f.A_0, sW = def_sW, sA = def_sA,
                      Qform = "Y ~ sum.netW3 + sum.netAW2",
                      hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                      hform.gstar = "netA ~ sum.netW3",
                      IDnode = "IDs", NETIDnode = "Net_str",
                      optPars = list(runTMLE = "tmle.covariate",n_MCsims = 1))

  # Test estimates:
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[tmle_idx] - 0.5053725) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[h_iptw_idx] - 0.5065960) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_2$EY_gstar1$vars[tmle_idx] - 0.0009265246) < 10^(-06))
  checkTrue(abs(res_K6_2$EY_gstar1$vars[h_iptw_idx] - 0.0021023317) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_2$EY_gstar1$CIs[tmle_idx][1] - 0.4457134) < 10^(-06)) &
              (abs(res_K6_2$EY_gstar1$CIs[tmle_idx][2] - 0.5650315) < 10^(-06)))
  checkTrue((abs(res_K6_2$EY_gstar1$CIs[h_iptw_idx][1] - 0.4167293) < 10^(-06)) &
              (abs(res_K6_2$EY_gstar1$CIs[h_iptw_idx][2] - 0.5964627) < 10^(-06)))
  # res_K6$EY_gstar1$other.vars

  # ================================================================
  # COMPARING OLD vs NEW OUTPUT
  # N=1,000
  # NEW INTERFACE FINAL RESULTS WITH MC EVAL (FULL MATCH TO OLD)
  # gIPTW and TMLE_gIPTW AREN'T IMPLEMENTED
  # ================================================================
  # with h method is subsetting (excluding degenerate outcome sA):
                               # old:   # new:
  # epsilon (covariate)      0.02549743 0.02549743
  # alpha (intercept)        0.05410938 0.05410938
  # iptw epsilon (covariate) 0.03556655 0.03556655

                # old:  # new:
  # tmle_A     0.5053725 0.5053725
  # tmle_B     0.5051903 0.5051903
  # iid.tmle_B 0.4475714 0.4475714
  # tmle_iptw  0.5123310 0.5123310
  # iptw_h     0.5065960 0.5065960
  # iptw       0.4910014 0.4910014
  # iid.iptw   0.4429414 0.4429414
  # mle        0.4970377 0.4970377

  #   fWi_init_A   fWi_init_B   fWi_star_A   fWi_star_B
  # -0.008334713 -0.008152518 -0.505372462 -0.505190268

  # NEW VARs:
  # gIPTW and TMLE_gIPTW AREN'T YET IMPLEMENTED
  # $EY_gstar1$vars
  #                      var
  # tmle_A      0.0009265246
  # tmle_B      0.0009268804
  # tmle_g_iptw 0.0069826701
  # h_iptw      0.0021023317
  # g_iptw      0.0000000000
  # mle         0.0000000000

  # NEW CIs
  # $EY_gstar1$CIs
  #             LBCI_0.025 UBCI_0.975
  # tmle_A       0.4457134  0.5650315
  # tmle_B       0.4455197  0.5648608
  # tmle_g_iptw -0.1637792  0.1637792
  # h_iptw       0.4167293  0.5964627
  # g_iptw       0.0000000  0.0000000
  # mle          0.4970377  0.4970377
  # > tmlenet_K6out2$EY_gstar1$other.vars
  #    var_iid.tmle_B var_tmleiptw_2ndO     var_iptw_2ndO var_tmle_A_Q.init var_tmle_B_Q.init
  #      0.0004350965      0.0846013137      0.0000000000      0.0008532711      0.0008535512

  #----------------------------------------------------------------------------------
  # Same as Example 1, but specifying the network with NETIDmat: a matrix of friend row numbers from the input data
  #----------------------------------------------------------------------------------
  Net_str <- df_netKmax6[, "Net_str"]
  IDs_str <- df_netKmax6[, "IDs"]
  net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)
  net_ind_obj$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = ' ')
  NetInd_mat <- net_ind_obj$NetInd
  # save(NetInd_mat_Kmax6, file = "NetInd_mat_Kmax6.rda")

  nF <- net_ind_obj$nF
  print(head(NetInd_mat))
  print(head(nF))
  checkTrue(all.equal(df_netKmax6[,"nFriends"], nF))

  res_K6net <- tmlenet(data = df_netKmax6, Kmax = Kmax, f_gstar1 = f.A_0, sW = def_sW, sA = def_sA,
                      Qform = "Y ~ sum.netW3 + sum.netAW2",
                      hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                      hform.gstar = "netA ~ sum.netW3",
                      Anode = "A", Ynode = "Y",
                      NETIDmat = NetInd_mat,
                      optPars = list(runTMLE = "tmle.intercept",n_MCsims = 1))

  checkTrue(all.equal(res_K6net$EY_gstar1$estimates, res_K6_1a$EY_gstar1$estimates))
  checkTrue(all.equal(res_K6net$EY_gstar1$vars, res_K6_1a$EY_gstar1$vars))
  checkTrue(all.equal(res_K6net$EY_gstar1$CIs, res_K6_1a$EY_gstar1$CIs))
  checkTrue(all.equal(res_K6net$EY_gstar1$other.vars, res_K6_1a$EY_gstar1$other.vars))

  #----------------------------------------------------------------------------------
  # Same as Example 1, but using results of eval.summaries() as input to tmlenet
  #----------------------------------------------------------------------------------
  res <- eval.summaries(sW = def_sW, sA = def_sA, Kmax = Kmax, data = df_netKmax6, NETIDmat = NetInd_mat)
  res_K6_alteval <- tmlenet(DatNet.ObsP0 = res$DatNet.ObsP0, Anode = "A",
                            Qform = "Y ~ sum.netW3 + sum.netAW2",
                            hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                            hform.gstar = "netA ~ sum.netW3",
                            f_gstar1 = f.A_0, optPars = list(n_MCsims = 1))

  checkTrue(all.equal(res_K6_alteval$EY_gstar1$estimates, res_K6_1a$EY_gstar1$estimates))
  checkTrue(all.equal(res_K6_alteval$EY_gstar1$vars, res_K6_1a$EY_gstar1$vars))
  checkTrue(all.equal(res_K6_alteval$EY_gstar1$CIs, res_K6_1a$EY_gstar1$CIs))
  checkTrue(all.equal(res_K6_alteval$EY_gstar1$other.vars, res_K6_1a$EY_gstar1$other.vars))

  #----------------------------------------------------------------------------------
  # Example 3. Same as Example 1 but with with true f_g0
  # *** Note that since f_g0 depends on (W1, netW1, netW2, netW3), these covariates also need to be added to sW summary measure ***
  #----------------------------------------------------------------------------------
  # options(tmlenet.verbose = TRUE)
  def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
              def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE) +
              def.sW(netW1 = W1[[0:Kmax]], netW3 = W3[[1:Kmax]])

  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]]) * W2[[1:Kmax]]), replaceNAw0 = TRUE) +
              def.sA(netA = A[[0:Kmax]])

  # True exposure model (under g0):
  f.A_g0 <- function(data, ...) {
    k <- 6
    rbinom(n = nrow(data), size = 1, prob = f.A(k=k, data = data, ...))
  }

  set.seed(seed=12345)
  res_K6_3 <- tmlenet(data = df_netKmax6, Kmax = Kmax, f_gstar1 = f.A_0, sW = def_sW, sA = def_sA,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Anode = "A", Ynode = "Y",
                    NETIDmat = NetInd_mat,
                    optPars = list(f_g0 = f.A_g0, n_MCsims = 10))

  # Test estimates:
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[tmle_idx] - 0.5049234) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[h_iptw_idx] - 0.4754295) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$estimates[gcomp_idx] - 0.4970377) < 10^(-06))
  # Test asymptotic vars:
  checkTrue(abs(res_K6_3$EY_gstar1$vars[tmle_idx] - 0.0009049428) < 10^(-06))
  checkTrue(abs(res_K6_3$EY_gstar1$vars[h_iptw_idx] - 0.001762281) < 10^(-06))
  # Test CIs:
  checkTrue((abs(res_K6_3$EY_gstar1$CIs[tmle_idx][1] - 0.4459632) < 10^(-06)) &
            (abs(res_K6_3$EY_gstar1$CIs[tmle_idx][2] - 0.5638835) < 10^(-06)))
  checkTrue((abs(res_K6_3$EY_gstar1$CIs[h_iptw_idx][1] - 0.3931511) < 10^(-06)) &
            (abs(res_K6_3$EY_gstar1$CIs[h_iptw_idx][2] - 0.5577079) < 10^(-06)))
  # res_K6$EY_gstar1$other.vars


  #***************************************************************************************
  # EQUIVALENT WAYS TO SPECIFY INTERVENTIONS f_gstar1/f_gstar2.
  # LOWERING THE DIMENSIONALITY OF THE SUMMARY MEASURES.
  #***************************************************************************************
  def_sW <- def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)
  def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE)
  # can define intervention by function f.A_0 that sets everyone's A to constant 0:
  res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                      Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                      NETIDmat = NetInd_mat, optPars = list(n_MCsims = 1))

  checkTrue(abs(res_K6_1$EY_gstar1$estimates[tmle_idx] - 0.4955328) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[h_iptw_idx] - 0.5010294) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[gcomp_idx] - 0.4966022) < 10^(-06))

  # equivalent way to define intervention f.A_0 is to just set f_gstar1 to 0:
  res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                      Anode = "A", Ynode = "Y", f_gstar1 = 0L,
                      NETIDmat = NetInd_mat, optPars = list(n_MCsims = 1))

  checkTrue(abs(res_K6_1$EY_gstar1$estimates[tmle_idx] - 0.4955328) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[h_iptw_idx] - 0.5010294) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[gcomp_idx] - 0.4966022) < 10^(-06))

  # or set f_gstar1 to a vector of 0's of length nrow(data):
  res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                      Anode = "A", Ynode = "Y", f_gstar1 = rep_len(0L, nrow(df_netKmax6)),
                      NETIDmat = NetInd_mat, optPars = list(n_MCsims = 1))

  checkTrue(abs(res_K6_1$EY_gstar1$estimates[tmle_idx] - 0.4955328) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[h_iptw_idx] - 0.5010294) < 10^(-06))
  checkTrue(abs(res_K6_1$EY_gstar1$estimates[gcomp_idx] - 0.4966022) < 10^(-06))


  #***************************************************************************************
  # EXAMPLE WITH SIMULATED DATA FOR 2 FRIENDS AND 1 COVARIATE W1 (SIMULATION 1)
  #***************************************************************************************
  #----------------------------------------------------------------------------------
  # Mean population outcome under stochastic intervention P(A=1)=0.2
  #----------------------------------------------------------------------------------
  data(df_netKmax2)
  head(df_netKmax2)
  Kmax <- 2
  # Define the summary measures:
  def_sW <- def.sW(W1[[0:Kmax]])
  def_sA <- def.sA(A[[0:Kmax]])
  # Define the network matrix:
  net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax2), Kmax = Kmax)
  NetInd_mat <- net_ind_obj$makeNetInd.fromIDs(Net_str = df_netKmax2[, "Net_str"],
                                IDs_str = df_netKmax2[, "IDs"])$NetInd

  options(tmlenet.verbose = FALSE)
  set.seed(seed=123456)
  res_K2_1 <- tmlenet(data = df_netKmax2, Kmax = Kmax, sW = def_sW, sA = def_sA,
                      Anode = "A", Ynode = "Y", f_gstar1 = f.A_.2,
                      NETIDmat = NetInd_mat, optPars = list(n_MCsims = 100))

  res_K2_1$EY_gstar1$estimates
  checkTrue(abs(res_K2_1$EY_gstar1$estimates[tmle_idx] - 0.2933030) < 10^(-06))
  checkTrue(abs(res_K2_1$EY_gstar1$estimates[h_iptw_idx] - 0.2829228) < 10^(-06))
  checkTrue(abs(res_K2_1$EY_gstar1$estimates[gcomp_idx] - 0.2845043) < 10^(-06))

  #----------------------------------------------------------------------------------
  # Average treatment effect (ATE) for two interventions, A=1 vs A=0
  #----------------------------------------------------------------------------------
  res_K2_2 <- tmlenet(data = df_netKmax2, Kmax = Kmax, sW = def_sW, sA = def_sA,
                      Anode = "A", Ynode = "Y", f_gstar1 = 1,
                      NETIDmat = NetInd_mat, 
                      optPars = list(f_gstar2 = 0, n_MCsims = 1))
  # names(res_K2_2)
  # Estimates under f_gstar1:
  # res_K2_2$EY_gstar1$estimates
  # Estimates under f_gstar2:
  # res_K2_2$EY_gstar2$estimates
  # ATE estimates for f_gstar1-f_gstar2:
  # res_K2_2$ATE$estimates
  # (f_gstar1 = 1) - (f_gstar1 = 0):
  checkTrue(abs(res_K2_2$ATE$estimates[tmle_idx] - (0.6611584-0.1997903)) < 10^(-06))
  checkTrue(abs(res_K2_2$ATE$estimates[h_iptw_idx] - (0.6866785-0.1739867)) < 10^(-06))
  checkTrue(abs(res_K2_2$ATE$estimates[gcomp_idx] - (0.5901202-0.2015410)) < 10^(-06))

}
