#***************************************************************************************
data(df_netKmax6) # Load observed data
data(NetInd_mat_Kmax6) # Load the network ID matrix
Kmax <- 6 # Max number of friends in the network
#***************************************************************************************

#***************************************************************************************
# EXAMPLES OF INTERVENTION FUNCTIONS:
#***************************************************************************************
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
# Stochastically sets (100*x)% of population to A=1 with probability 0.2:
f.A_.2 <- make_f.gstar(x = 0.2)

#***************************************************************************************
# DEFINING SUMMARY MEASURES:
#***************************************************************************************
def_sW <- def.sW(netW2 = W2[[1:Kmax]]) +
          def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)

def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE) +
          def.sA(netA = A[[0:Kmax]])

#***************************************************************************************
# EVALUATING SUMMARY MEASURES BASED ON INPUT DATA:
#***************************************************************************************
# A helper function that can pre-evaluate the summary measures on (O)bserved data, 
# given one of two ways of specifying the network:
res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
                      NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# Specifying the clever covariate regressions hform.g0 and hform.gstar:
# Left side can consist of any summary names defined by def.sA (as linear terms)
# Right side can consist of any summary names defined by def.sW (as linear terms) & 'nF'
#***************************************************************************************
hform.g01 <- "netA ~ netW2 + sum.netW3 + nF"
hform.gstar1 <- "netA ~ sum.netW3"

# alternatives:
hform.g02 <- "netA + sum.netAW2 ~ netW2 + sum.netW3 + nF"
hform.g03 <- "sum.netAW2 ~ netW2 + sum.netW3"

#***************************************************************************************
# Specifying the outcome regression Qform:
# Left side is ignored (with a warning if not equal to Ynode)
# Right side can by any summary names defined by def.sW, def.sA & 'nF'
#***************************************************************************************
Qform1 <- "Y ~ sum.netW3 + sum.netAW2"
Qform2 <- "Y ~ netA + netW + nF"
Qform3 <- "blah ~ netA + netW + nF"

#***************************************************************************************
# Estimate mean population outcome under deterministic intervention A=0 with 6 friends:
# Estimation with regression formulas.
#***************************************************************************************
# Note that Ynode is optional when Qform is specified;
options(tmlenet.verbose = FALSE) # set to TRUE to print status messages
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", f_gstar1 = 0L,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    IDnode = "IDs", NETIDnode = "Net_str",
                    optPars = list(n_MCsims = 1))

res_K6_1$EY_gstar1$estimates
res_K6_1$EY_gstar1$vars
res_K6_1$EY_gstar1$CIs

#***************************************************************************************
# Run exactly the same estimators as above, 
# but using the output "DatNet.ObsP0" of eval.summaries as input to tmlenet:
#***************************************************************************************
res_K6_1b <- tmlenet(DatNet.ObsP0 = res$DatNet.ObsP0,
                    Anode = "A", f_gstar1 = 0L,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    optPars = list(n_MCsims = 1))

res_K6_1b$EY_gstar1$estimates
res_K6_1b$EY_gstar1$vars
res_K6_1b$EY_gstar1$CIs

#***************************************************************************************
# Same as above but for covariate-based TMLE.
#***************************************************************************************
res_K6_2 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", f_gstar1 = 0L,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    IDnode = "IDs", NETIDnode = "Net_str",
                    optPars = list(runTMLE = "tmle.covariate", n_MCsims = 1))

res_K6_2$EY_gstar1$estimates
res_K6_2$EY_gstar1$vars
res_K6_2$EY_gstar1$CIs

#***************************************************************************************
# SPECIFYING THE NETWORK AS A MATRIX OF FRIEND ROW NUMBERS:
#***************************************************************************************
net_ind_obj <- simcausal::NetIndClass$new(nobs = nrow(df_netKmax6), Kmax = Kmax)
# generating the network matrix from input data:
NetInd_mat <- net_ind_obj$makeNetInd.fromIDs(Net_str = df_netKmax6[, "Net_str"],
                                              IDs_str = df_netKmax6[, "IDs"])$NetInd
nF <- net_ind_obj$nF # number of friends

data(NetInd_mat_Kmax6)
all.equal(NetInd_mat, NetInd_mat_Kmax6) # TRUE

print(head(NetInd_mat))
print(head(nF))
print(all.equal(df_netKmax6[,"nFriends"], nF))

res_K6_net1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                        Anode = "A", f_gstar1 = f.A_0,
                        Qform = "Y ~ sum.netW3 + sum.netAW2",
                        hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                        hform.gstar = "netA ~ sum.netW3",
                        NETIDmat = NetInd_mat,
                        optPars = list(runTMLE = "tmle.intercept", n_MCsims = 1))

all.equal(res_K6_net1$EY_gstar1$estimates, res_K6_1$EY_gstar1$estimates)
all.equal(res_K6_net1$EY_gstar1$vars, res_K6_1$EY_gstar1$vars)
all.equal(res_K6_net1$EY_gstar1$CIs, res_K6_1$EY_gstar1$CIs)

#***************************************************************************************
# EQUIVALENT WAYS OF SPECIFYING INTERVENTIONS f_gstar1/f_gstar2.
# LOWERING THE DIMENSIONALITY OF THE SUMMARY MEASURES.
#***************************************************************************************
def_sW <- def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0=TRUE)
def_sA <- def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0=TRUE)

# can define intervention by function f.A_0 that sets everyone's A to constant 0:
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = f.A_0,
                    NETIDmat = NetInd_mat)
res_K6_1$EY_gstar1$estimates

# equivalent way to define intervention f.A_0 is to just set f_gstar1 to 0:
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = 0L,
                    NETIDmat = NetInd_mat)
res_K6_1$EY_gstar1$estimates

# or set f_gstar1 to a vector of 0's of length nrow(data):
res_K6_1 <- tmlenet(data = df_netKmax6, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = rep_len(0L, nrow(df_netKmax6)),
                    NETIDmat = NetInd_mat)
res_K6_1$EY_gstar1$estimates

#***************************************************************************************
# EXAMPLE WITH SIMULATED DATA FOR AT MOST 2 FRIENDS AND 1 BASELINE COVARIATE W1
#***************************************************************************************
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

#***************************************************************************************
# Mean population outcome under stochastic intervention P(A=1)=0.2
#***************************************************************************************
# Evaluates the target parameter by sampling exposures A from f_gstar1;
# Setting optPars$n_MCsims=100 means that A will be sampled 1000 times (
# for a total sample size of nrow(data)*n_MCsims under f_gstar1)
options(tmlenet.verbose = TRUE)
res_K2_1 <- tmlenet(data = df_netKmax2, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y", f_gstar1 = f.A_.2,
                    NETIDmat = NetInd_mat, optPars = list(n_MCsims = 100))
res_K2_1$EY_gstar1$estimates
res_K2_1$EY_gstar1$vars
res_K2_1$EY_gstar1$CIs

#***************************************************************************************
# Estimating the average treatment effect (ATE) for two static interventions:
# f_gstar1 sets everyone A=1 vs f_gstar2 sets everyone A=0;
#***************************************************************************************
res_K2_2 <- tmlenet(data = df_netKmax2, Kmax = Kmax, sW = def_sW, sA = def_sA,
                    Anode = "A", Ynode = "Y",
                    f_gstar1 = 1,
                    NETIDmat = NetInd_mat,
                    optPars = list(f_gstar2 = 0, n_MCsims = 1))
names(res_K2_2)

# Estimates under f_gstar1:
res_K2_2$EY_gstar1$estimates
res_K2_2$EY_gstar1$vars
res_K2_2$EY_gstar1$CIs

# Estimates under f_gstar2:
res_K2_2$EY_gstar2$estimates
res_K2_2$EY_gstar2$vars
res_K2_2$EY_gstar2$CIs

# ATE estimates for f_gstar1-f_gstar2:
res_K2_2$ATE$estimates
res_K2_2$ATE$vars
res_K2_2$ATE$CIs



