#***************************************************************************************
# Define some summary measures for sW
#***************************************************************************************
def_sW <- def.sW(W1, W2, W3) +
          def.sW(netW1 = W1[[1:Kmax]]) +
          def.sW(netW2 = W2[[1:Kmax]]) +
          def.sW(mean.netW2 = mean(W2[[1:Kmax]]), replaceNAw0 = TRUE) +
          def.sW(sum.netW3 = sum(W3[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Define some summary measures for sA
#***************************************************************************************
def_sA <- def.sA(netA = A[[0:Kmax]]) +
          def.sA(sum.netAW2 = sum((1-A[[1:Kmax]])*W2[[1:Kmax]]), replaceNAw0 = TRUE)

#***************************************************************************************
# Evaluate the summary measures applied to the  (O)bserved data (data.frame) and network
#***************************************************************************************
# load observed data:
data(df_netKmax6)
# load the network ID matrix:
data(NetInd_mat_Kmax6)
res <- eval.summaries(sW = def_sW, sA = def_sA,  Kmax = 6, data = df_netKmax6,
  NETIDmat = NetInd_mat_Kmax6, verbose = TRUE)

#***************************************************************************************
# Contents of the list returned by eval.summaries():
#***************************************************************************************
names(res)
# matrix of sW summary measures:
head(res$sW.matrix)
# matrix of sA summary measures:
head(res$sA.matrix)
# matrix of network IDs:
head(res$NETIDmat)
# Observed data summary measures (sW,sA) and network 
# stored as "DatNet.sWsA" R6 class object:
res$DatNet.ObsP0
class(res$DatNet.ObsP0)

#***************************************************************************************
# Using DatNet.ObsP0 as input to tmlenet():
#***************************************************************************************
options(tmlenet.verbose = FALSE) # set to TRUE to print status messages
res_K6 <- tmlenet(DatNet.ObsP0 = res$DatNet.ObsP0,
                    Qform = "Y ~ sum.netW3 + sum.netAW2",
                    hform.g0 = "netA ~ netW2 + sum.netW3 + nF",
                    hform.gstar = "netA ~ sum.netW3",
                    Anode = "A", Ynode = "Y", f_gstar1 = 0L)

res_K6$EY_gstar1$estimates
res_K6$EY_gstar1$vars
res_K6$EY_gstar1$CIs



