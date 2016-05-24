# ------------------------------------------------------------------------------------------
# TEST SET 5 FOR CATEGORICAL TREATMENT sA (with network net.sA=sum(sA[[1:Kmax]]/nF))
# ------------------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)
# ------------------------------------------------------------------------------------------
# Simulate network data:
# ------------------------------------------------------------------------------------------
get.net.densityOdat <- function(nsamp = 100000, rndseed = NULL, Kmax = 5, shift = 1, sAmax = 7) {
  require(simcausal)
  options(simcausal.verbose=FALSE)
  `%+%` <- function(a, b) paste0(a, b)
  print("shift: " %+% shift)
  #------------------------------------------------------------------------------------------------------------
  # The user-defined network sampler(s) from igraph (regular graph model)
  # Generate regular random graphs with same degree for each node
  # Kmax - degree of each node
  generate.igraph.k.regular <- function(n, Kmax, ...) {
    if (n < 20) Kmax <- 5
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = n, k = Kmax, directed = TRUE, multiple = FALSE)
    # From igraph object to sparse adj. matrix:
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    # From igraph object to simcausal/tmlenet input (NetInd_k, nF, Kmax):
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    if (Kmax < NetInd_out$Kmax) message("new network has larger Kmax value than requested, new Kmax = " %+% NetInd_out$Kmax)
    message("done generating network")
    return(NetInd_out$NetInd_k)
  }
  #------------------------------------------------------------------------------------------------------------
  # Define DAG:
  D <- DAG.empty()
  # Regular network (everyone gets Kmax friends):
  D <- D + network("NetInd_k", Kmax = Kmax, netfun = "generate.igraph.k.regular")
  D <- D + node("nF", distr = "rconst", const = nF)
  # W1 - iid categorical confounder (5 categories, 0-4):
  nW1cat <- 6
  rbinom2 <- function(n, size, prob) rbinom(n, size = size, prob = prob[1,])
  D <- D + node("W1", distr = "rbinom2", size = (nW1cat-1), prob = c(0.4, 0.5, 0.7, 0.4))
  # W2 - positively correlated with W1:
  prob_W2 <- seq.int(0.45, 0.8, length.out = nW1cat)
  D <- D + node("W2", distr = "rbern",
              prob = (W1==0)*.(prob_W2[1]) + (W1==1)*.(prob_W2[2]) + (W1==2)*.(prob_W2[3]) +
                     (W1==3)*.(prob_W2[4]) + (W1==4)*.(prob_W2[5]) + (W1==5)*.(prob_W2[6])) +
           node("net.W2", distr = "rconst", const = ifelse(nF > 0, sum(W2[[1:Kmax]])/nF, 0)) +
           node("W3", distr = "rbern", prob = 0.6)

  normprob <- function(x) x / sum(x)
  k_arr <-c(1:(sAmax-1))
  pN_0 <- 0.02
  prob_Ni_W1_0 <- normprob(c(pN_0, plogis(-3 - 0 - k_arr / 2)))    # W1=0 probabilities of |sA_i|
  prob_Ni_W1_1 <- normprob(c(pN_0, plogis(-1.5 - 0 - k_arr / 3)))  # W1=1 probabilities of |sA_i|
  prob_Ni_W1_2 <- normprob(c(pN_0, pnorm(-2*abs(2 - k_arr) / 5)))  # W1=2 probabilities of |sA_i|
  prob_Ni_W1_3 <- normprob(c(pN_0, pnorm(-2*abs(3 - k_arr) / 5)))  # W1=3 probabilities of |sA_i|
  prob_Ni_W1_4 <- normprob(c(pN_0, plogis(-4 + 2 * (k_arr - 2))))  # W1=4 probabilities of |sA_i|
  prob_Ni_W1_5 <- normprob(c(pN_0, plogis(-4 + 2 * (k_arr - 3))))  # W1=5 probabilities of |sA_i|

  # Define categorical exposure sA, 1:sAmax, where each sA is influenced by categorical W1
  D <- D + node("sA", distr = "rcategor.int",
                probs = (W1 == 0)*.(prob_Ni_W1_0) + (W1 == 1)*.(prob_Ni_W1_1) +
                        (W1 == 2)*.(prob_Ni_W1_2) + (W1 == 3)*.(prob_Ni_W1_3) +
                        (W1 == 4)*.(prob_Ni_W1_4) + (W1 == 5)*.(prob_Ni_W1_5))
  # Define network part of exposure as mean of sA over nF:
  D <- D + node("net.sA", distr = "rconst", const = ifelse(nF > 0, sum(sA[[1:Kmax]])/nF, 0), replaceNAw0 = TRUE)
  # Define the intervention gstar as sA+shift, if sA+shift <= sAmax, otherwise leave sA as is:
  D <- D +node("sA.gstar",  distr = "rconst", const = ifelse(sA+shift <= sAmax, sA+shift, sA))
  D <- D +node("net.sA.gstar",  distr = "rconst", const = ifelse(nF > 0, sum(sA.gstar[[1:Kmax]])/nF, 0),replaceNAw0 = TRUE)

  D <- D +node("probY", distr = "rconst",
            const = plogis( -0.10*sA +
                            -0.20*net.sA +
                            -0.5*net.W2 +
                            +0.5*W1 - 0.58*W2),
                             # - 0.58*W2 - 0.33*W3),
            replaceNAw0 = TRUE) +
          node("Y", distr = "rbern", prob = probY) +
          node("probY.gstar", distr = "rconst",
            const = plogis( -0.10*sA.gstar +
                            -0.20*net.sA.gstar +
                            -0.5*net.W2 +
                            +0.5*W1 - 0.58*W2),
                             # - 0.58*W2 - 0.33*W3),
            replaceNAw0 = TRUE) +
          node("Y.gstar", distr = "rbern", prob = probY.gstar)

  D <- set.DAG(D)

  # nsamp <- 10000
  # rndseed <- 12345
  datO <- sim(D, n = nsamp, rndseed = rndseed)
  # head(datO)
  # table(datO$W1)
  # #   0    1    2    3    4    5
  # # 467 1706 2895 2721 1694  517
  # table(datO$sA)
  # #   1    2    3    4    5    6    7
  # # 264 1569 1947 1801 1590 1456 1373
  # table(datO$net.sA)
  # # 2 2.2 2.4 2.6 2.8   3 3.2 3.4 3.6 3.8   4 4.2 4.4 4.6 4.8   5 5.2 5.4 5.6 5.8   6 6.2 6.4 6.6 6.8
  # # 7  23  62 100 184 283 387 589 737 833 971 948 932 964 851 666 538 351 244 166  88  39  27   8   2
  # table(datO$sA.gstar)
  # #  2    3    4    5    6    7
  # # 264 1569 1947 1801 1590 2829
  # table(datO$net.sA.gstar)
  # # 3  3.2  3.4  3.6  3.8    4  4.2  4.4  4.6  4.8    5  5.2  5.4  5.6  5.8    6  6.2  6.4  6.6  6.8    7
  # #  7   23   65  101  195  315  446  662  857  989 1077 1085 1044  978  786  573  374  241  121   48   13
  # mean(datO$Y) # [1] 0.3446
  # mean(datO$Y.gstar) # [1] 0.2819

  psi0 <- mean(datO$Y.gstar)
  print("mean(datO$Y): " %+% mean(datO$Y)); print("psi0: " %+% psi0)
  return(list(psi0 = psi0, datO = datO, netind_cl = attributes(datO)$netind_cl, NetInd_mat = attributes(datO)$netind_cl$NetInd))
}

test.catnet.fit.density.iptw <- function() {
  def.nodeojb.net <- function(Kmax, datO, NetInd_mat, gstar = FALSE) {
    if (gstar) {
      Anode <- "sA.gstar"
      def_sA <- def.sA(sA = sA.gstar,
                      net.sA = ifelse(nF > 0, rowSums(sA.gstar[[1:Kmax]])/nF, 0),
                      replaceNAw0 = TRUE)
    } else {
      Anode <- "sA"
      def_sA <- def.sA(sA = sA,
                      net.sA = ifelse(nF > 0, rowSums(sA[[1:Kmax]])/nF, 0),
                      replaceNAw0 = TRUE)
    }
    nodes <- list(Anode = Anode, Wnodes = c("W1", "W2"))
    def_sW <- def.sW(W1 = "W1", W2 = "W2",
                    net.W2 = ifelse(nF > 0, rowSums(W2[[1:Kmax]])/nF, 0),
                    replaceNAw0 = TRUE)
    # directly assign already existing network:
    netind_cl <- simcausal::NetIndClass$new(nobs = nrow(datO), Kmax = Kmax)
    netind_cl$NetInd <- NetInd_mat
    # Define datNetObs:
    datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sW)
    datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sA)
    datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
    return(list(datNetObs = datNetObs, netind_cl = netind_cl, def_sA = def_sA, def_sW = def_sW, nodes = nodes))
  }

  # -------------------------------------------------------------------------------------------
  # simulation parameters:
  nsamp <- 10000
  shift <- 1
  Kmax <- 5
  sAmax <- 7

  # -------------------------------------------------------------------------------------------
  # simulating data:
  tsimdat <- system.time(DAGobj <- get.net.densityOdat(nsamp = nsamp, rndseed = 12345, Kmax = Kmax, shift = shift, sAmax = sAmax))
  print(tsimdat)
  NetInd_mat <- DAGobj$NetInd_mat
  datO <- DAGobj$datO
  # [1] "mean(datO$Y): 0.3446"
  # [1] "psi0: 0.2819"
  nrow(datO)
  length(unique(datO$sA))
  length(unique(datO$net.sA))

  # -------------------------------------------------------------------------------------------
  # defining summary measures def.sA, def.sW and DatNet objects for g0 AND g_star:
  nodeobjs.g0 <- def.nodeojb.net(Kmax = Kmax, datO = datO, NetInd_mat = NetInd_mat)
  nodeobjs.gstar <- def.nodeojb.net(Kmax = Kmax, datO = datO, NetInd_mat = NetInd_mat, gstar = TRUE)
  # head(nodeobjs.g0$datNetObs$mat.sVar); head(nodeobjs.gstar$datNetObs$mat.sVar)
  # sW:
  testm.sW <- nodeobjs.g0$def_sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.g0$netind_cl)
  print(head(datO))
  print("testm.sW"); print(head(testm.sW)); print("testm.sW map"); print(nodeobjs.g0$def_sW$sVar.names.map)
  # sA under g0:
  testm.sA <- nodeobjs.g0$def_sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.g0$netind_cl)
  print("testm.sA"); print(head(testm.sA)); print("testm.sA map"); print(nodeobjs.g0$def_sA$sVar.names.map)
  # sA under gstar:
  testm.sA.gstar <- nodeobjs.gstar$def_sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.gstar$netind_cl)
  print("testm.sW.gstar"); print(head(testm.sA.gstar)); print("testm.sA.gstar map"); print(nodeobjs.gstar$def_sA$sVar.names.map)

  # -------------------------------------------------------------------------------------------
  # Define regression parameters and RegressionClass object that defines regressions sA ~ sW:
  reg.sVars <- list(outvars = c("sA", "net.sA"), predvars = c("W1", "W2", "net.W2"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- nodeobjs.g0$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  regclass.obj <- RegressionClass$new(outvar.class = sA_class,
                                      outvar = reg.sVars$outvars,
                                      predvars = reg.sVars$predvars,
                                      subset = subset_vars,
                                      bin_bymass = TRUE,
                                      max_nperbin = 500,
                                      pool_cont = FALSE,
                                      useglm = FALSE,
                                      parfit = FALSE)
  # -------------------------------------------------------------------------------------------
  # estimating h_g0:
  summeas.g0 <- SummariesModel$new(reg = regclass.obj, DatNet.sWsA.g0 = nodeobjs.g0$datNetObs)
  summeas.g0$fit(data = nodeobjs.g0$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs.g0$datNetObs)
  print("mean(h_gN): " %+% round(mean(h_gN), 5)) # [1] 0.07417"
  # -------------------------------------------------------------------------------------------
  # estimating h_gstar (cont sVar intervals based on observed sA under g0):
  summeas.gstar <- SummariesModel$new(reg = regclass.obj, DatNet.sWsA.g0 = nodeobjs.g0$datNetObs)
  summeas.gstar$fit(data = nodeobjs.gstar$datNetObs)
  h_gstar_obs.sA <- summeas.gstar$predictAeqa(newdata = nodeobjs.g0$datNetObs)
  print("mean(h_gstar_obs.sA): " %+% round(mean(h_gstar_obs.sA), 5)) # 0.0524

  # -------------------------------------------------------------------------------------------
  # IPTW:
  print(max(h_gN)); print(min(h_gN))
  # [1] 0.287992
  # [1] 0.0001123328
  print(max(h_gstar_obs.sA)); print(min(h_gstar_obs.sA));
  # [1] 0.4795622
  # [1] 8.330761e-24

  trim_wt <- 130
  wts <- h_gstar_obs.sA / h_gN
  wts[is.nan(wts)] <- 0
  wts[wts > trim_wt] <- trim_wt
  print(summary(h_gstar_obs.sA/h_gN))
   #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 0.0000  0.1112  0.4062  1.0080  1.1100 26.9200
  print(summary(wts))
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.0000  0.1112  0.4062  1.0080  1.1100 26.9200
  iptw_untrimmed <- mean(datO[,"Y"] * (h_gstar_obs.sA/h_gN))
  iptw_trimmed <- mean(datO[,"Y"] * (wts))

  psi0 <- mean(datO$Y.gstar)
  print("true psi0: " %+% psi0)
  print("iptw (untrimmed): " %+% round(iptw_untrimmed, 6))
  print("iptw (wts trimmed by " %+% trim_wt %+% "): " %+% round(iptw_trimmed, 6))

  # TESTS (For 50K w/ rndseed = 12345):
  checkTrue(abs(psi0 - 0.2819) < 10^-6)
  checkTrue(abs(iptw_untrimmed - 0.292515) < 10^-6)
  checkTrue(abs(iptw_trimmed - 0.292515) < 10^-6)

  # ------------------------------------------------------
  # Benchmark for 50K:
  # ------------------------------------------------------
  # >   summary(h_gstar_obs.sA/h_gN)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #  0.0000  0.1132  0.3799  1.0040  1.0980 42.3200
  # >   summary(wts)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #  0.0000  0.1132  0.3799  1.0040  1.0980 42.3200
  # [1] "true psi0: 0.2901"
  # [1] "iptw (untrimmed): 0.291111"
  # [1] "iptw (wts trimmed by 130): 0.291111"

}


  