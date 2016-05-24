# ---------------------------------------------------------------------------------
# TEST SET 4
# ---------------------------------------------------------------------------------
# network TMLE fit with continous sA 
# Run one TMLE simulation for network data generated igraph::sample_k_regular (all nodes have the same degree k)
# estimate psi0 under trunced g.star (same as in iid)
# sA is the same as in iid example w sAnet=(sA, rowMeans(sA[[1:Kmax]])).
# ---------------------------------------------------------------------------------------------------------
# TMLE for causal effect in network data with continuous exposure under continuous stochastic intervention;
# intervention g.star is defined by shifting the normal density of observed sA until g.star/g.0 >= 10,
# then its truncated to be equal to g.0
# ---------------------------------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)
run.net.1sim.tmlenet <- function(datO, NetInd_mat, def_sW, def_sA, Kmax, Qform, f.gstar, psi0) {
  datO_input <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  res <- tmlenet(data = datO_input, Anode = "sA", Ynode = "Y",
                  Kmax = Kmax,
                  NETIDmat = NetInd_mat,
                  f_gstar1 = f.gstar,
                  sW = def_sW, sA = def_sA,
                  Qform = Qform,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1))
                  # correct Q:
                  # Qform = "Y ~ W1 + W2 + W3 + sA + net.mean.sA",
                  # misspecified Q:
                  # Qform = "Y ~ W2 + W3 + net.mean.sA",

  CIs <- res$EY_gstar1$CIs
  print("CIs: "); print(CIs)
  (tmle_B.CI <- CIs[rownames(CIs)%in%"tmle_B",])
  (h_iptw.CI <- CIs[rownames(CIs)%in%"h_iptw",])
  cover.tmle_B <- ((psi0 <= tmle_B.CI[2]) && (psi0 >= tmle_B.CI[1]))
  cover.h_iptw <- ((psi0 <= h_iptw.CI[2]) && (psi0 >= h_iptw.CI[1]))

  est_mat <- res$EY_gstar1$estimates
  est <- as.vector(est_mat)
  names(est) <- rownames(est_mat)
  cover <- apply(CIs, 1, function(row) ((psi0 <= row[2]) && (psi0 >= row[1])))
  cover2 <- c(tmle_B = cover.tmle_B, h_iptw = cover.h_iptw)

  return(list(est = est, cover = cover))
}

#------------------------------------------------------------------------------------------------------------
# The user-defined network sampler(s) from igraph (regular graph model)
# Generate regular random graphs with same degree for each node
# Kmax - degree of each node
make_netDAG <- function(Kmax, trunc.const, shift.const) {
  generate.igraph.k.regular <- function(n, Kmax, ...) {
    if (n < 20) Kmax <- 5
    igraph.reg <- igraph::sample_k_regular(no.of.nodes = n, k = Kmax, directed = TRUE, multiple = FALSE)
    # From igraph object to sparse adj. matrix:
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(igraph.reg)
    # From igraph object to simcausal/tmlenet input (NetInd_k, nF, Kmax):
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat)
    print("old Kmax:" %+% Kmax)
    print("new Kmax:" %+% NetInd_out$Kmax)
    print("NetInd_k"); print(head(NetInd_out$NetInd_k))
    if (Kmax < NetInd_out$Kmax) message("new network has larger Kmax value than requested, new Kmax = " %+% NetInd_out$Kmax)
    return(NetInd_out$NetInd_k)
  }
  # graph <- igraph::sample_k_regular(no.of.nodes = 50, k = 10, directed = TRUE, multiple = FALSE)
  # par(mar=c(.1,.1,.1,.1))
  # igraph::plot.igraph(graph,
  #     layout=igraph::layout.fruchterman.reingold,
  #     vertex.size=7,
  #     vertex.label.cex=.5,
  #     edge.arrow.size=.5)

  D <- DAG.empty()
  # Adding the ER model network generator from igraph:
  D <- D + network("NetInd_k", Kmax = Kmax, netfun = "generate.igraph.k.regular")
  D <- D +
      node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("net.mean.sA", distr = "rconst", const = mean(sA[[1:Kmax]])) +
      node("untrunc.sA.gstar",  distr = "rconst", const = sA + shift.const) +
      node("r.new.sA",  distr = "rconst", const =
            exp(shift.const * (untrunc.sA.gstar - sA.mu - shift.const / 2))) +
      node("tr.sA.gstar",  distr = "rconst", const = ifelse(r.new.sA > trunc.const, sA, untrunc.sA.gstar)) +
      node("probY", distr = "rconst", const =
            plogis(-0.35 * sA - 0.20 * mean(sA[[1:Kmax]]) - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y", distr = "rbern", prob = probY) +
      node("probY.gstar", distr = "rconst", const =
            plogis(-0.35 * tr.sA.gstar - 0.20 * mean(tr.sA.gstar[[1:Kmax]]) - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y.gstar", distr = "rbern", prob = probY.gstar)
  Dset <- set.DAG(D)
  return(Dset)
}

# Function that returns a stochastic intervention function intervening on sA, for given shift.const
create_f.gstar <- function(shift, trunc.const) {
  shift.const <- shift
  trunc.const <- trunc.const
  f.gstar <- function(data, ...) {
    sA.mu <- 0.98 * data[,"W1"] + 0.58 * data[,"W2"] + 0.33 * data[,"W3"]
    sA <- data[,"sA"]
    untrunc.sA.gstar <- sA + shift.const
    # ratio of P_g^*(sA=sa|W)/P_g0(sA=sa|W) for sa=sA generated under g^*:
    r.new.sA <- exp(shift.const * (untrunc.sA.gstar - sA.mu - shift.const / 2))
    trunc.sA.gstar <- ifelse(r.new.sA > trunc.const, sA, untrunc.sA.gstar)
    return(trunc.sA.gstar)
  }
  return(f.gstar)
}

test.onesim.net.tmlefit <- function() {
  require(simcausal)
  #------------------------------------------------------------------------------------------------------------
  Kmax <- 5
  trunc.const <- 4
  shift.const <- 1
  Dset <- make_netDAG(Kmax, trunc.const, shift.const)
  #------------------------------------------------------------------------------------------------------------
  # # plotDAG(Dset)
  # rndseed2 <- 54321
  # nsamp <- 50000
  # datFull <- sim(Dset, n = nsamp, rndseed = rndseed2)
  # print(head(datFull, 50))
  # # netind_cl <- attributes(datFull)$netind_cl # to get the network object from sim data:
  # # NetInd_mat <- attributes(datFull)$netind_cl$NetInd # to get the network matrix from sim data:
  # print("mean(datFull$Y): " %+% mean(datFull$Y));
  # # [1] "mean(datFull$Y): 0.301425"
  # psi0 <- mean(datFull$Y.gstar)
  # print("psi0: " %+% psi0)
  # # [1] "psi0: 0.22062" for 50K sample (trunc.const = 4, shift = 1, Kmax = 5)
  # # [1] "psi0: 0.19184" for 50K sample (trunc.const = 7, shift = 1, Kmax = 10)

  f.gstar <- create_f.gstar(shift = shift.const, trunc.const = trunc.const)
  # DEFINE SUMMARY MEASURES:
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA", net.mean.sA = rowMeans(sA[[1:Kmax]]), replaceNAw0 = TRUE)

  seed <- 12345
  datO <- sim(Dset, n = 10000, rndseed = seed)
  mean(datO$Y.gstar) # [1] 0.2136
  checkTrue(abs(mean(datO$Y.gstar) - 0.2136) < 10^-6)

  netind_cl <- attributes(datO)$netind_cl
  NetInd_mat <- attributes(datO)$netind_cl$NetInd
  dim(NetInd_mat)

  getOption("tmlenet.verbose")
  options(tmlenet.verbose = FALSE)
  # tmlenet_options(binByMass = FALSE, useglm = TRUE)
  # tmlenet_options(poolContinVar = TRUE, useglm = FALSE) # to pool by contin outcome:
  print_tmlenet_opts()

  #------------------------------------------------------------------------------------------------------------
  # TEST MISSPECIFIED Q:
  #------------------------------------------------------------------------------------------------------------
  tmlenet_options(maxNperBin = 1000, bin.method="equal.mass")
  print_tmlenet_opts()
  Qform.mis <- "Y ~ W2 + W3 + net.mean.sA" # # misspecified Q:
  timerun <- system.time(
    estres <- run.net.1sim.tmlenet(datO = datO, NetInd_mat = NetInd_mat,
                                    def_sW = def_sW, def_sA = def_sA, Kmax = Kmax,
                                    Qform = Qform.mis, f.gstar = f.gstar, psi0 = 0)
  )
  timerun
  #  user  system elapsed
  # 1.659   0.121   1.791
  estres
  #      tmle    h_iptw     gcomp
  # 0.2159311 0.2181338 0.2648796   # [1] "CIs: "
  #        LBCI_0.025 UBCI_0.975
  # tmle    0.1930713  0.2387909
  # h_iptw  0.1886270  0.2476407
  # gcomp          NA         NA

  # test 1:
  checkTrue(abs(estres$est["tmle"] - 0.2159311) < 10^-6)
  # test 2:
  checkTrue(abs(estres$est["h_iptw"] - 0.2181338) < 10^-6)
  # test 3:
  checkTrue(abs(estres$est["gcomp"] - 0.2648796) < 10^-6)

  #------------------------------------------------------------------------------------------------------------
  # TEST EQUAL LENGTH INTERVAL FITTING METHOD (bin.method = "equal.len") (MISSPECIFIED Q):
  #------------------------------------------------------------------------------------------------------------
  Qform.mis <- "Y ~ W2 + W3 + net.mean.sA" # # misspecified Q:
  tmlenet_options(maxNperBin = 1000, bin.method = "equal.len", nbins = 10)
  # tmlenet_options(binByMass = FALSE, useglm = TRUE)
  # tmlenet_options(poolContinVar = TRUE, useglm = FALSE) # to pool by contin outcome:
  print_tmlenet_opts()
  timerun <- system.time(
    estres_eqlen <- run.net.1sim.tmlenet(datO = datO, NetInd_mat = NetInd_mat,
                                    def_sW = def_sW, def_sA = def_sA, Kmax = Kmax,
                                    Qform = Qform.mis, f.gstar = f.gstar, psi0 = 0)
  )
  timerun
  estres_eqlen
  #      tmle    h_iptw     gcomp
  # 0.2257509 0.2299967 0.2648796
  #        LBCI_0.025 UBCI_0.975
  # tmle    0.2001600  0.2513419
  # h_iptw  0.1961125  0.2638809
  # test 1:
  checkTrue(abs(estres_eqlen$est["tmle"] - 0.2257509) < 10^-6)
  # test 2:
  checkTrue(abs(estres_eqlen$est["h_iptw"] - 0.2299967) < 10^-6)
  # test 3:
  checkTrue(abs(estres_eqlen$est["gcomp"] - 0.2648796) < 10^-6)

  #------------------------------------------------------------------------------------------------------------
  # TEST CORRECT Q:
  #------------------------------------------------------------------------------------------------------------
  tmlenet_options(maxNperBin = 1000, bin.method="equal.mass")
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA + net.mean.sA"
  estres2 <- run.net.1sim.tmlenet(datO = datO, NetInd_mat = NetInd_mat,
                                  def_sW = def_sW, def_sA = def_sA, Kmax = Kmax,
                                  Qform = Qform.corr, f.gstar = f.gstar, psi0 = 0)
  estres2
  #      tmle    h_iptw     gcomp
  # 0.2150487 0.2181338 0.2140926
  # [1] "CIs: "
  #        LBCI_0.025 UBCI_0.975
  # tmle    0.1925783  0.2375190
  # h_iptw  0.1886270  0.2476407
  # gcomp          NA         NA

  # test 1:
  checkTrue(abs(estres2$est["tmle"] - 0.2150487) < 10^-6)
  # test 2:
  checkTrue(abs(estres2$est["h_iptw"] - 0.2181338) < 10^-6)
  # test 3:
  checkTrue(abs(estres2$est["gcomp"] - 0.2140926) < 10^-6)

  #------------------------------------------------------------------------------------------------------------
  # TEST PARALLEL FITTING METHOD FOR LOGISTIC REGRESSIONS (MISSPECIFIED Q):
  # FAILS ON win build (COMMENTING OUT):
  #------------------------------------------------------------------------------------------------------------
  # require(foreach)
  # require(doParallel)
  # registerDoParallel(cores = 2)
  # # registerDoParallel(cores = 1)
  # tmlenet_options(maxNperBin = 1000, bin.method="equal.mass", parfit = TRUE)
  # # tmlenet_options(binByMass = FALSE, useglm = TRUE)
  # # tmlenet_options(poolContinVar = TRUE, useglm = FALSE) # to pool by contin outcome:
  # print_tmlenet_opts()
  # estres_par <- run.net.1sim.tmlenet(datO = datO, NetInd_mat = NetInd_mat,
  #                                   def_sW = def_sW, def_sA = def_sA, Kmax = Kmax,
  #                                   Qform = Qform.mis, f.gstar = f.gstar, psi0 = 0)
  # # test 1:
  # checkTrue(abs(estres_par$est["tmle"] - 0.2159311) < 10^-6)
  # # test 2:
  # checkTrue(abs(estres_par$est["h_iptw"] - 0.2181338) < 10^-6)
  # # test 3:
  # checkTrue(abs(estres_par$est["gcomp"] - 0.2648796) < 10^-6)












}