# ---------------------------------------------------------------------------------
# TEST SET 6
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

test.shallowdeep.copy <- function() {
  require(simcausal)
  #------------------------------------------------------------------------------------------------------------
  Kmax <- 5
  trunc.const <- 4
  shift.const <- 1
  Dset <- make_netDAG(Kmax, trunc.const, shift.const)
  #------------------------------------------------------------------------------------------------------------
  f.gstar <- create_f.gstar(shift = shift.const, trunc.const = trunc.const)
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA", net.mean.sA = rowMeans(sA[[1:Kmax]]), replaceNAw0 = TRUE)
  seed <- 12345
  datO <- sim(Dset, n = 100, rndseed = seed)
  mean(datO$Y.gstar) # [1] 0.2136
  checkTrue(abs(mean(datO$Y.gstar) - 0.2) < 10^-4)
  netind_cl <- attributes(datO)$netind_cl
  NetInd_mat <- attributes(datO)$netind_cl$NetInd
  #------------------------------------------------------------------------------------------------------------
  # RUN TMLE (bin by equal.mass):
  #------------------------------------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)
  tmlenet_options(maxNperBin = 100, bin.method="equal.mass")
  print_tmlenet_opts()
  Qform.mis <- "Y ~ W2 + W3 + net.mean.sA" # # misspecified Q:
  datO_input <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  res <- tmlenet(data = datO_input, Kmax = Kmax, sW = def_sW, sA = def_sA, Anode = "sA", Ynode = "Y", f_gstar1 = f.gstar,
                  NETIDmat = NetInd_mat,
                  Qform = Qform.mis,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1))
  #------------------------------------------------------------------------------------------------------------
  # save the fits for h_g0 and h_star:
  #------------------------------------------------------------------------------------------------------------
  h_g0_SummariesModel <- res$EY_gstar1$h_g0_SummariesModel
  h_gstar_SummariesModel <- res$EY_gstar1$h_gstar_SummariesModel
  #-----------------------------------------------------------------------------------------
  # PERFORMING SHALLOW AND A DEEP copy of the objects in SummariesModel:
  # the field private$PsAsW.models is a list of R6 objects -> cloned manually
  #-----------------------------------------------------------------------------------------
  shallow_copy <- h_g0_SummariesModel$clone()
  deep_copy <- h_g0_SummariesModel$clone(deep=TRUE)
  #-----------------------------------------------------------------------------------------
  # 1st layer of saved list of R6 objects:
  #-----------------------------------------------------------------------------------------
  model <- h_g0_SummariesModel$getPsAsW.models()[[1]]
  model$intrvls
  model$intrvls <- NULL
  checkTrue(is.null(h_g0_SummariesModel$getPsAsW.models()[[1]]$intrvls))
  checkTrue(is.null(shallow_copy$getPsAsW.models()[[1]]$intrvls))
  checkTrue(!is.null(deep_copy$getPsAsW.models()[[1]]$intrvls))
  #-----------------------------------------------------------------------------------------
  # RegressionClass pointer in ContinSummaryModel:
  #-----------------------------------------------------------------------------------------
  regchild <- h_g0_SummariesModel$getPsAsW.models()[[1]]$reg
  regchild$intrvls
  regchild$intrvls <- NULL
  checkTrue(is.null(h_g0_SummariesModel$getPsAsW.models()[[1]]$reg$intrvls))
  checkTrue(is.null(shallow_copy$getPsAsW.models()[[1]]$reg$intrvls))
  checkTrue(!is.null(deep_copy$getPsAsW.models()[[1]]$reg$intrvls))
  #-----------------------------------------------------------------------------------------
  # 2nd layer of saved list of R6 objects (BinOutModel):
  #-----------------------------------------------------------------------------------------
  BinOutModel <- h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]
  BinOutModel$bw.j <- NULL
  checkTrue(is.null(h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bw.j))
  checkTrue(is.null(shallow_copy$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bw.j))
  checkTrue(!is.null(deep_copy$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bw.j))
  #-----------------------------------------------------------------------------------------
  # Pointer in BinOutModel to BinDat:
  #-----------------------------------------------------------------------------------------
  bindat <- h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bindat
  bindat$nbins <- NULL
  checkTrue(is.null(h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bindat$nbins))
  checkTrue(is.null(shallow_copy$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bindat$nbins))
  checkTrue(!is.null(deep_copy$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bindat$nbins))
}

test.usefitted.h <- function() {
  require(simcausal)
  #------------------------------------------------------------------------------------------------------------
  Kmax <- 5
  trunc.const <- 4
  shift.const <- 1
  Dset <- make_netDAG(Kmax, trunc.const, shift.const)
  #------------------------------------------------------------------------------------------------------------
  f.gstar <- create_f.gstar(shift = shift.const, trunc.const = trunc.const)
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA", net.mean.sA = rowMeans(sA[[1:Kmax]]), replaceNAw0 = TRUE)

  seed <- 12345
  datO <- sim(Dset, n = 100, rndseed = seed)
  mean(datO$Y.gstar) # [1] 0.2136
  checkTrue(abs(mean(datO$Y.gstar) - 0.2) < 10^-4)
  netind_cl <- attributes(datO)$netind_cl
  NetInd_mat <- attributes(datO)$netind_cl$NetInd

  #------------------------------------------------------------------------------------------------------------
  # RUN TMLE (bin by equal.mass):
  #------------------------------------------------------------------------------------------------------------
  options(tmlenet.verbose = FALSE)
  tmlenet_options(maxNperBin = 100, bin.method="equal.mass")
  print_tmlenet_opts()
  Qform.mis <- "Y ~ W2 + W3 + net.mean.sA" # # misspecified Q:
  datO_input <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  res <- tmlenet(data = datO_input, Kmax = Kmax, sW = def_sW, sA = def_sA, Anode = "sA", Ynode = "Y", f_gstar1 = f.gstar,
                  NETIDmat = NetInd_mat,
                  Qform = Qform.mis,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1))
  print("res$EY_gstar1$estimates: "); print(res$EY_gstar1$estimates)
  # tmle   0.2230316
  # h_iptw 0.2464948
  # gcomp  0.2264758

  #------------------------------------------------------------------------------------------------------------
  # save the fits for h_g0 and h_star:
  #------------------------------------------------------------------------------------------------------------
  h_g0_SummariesModel <- res$EY_gstar1$h_g0_SummariesModel
  h_gstar_SummariesModel <- res$EY_gstar1$h_gstar_SummariesModel

  #-----------------------------------------------------------------------------------------
  # PERFORMING SHALLOW AND A DEEP copy of the objects in SummariesModel:
  # the field private$PsAsW.models is a list of R6 objects -> cloned manually
  #-----------------------------------------------------------------------------------------
  h_g0_shallow_copy <- h_g0_SummariesModel$clone()
  h_g0_deep_copy <- h_g0_SummariesModel$clone(deep=TRUE)

  #-----------------------------------------------------------------------------------------
  # 1st layer of saved list of R6 objects:
  model <- h_g0_SummariesModel$getPsAsW.models()[[1]]
  model$intrvls <- NULL
  # RegressionClass pointer in ContinSummaryModel:
  regchild <- h_g0_SummariesModel$getPsAsW.models()[[1]]$reg
  regchild$intrvls <- NULL
  # 2nd layer of saved list of R6 objects (BinOutModel):
  BinOutModel <- h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]
  BinOutModel$bw.j <- NULL
  # Pointer in BinOutModel to BinDat:
  bindat <- h_g0_SummariesModel$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]$bindat
  bindat$nbins <- NULL

  #------------------------------------------------------------------------------------------------------------
  # Run tmle by with different modeling settings for fitting h_g0 and h_star:
  #------------------------------------------------------------------------------------------------------------
  tmlenet_options(maxNperBin = 100, bin.method="equal.len")
  print_tmlenet_opts()
  Qform.mis <- "Y ~ W2 + W3 + net.mean.sA" # # misspecified Q:
  datO_input <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  
  res2a <- tmlenet(data = datO_input, Kmax = Kmax, sW = def_sW, sA = def_sA, Anode = "sA", Ynode = "Y",
                  f_gstar1 = f.gstar,
                  NETIDmat = NetInd_mat,
                  Qform = Qform.mis,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1))
  print("res2a$EY_gstar1$estimates: "); print(res2a$EY_gstar1$estimates)
  # tmle   0.2264758
  # h_iptw 0.2800000
  # gcomp  0.2264758
  checkTrue((res2a$EY_gstar1$estimates-res$EY_gstar1$estimates)[2,]!=0.0)

  #------------------------------------------------------------------------------------------------------------
  # RUN TMLE WITH EXISTING MODEL FITS FOR h_g0 and h_star:
  #------------------------------------------------------------------------------------------------------------
  # NULL'ed some fields in h_g0_SummariesModel -> will be NULL in the shallow copy => has to fail:
  checkException(
  res2b <- tmlenet(data = datO_input, Kmax = Kmax, sW = def_sW, sA = def_sA, Anode = "sA", Ynode = "Y",
                  f_gstar1 = f.gstar,
                  NETIDmat = NetInd_mat,
                  Qform = Qform.mis,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1,
                    h_g0_SummariesModel = h_g0_shallow_copy,
                    h_gstar_SummariesModel = h_gstar_SummariesModel)))
  # However, the deep copy should be allright and return results equivalent to results from the first tmlenet call:
  res2b <- tmlenet(data = datO_input, Kmax = Kmax, sW = def_sW, sA = def_sA, Anode = "sA", Ynode = "Y",
                  f_gstar1 = f.gstar,
                  NETIDmat = NetInd_mat,
                  Qform = Qform.mis,
                  hform.g0 = "sA + net.mean.sA ~ W1 + W2 + W3",
                  hform.gstar = "sA + net.mean.sA ~ W1 + W2 + W3",
                  optPars = list(n_MCsims = 1,
                    h_g0_SummariesModel = h_g0_deep_copy,
                    h_gstar_SummariesModel = h_gstar_SummariesModel))
  print("res2b$EY_gstar1$estimates: "); print(res2b$EY_gstar1$estimates)
  #         estimate
  # tmle   0.2230316
  # h_iptw 0.2464948
  # gcomp  0.2264758
  checkTrue(all.equal(res2b$EY_gstar1$estimates, res$EY_gstar1$estimates))
  #------------------------------------------------------------------------------------------------------------
  # TEST THAT tmlenet always makes a deep copy of input h model fits:
  #------------------------------------------------------------------------------------------------------------
  new_h_g0_copy <- res2b$EY_gstar1$h_g0_SummariesModel
  model <- h_g0_deep_copy$getPsAsW.models()[[1]]
  model$intrvls <- NULL
  checkTrue(is.null(h_g0_deep_copy$getPsAsW.models()[[1]]$intrvls))
  checkTrue(!is.null(new_h_g0_copy$getPsAsW.models()[[1]]$intrvls))


}
















