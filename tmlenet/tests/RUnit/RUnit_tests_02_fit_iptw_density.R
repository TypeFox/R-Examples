# ------------------------------------------------------------------------------------------
# TEST SET 2
# The IPTW estimator that fits continuous density
# ------------------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)
# ------------------------------------------------------------------------------------------
# takes (directed) igraph network and trim IN edges for each vertex with more than Kmax friends
# removes a random sample of IN edges for highly connected vertices with degree(g, mode = "in") > Kmax
# the number of random edges removed is equal to degree(g, mode = "in")-Kmax for each vertex in g with degree > Kmax
trim.igraphnet <- function(g, Kmax) {
  require(data.table)
  degs <- degree(g, mode = "in") # number of edges going IN, for each vertex:
  high_v <- which(degs > Kmax) # vertex IDs that have too many friends (deg > Kmax)
  # only trim edges when needed:
  if (length(high_v)>0) {
    n_edges_to_del <- degs[high_v] - Kmax # number of edges going IN each high_v that need to be deleted (for each high_v):
    del_edges <- vector(mode="integer")
    del_idx_byvert <- vector(mode="list", length=length(high_v))
    for (i in seq_along(high_v)) {
      del_idx <- igraph::sample_seq(1, degs[high_v][i], n_edges_to_del[i])
      del_idx_byvert[[i]] <- del_idx
    }
    names(del_idx_byvert) <- as.character(high_v)
    edgelist_mat <- as_edgelist(g)
    colnames(edgelist_mat) <- c("from", "to")
    edgelist_mat <- cbind(id=(1:nrow(edgelist_mat)), edgelist_mat)
    edgelist_mat <- edgelist_mat[edgelist_mat[,"to"]%in%high_v,]
    # print(nrow(edgelist_mat)==length(E(g)[to(high_v)])) # should be equal
    edgelist_df <- data.table(edgelist_mat)
    setkeyv(edgelist_df, "to")
    # use the list del_idx_byvert to index the row numbers for edge IDs that should be deleted
    # del_idx_byvert enumerates all edges that should be deleted indexed in 1:n_edges_to_del[v_i]
    delrows <- edgelist_df[, list(delrows=.I[del_idx_byvert[[.GRP]]]), by=to]
    id <- edgelist_df[delrows[["delrows"]], ][["id"]] # obtain IDs of the edges that need to be removed from the network:
    g <- g - E(g)[id] # delete random sample of edges going IN high_v:
  }
  return(g)
}

# ------------------------------------------------------------------------------------------
# Simulate network data
get.net.densityOdat <- function(nsamp = 100000, rndseed = NULL, Kmax = 10, trunc.const = 10, shift = 1) {
  require(simcausal)
  options(simcausal.verbose=FALSE)
  `%+%` <- function(a, b) paste0(a, b)
  print("shift: " %+% shift)
  print("trunc.const: " %+% trunc.const)
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
  # trimmed preferential attachment (BA) model (power law deg distribution):
  generate.igraph.prefattach <- function(n, Kmax, power, zero.appeal, m, ...) {
    g <- sample_pa(n, power = power[1], zero.appeal = zero.appeal[1], m = m[1])
    g <- trim.igraphnet(g, Kmax)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(g) # From igraph object to sparse adj. matrix:
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat) # From igraph object to simcausal/tmlenet input (NetInd_k, nF, Kmax):
    if (Kmax < NetInd_out$Kmax) message("new network has larger Kmax value than requested, new Kmax = " %+% NetInd_out$Kmax)
    message("done generating network")
    return(NetInd_out$NetInd_k)
  }
  #------------------------------------------------------------------------------------------------------------
  # trimmed small world (Watts-Strogatz network) model:
  generate.igraph.smallwld <- function(n, Kmax, dim, nei, p, ...) {
    g <- sample_smallworld(dim = 1, size = n, nei = nei[1], p = p[1], loops = FALSE, multiple = FALSE)
    g <- as.directed(g, mode = c("mutual"))
    g <- trim.igraphnet(g, Kmax)
    sparse_AdjMat <- simcausal::igraph.to.sparseAdjMat(g) # From igraph object to sparse adj. matrix:
    NetInd_out <- simcausal::sparseAdjMat.to.NetInd(sparse_AdjMat) # From igraph object to simcausal/tmlenet input (NetInd_k, nF, Kmax):
    if (Kmax < NetInd_out$Kmax) message("new network has larger Kmax value than requested, new Kmax = " %+% NetInd_out$Kmax)
    message("done generating network")
    return(NetInd_out$NetInd_k)
  }

  #------------------------------------------------------------------------------------------------------------
  # Define DAG:
  D <- DAG.empty()
  # Regular network generator (everyone gets Kmax friends):
  D <- D + network("NetInd_k", Kmax = Kmax, netfun = "generate.igraph.k.regular")
  # Trimmed pref. attachment model network generator (trims all high order nodes to have at most Kmax friends):
  # D <- D + network("NetInd_k", netfun = "generate.igraph.prefattach", Kmax = Kmax, power = 0.5, zero.appeal = 1, m = 10)
  # Trimmed pref. attachment model network generator with low m (m=1):
  # D <- D + network("NetInd_k", netfun = "generate.igraph.prefattach", Kmax = Kmax, power = 0.5, zero.appeal = 1, m = 1)
  # Trimmed small world model network generator:
  # D <- D + network("NetInd_k", netfun = "generate.igraph.smallwld", Kmax = Kmax, dim = 1, nei = 9, p = 0.1)
  # Trimmed small world model network generator:
  # D <- D + network("NetInd_k", netfun = "generate.igraph.smallwld", Kmax = Kmax, dim = 1, nei = 4, p = 0.05)      
  D <- D +
      node("nF", distr = "rconst", const = nF) +
      node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("untrunc.sA.gstar",  distr = "rconst", const = sA + shift) +
      node("r.new.sA",  distr = "rconst", const = exp(shift * (untrunc.sA.gstar - sA.mu - shift / 2))) +
      node("trunc.sA.gstar",  distr = "rconst", const = ifelse(r.new.sA > trunc.const, sA, untrunc.sA.gstar)) +
      node("probY", distr = "rconst",
        const = plogis( -0.35*sA +
                        -0.5*ifelse(nF > 0, sum(sA[[1:Kmax]])/nF, 0) +
                        -0.5*ifelse(nF > 0, sum(W1[[1:Kmax]])/nF, 0) +
                        -0.5*W1 - 0.58*W2 - 0.33*W3),
        replaceNAw0 = TRUE) +
      node("Y", distr = "rbern", prob = probY) +

      node("probY.gstar", distr = "rconst",
        const = plogis( -0.35*trunc.sA.gstar +
                        -0.5*ifelse(nF > 0, sum(trunc.sA.gstar[[1:Kmax]])/nF, 0) +
                        -0.5*ifelse(nF > 0, sum(W1[[1:Kmax]])/nF, 0) +
                        -0.5*W1 - 0.58*W2 - 0.33*W3),
        replaceNAw0 = TRUE) +
      node("Y.gstar", distr = "rbern", prob = probY.gstar)

  D <- set.DAG(D)

  datO <- sim(D, n = nsamp, rndseed = rndseed)
  psi0 <- mean(datO$Y.gstar)
  print("mean(datO$Y): " %+% mean(datO$Y)); print("psi0: " %+% psi0)
  return(list(psi0 = psi0, datO = datO, netind_cl = attributes(datO)$netind_cl, NetInd_mat = attributes(datO)$netind_cl$NetInd))
}

test.net.fit.density.iptw <- function() {
  def.nodeojb.net <- function(Kmax, datO, NetInd_mat, gstar = FALSE) {
    if (gstar) {
      Anode <- "trunc.sA.gstar"
      def_sA <- def.sA(sA = trunc.sA.gstar,
                      net.mean.sA = ifelse(nF > 0, rowSums(trunc.sA.gstar[[1:Kmax]])/nF, 0),
                      replaceNAw0 = TRUE)
    } else {
      Anode <- "sA"
      def_sA <- def.sA(sA = sA,
                      net.mean.sA = ifelse(nF > 0, rowSums(sA[[1:Kmax]])/nF, 0),
                      replaceNAw0 = TRUE)
    }

    nodes <- list(Anode = Anode, Wnodes = c("W1", "W2", "W3"))

    def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3",
                    net.mean.W1 = ifelse(nF > 0, rowSums(W1[[1:Kmax]])/nF, 0),
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

  nsamp <- 10000
  trunc.const <- 10
  shift <- 1
  # Kmax <- 1
  Kmax <- 5
  # Kmax <- 10
  # Kmax_vec <- c(1:20)
  # Kmax <- Kmax_vec[1]
  tsimdat <- system.time(DAGobj <- get.net.densityOdat(nsamp = nsamp, rndseed = 12345, Kmax = Kmax, trunc.const = trunc.const, shift = shift))
  print(tsimdat)
  NetInd_mat <- DAGobj$NetInd_mat
  datO <- DAGobj$datO

  nodeobjs.g0 <- def.nodeojb.net(Kmax = Kmax, datO = datO, NetInd_mat = NetInd_mat)
  nodeobjs.gstar <- def.nodeojb.net(Kmax = Kmax, datO = datO, NetInd_mat = NetInd_mat, gstar = TRUE)

  # g0:
  testm.sW <- nodeobjs.g0$def_sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.g0$netind_cl)
  # names(nodeobjs.g0)
  # head(nodeobjs.g0$datNetObs$mat.sVar)
  print("testm.sW"); print(head(testm.sW)); print("testm.sW map"); print(nodeobjs.g0$def_sW$sVar.names.map); print(head(datO))
  testm.sA <- nodeobjs.g0$def_sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.g0$netind_cl)
  print("testm.sA"); print(head(testm.sA)); print("testm.sA map"); print(nodeobjs.g0$def_sA$sVar.names.map); print(head(datO))

  # gstar:
  testm.sW.gstar <- nodeobjs.gstar$def_sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.gstar$netind_cl)
  testm.sA.gstar <- nodeobjs.gstar$def_sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs.gstar$netind_cl)

  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA", "net.mean.sA"), predvars = c("W1", "W2", "W3", "net.mean.W1"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- nodeobjs.g0$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  # Intervals derived from already created summeas.g0 model
  # intrvl1 <- as.vector(summeas.g0$getPsAsW.models()[[1]]$intrvls[-c(1,23)])
  #manually enter above vector:
  intrvl1 <- c(-3.62548255, -1.12785162, -0.74009659, -0.44554917, -0.21835575, -0.02732957, 0.14833284, 0.31201017, 0.46825334, 0.62713418,
              0.78276456, 0.92997633, 1.08788579, 1.23596099, 1.39740518, 1.56734601, 1.76611869, 1.97482902, 2.26222101,
              2.68823295, 4.61699406)
  # intrvl2 <- as.vector(summeas.g0$getPsAsW.models()[[2]]$intrvls[-c(1,23)])
  intrvl2 <- c(-1.06075521, -0.08038125, 0.10973948, 0.23777585, 0.34791228, 0.43103335, 0.50965150, 0.57614829, 0.64389548, 0.71249435,
              0.77465206, 0.83962149, 0.90341750, 0.97158744, 1.04194207, 1.11968210, 1.20547029, 1.29951751, 1.41824629,
              1.61717215, 2.56542624)
  intervals <- list(intrvl1, intrvl2)
  names(intervals) <- reg.sVars$outvars

  # RegressionClass object that defines the regression for sA ~ sW:
  # regclass.obj <- RegressionClass$new(outvar.class = sA_class,
  #                                                 outvar = reg.sVars$outvars,
  #                                                 predvars = reg.sVars$predvars,
  #                                                 subset = subset_vars,
  #                                                 bin_bymass = FALSE,
  #                                                 nbins = 1000,
  #                                                 # nbins = 500,
  #                                                 # nbins = 100,
  #                                                 # nbins = 50,
  #                                                 # nbins = 20,
  #                                                 useglm = FALSE,
  #                                                 parfit = TRUE
  #                                                 )
  regclass.obj <- RegressionClass$new(outvar.class = sA_class,
                                                outvar = reg.sVars$outvars,
                                                predvars = reg.sVars$predvars,
                                                subset = subset_vars,
                                                bin_bymass = TRUE,
                                                max_nperbin = 500,
                                                pool_cont = FALSE,
                                                useglm = FALSE,
                                                parfit = FALSE,
                                                intrvls = intervals
                                                )

  # regclass.obj <- RegressionClass$new(outvar.class = sA_class,
  #                                               outvar = reg.sVars$outvars,
  #                                               predvars = reg.sVars$predvars,
  #                                               subset = subset_vars,
  #                                               bin_bydhist = TRUE,
  #                                               nbins = 50,
  #                                               # nbins = 70,
  #                                               pool_cont = FALSE,
  #                                               useglm = FALSE,
  #                                               parfit = TRUE
  #                                               )

  # -------------------------------------------------------------------------------------------
  # estimating h_g0
  # -------------------------------------------------------------------------------------------
  summeas.g0 <- SummariesModel$new(reg = regclass.obj, DatNet.sWsA.g0 = nodeobjs.g0$datNetObs)
  summeas.g0$fit(data = nodeobjs.g0$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs.g0$datNetObs)

  # -------------------------------------------------------------------------------------------
  # estimating h_gstar
  # -------------------------------------------------------------------------------------------
  summeas.gstar <- SummariesModel$new(reg = regclass.obj, DatNet.sWsA.g0 = nodeobjs.g0$datNetObs)
  # for intervals based on observed sA under gstar:
  # summeas.gstar <- SummariesModel$new(reg = regclass.obj, DatNet.sWsA.g0 = nodeobjs.gstar$datNetObs)
  summeas.gstar$fit(data = nodeobjs.gstar$datNetObs)
  h_gstar_obs.sA <- summeas.gstar$predictAeqa(newdata = nodeobjs.g0$datNetObs)

  # -------------------------------------------------------------------------------------------
  max(h_gN); min(h_gN)
  max(h_gstar_obs.sA); min(h_gstar_obs.sA);
  trim_wt <- 130
  wts <- h_gstar_obs.sA / h_gN
  wts[is.nan(wts)] <- 0
  wts[wts > trim_wt] <- trim_wt
  summary(h_gstar_obs.sA/h_gN)
  summary(wts)
  iptw_untrimmed <- mean(datO[,"Y"] * (h_gstar_obs.sA/h_gN))
  iptw_trimmed <- mean(datO[,"Y"] * (wts))

  psi0 <- mean(datO$Y.gstar)
  print("true psi0: " %+% psi0)
  print("iptw (untrimmed): " %+% round(iptw_untrimmed, 6))
  print("iptw (wts trimmed by " %+% trim_wt %+% "): " %+% round(iptw_trimmed, 6))

  # test 1:
  checkTrue(abs(psi0 - 0.1138) < 10^-6)
  # test 2:
  checkTrue(abs(iptw_untrimmed - 0.1137114) < 10^-6)
  # test 3:
  checkTrue(abs(iptw_trimmed - 0.1137114) < 10^-6)

  # ------------------------------------------------------
  # Benchmark for 10K:
  # ------------------------------------------------------
  # >   # summeas.gstar$getPsAsW.models()[[2]]$intrvls.width
  # >   max(h_gN); min(h_gN)
  # [1] 0.4679077
  # [1] 3.139113e-05
  # >   max(h_gstar_obs.sA); min(h_gstar_obs.sA);
  # [1] 0.4188395
  # [1] 4.281069e-16
  # > 
  # >   wts <- h_gstar_obs.sA / h_gN
  # >   wts[is.nan(wts)] <- 0
  # >   # wts[wts > 500] <- 500
  # > 
  # >   summary(h_gstar_obs.sA/h_gN)
  #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
  #   0.00000   0.01102   0.06283   1.03100   0.30460 196.60000
  # >   summary(wts)
  #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
  #   0.00000   0.01102   0.06283   1.01900   0.30460 130.00000
  # >   (iptw <- mean(datO[,"Y"] * (wts)))
  # [1] "true psi0: 0.1138"
  # [1] "iptw (untrimmed): 0.113711"
  # [1] "iptw (wts trimmed by 130): 0.113711"
}


  