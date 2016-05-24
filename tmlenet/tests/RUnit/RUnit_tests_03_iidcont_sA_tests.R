# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE sA IN IID DATA
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by  binning, conditional on covariates
# Overall exposure g0 (sA) is a mixture of 3 normals, 
# individual exposure is normal with mu for each observation being a function of (W1,W2,W3), sd = 1;
# ---------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
allNA <- function(x) all(is.na(x))

# ---------------------------------------------------------------------------------
# Test 1. Directly fit a joint density for sA, sA2 ~ W, for sA - continuous (with speedglm and glm.fit)
# ---------------------------------------------------------------------------------
test.simple.fit.density.sA <- function() {
  get.density.sAdat <- function(nsamp = 100000) {
    require(simcausal)
    D <- DAG.empty()
    D <-
    D + node("W1", distr = "rbern", prob = 0.5) +
        node("W2", distr = "rbern", prob = 0.3) +
        node("W3", distr = "rbern", prob = 0.3) +
        node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
        node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
    D <- set.DAG(D)
    datO <- sim(D, n = nsamp, rndseed = 12345)
  }

  get.setW.sAdat <- function(setWvals, nsamp = 100000) {
    require(simcausal)
    W1val <- setWvals[1]; W2val <- setWvals[2]; W3val <- setWvals[3]
    D <- DAG.empty()
    D <-
    D + node("W1", distr = "rconst", const = .(W1val)) +
        node("W2", distr = "rconst", const = .(W2val)) +
        node("W3", distr = "rconst", const = .(W3val)) +
        node("sA", distr = "rnorm", mean = (0.98 * W1 + 0.58 * W2 + 0.33 * W3), sd = 1)
    D <- set.DAG(D)
    datWset <- sim(D, n = nsamp, rndseed = 12345)
    setWmat <- as.matrix(data.frame(W1 = W1val, W2 = W2val, W3 = W3val, sA = seq(-4, 4, by = 0.2)))
    return(list(setWsA = datWset$sA, setWmat = setWmat))
  }

  def.nodeojb <- function(datO) {
    Kmax <- 1
    nodes <- list(Anode = "sA", Wnodes = c("W1", "W2", "W3"), nFnode = "nF")
    def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
    def_sA <- def.sA(sA = "sA")
    netind_cl <- simcausal::NetIndClass$new(nobs = nrow(datO))
    # Define datNetObs:
    datnetW <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sW)
    checkTrue(tmlenet:::is.DatNet(datnetW))
    datnetA <- DatNet$new(netind_cl = netind_cl, nodes = nodes)$make.sVar(Odata = datO, sVar.object = def_sA)
    datNetObs <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
    return(list(datNetObs = datNetObs, netind_cl = netind_cl, def_sA = def_sA, def_sW = def_sW, nodes = nodes))
  }

  nsamp <- 10000
  # nsamp <- 100000
  datO <- get.density.sAdat(nsamp)
  nodeobjs <- def.nodeojb(datO)
  testm.sW <- nodeobjs$def_sW$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)
  testm.sA <- nodeobjs$def_sA$eval.nodeforms(data.df = datO, netind_cl = nodeobjs$netind_cl)

  # Define est_params_list:
  reg.sVars <- list(outvars = c("sA"), predvars = c("W1", "W2", "W3"))
  subset_vars <- lapply(reg.sVars$outvars, function(var) {var})
  sA_class <- nodeobjs$datNetObs$datnetA$type.sVar[reg.sVars$outvars]

  # Put all est_params in RegressionClass (regression with speedglm package)
  print("fitting h_gN based equal.len intervals (default) and speedglm (default): ")
  regclass <- RegressionClass$new(outvar.class = sA_class,
                                  outvar = reg.sVars$outvars,
                                  predvars = reg.sVars$predvars,
                                  subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[3]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[4]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[5]]))

  # (intrvls <- out_ContinSummaryModel$intrvls)
  # (intrvls.width <- diff(intrvls))
  # length(intrvls.width)

  # ord.sVar <- nodeobjs$datNetObs$discretize.sVar(name.sVar = "sA", intervals = out_ContinSummaryModel$intrvls)
  # ord.sVar_bw <- intrvls.width[ord.sVar]
  # print(head(cbind(sA = nodeobjs$datNetObs$dat.sVar[, "sA"], ord.sVar, bw = ord.sVar_bw, nodeobjs$datNetObs$dat.bin.sVar), 5))
  # print("freq count for transformed ord.sVar: "); print(table(ord.sVar))
  # plot(density(ord.sVar))
  # hist(ord.sVar)
  summeas.g0$predict(newdata = nodeobjs$datNetObs)  # summeas.g0$sA_nms
  # Get P(sA|W) for the observed data (W,sA):
  # SHOULD BE SIMILAR TO THE OBSERVED DENSITY OF s.A (but discretized)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  
  print("h_gN fit under speedglm: " %+% mean(h_gN)) # [1] 0.2718823
  checkTrue(abs(mean(h_gN)-0.2718823) < 10^-4)
  # ---------------------------------------------------------------------------------------------------------
  # Plot predicted discretized probs conditional on some values of W's
  # ---------------------------------------------------------------------------------------------------------
  # setWvals <- c(W1 = 0, W2 = 0, W3 = 1)
  # subs <- (datO$W1==setWvals["W1"] & datO$W2==setWvals["W2"] & datO$W3==setWvals["W3"])
  # sum(subs)
  # setWdat_res <- get.setW.sAdat(setWvals, nsamp)
  # plot(density(setWdat_res$setWsA))
  # lines(datO[subs,"sA"], h_gN[subs], type = "p", cex = .3, col = "red")
  # plot(datO[subs,"sA"], h_gN[subs], type = "p", cex = .3, col = "red")
  # lines(density(setWdat_res$setWsA))

  # ---------------------------------------------------------------------------------------------------------
  # Plot all predicted discretized probs together (without conditioning on particular subset of W's)
  # ---------------------------------------------------------------------------------------------------------
  # plot(datO[,"sA"], h_gN, type = "p", cex = .3, col = "red")
  # lines(density(datO[,"sA"]))

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing regressions with stats::glm.fit ****
  # ---------------------------------------------------------------------------------------------------------
  print("fitting h_gN based equal.len intervals (default) and glm.fit: ")
  regclass.gml <- RegressionClass$new(useglm = TRUE, outvar.class = sA_class,
                                      outvar = reg.sVars$outvars,
                                      predvars = reg.sVars$predvars,
                                      subset = subset_vars)
  summeas.g0.glm <- SummariesModel$new(reg = regclass.gml, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0.glm$fit(data = nodeobjs$datNetObs)

  # Test the coef and summary functions for binoutmodel class:
  out_ContinSummaryModel <- summeas.g0.glm$getPsAsW.models()$`P(sA|sW).1`
  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[3]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[4]]))
  # print(tmlenet:::summary.BinOutModel(out_BinModels[[5]]))

  summeas.g0.glm$predict(newdata = nodeobjs$datNetObs)
  h_gN.glm <- summeas.g0.glm$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***

  print("h_gN fit under glm: " %+% mean(h_gN.glm)) # [1] 0.2718823
  checkTrue(abs(mean(h_gN.glm)-0.2718823) < 10^-4)


  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but doing binnning by mass intervals & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  print("fitting h_gN based on bin_bymass = TRUE and speedglm (default): ")
  regclass.binmass <- RegressionClass$new(useglm = FALSE,
                                          bin_bymass = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.binmass, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)
  summeas.g0$predict(newdata = nodeobjs$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  mean(h_gN) # [1] 0.2668144
  checkTrue(abs(mean(h_gN)-0.2668144) < 10^-4)

  # ---------------------------------------------------------------------------------------------------------
  # **** Same fit as before but binning using "dhist" function & regressions with speedglm ****
  # ---------------------------------------------------------------------------------------------------------
  print("fitting h_gN based on dhist intervals and speedglm (default): ")
  regclass.bidhist <- RegressionClass$new(useglm = FALSE,
                                          bin_bymass = FALSE,
                                          bin_bydhist = TRUE,
                                          max_nperbin = 1000,
                                          outvar.class = sA_class,
                                          outvar = reg.sVars$outvars,
                                          predvars = reg.sVars$predvars,
                                          subset = subset_vars)
  summeas.g0 <- SummariesModel$new(reg = regclass.bidhist, DatNet.sWsA.g0 = nodeobjs$datNetObs)
  summeas.g0$fit(data = nodeobjs$datNetObs)

  out_ContinSummaryModel <- summeas.g0$getPsAsW.models()$`P(sA|sW).1`
  print("Intervals for dhist: ")
  print(out_ContinSummaryModel$intrvls)
  # [1] -1003.5240566    -3.5240566    -1.9562682    -0.8766233    -0.2052710     0.2963970     0.7404754     1.2078705
  # [9]     1.6959939     2.3455552     3.3931981     4.9362651  1004.9362651

  out_BinModels <- out_ContinSummaryModel$getPsAsW.models()
  print(tmlenet:::coef.BinOutModel(out_BinModels[[1]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[2]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[3]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[4]]))
  print(tmlenet:::summary.BinOutModel(out_BinModels[[5]]))

  summeas.g0$predict(newdata = nodeobjs$datNetObs)
  h_gN <- summeas.g0$predictAeqa(newdata = nodeobjs$datNetObs) # *** DatNet.sWsA$O.datnetA IS TO BE RENAMED TO $O.O.datnetA for clarity ***
  print("regclass.bidhist mean(h_gN) under speedglm: " %+% mean(h_gN))
  # [1] 0.276875
  # [1] 0.27687500785635
  checkTrue(abs(mean(h_gN)-0.276875) < 10^-4)
}


# ---------------------------------------------------------------------------------------------------------
# Test 2. Running iid TMLE fit for continous sA
# TMLE for causal effect in i.i.d. data with continuous exposure under continuous stochastic intervention;
# intervention g.star is defined by shifting the normal density of observed sA until g.star/g.0 >= 10, 
# then its truncated to be equal to g.0
# ---------------------------------------------------------------------------------------------------------
# Run one TMLE simulation for iid data sampled from get.iid.densityOdat, estimating psi0 under trunced g.star
# ---------------------------------------------------------------------------------------------------------
get.iid.densityOdat <- function(nsamp = 100000, rndseed = NULL, trunc.const = 10, shift.const = 2) {
  require(simcausal)
  # nsamp = 100000
  # rndseed = 12345
  # trunc.const <- 10
  # shift.const <- 2
  D <- DAG.empty()
  D <-
  D + node("W1", distr = "rbern", prob = 0.5) +
      node("W2", distr = "rbern", prob = 0.3) +
      node("W3", distr = "rbern", prob = 0.3) +
      node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
      node("shift",  distr = "rconst", const = .(shift.const)) +
      node("sA", distr = "rnorm", mean = sA.mu, sd = 1) +
      node("r.obs.sA",  distr = "rconst", const = exp(shift * (sA - sA.mu - shift / 2))) +
      node("trunc.c",  distr = "rconst", const = .(trunc.const)) +
      node("untrunc.sA.gstar",  distr = "rconst", const = sA + shift) +
      # node("p.gstar.sA", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - sA.mu)^2)) +
      # node("p.gstar.sA.gstar", distr = "rconst", const = (1/sqrt(2*.(pi))) * exp((-1/2) * (untrunc.sA.gstar - (sA.mu + shift))^2)) +
      node("r.new.sA",  distr = "rconst", const = exp(shift * (untrunc.sA.gstar - sA.mu - shift / 2))) +
      node("trunc.sA.gstar",  distr = "rconst", const = ifelse(r.new.sA > trunc.c, sA, untrunc.sA.gstar)) +
      node("probY", distr = "rconst", const = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y", distr = "rbern", prob = plogis(-0.45 * sA - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("probY.gstar", distr = "rconst", const = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3)) +
      node("Y.gstar", distr = "rbern", prob = plogis(-0.45 * trunc.sA.gstar - 0.5 * W1 - 0.58 * W2 - 0.33 * W3))

  D <- set.DAG(D)
  datO <- sim(D, n = nsamp, rndseed = rndseed)

  head(datO, 100)
  print("mean(datO$Y): " %+% mean(datO$Y));
  # [1] 0.31625
  psi0 <- mean(datO$Y.gstar)
  print("psi0: " %+% psi0)
  # [1] 0.22378

  return(list(psi0 = psi0, datO = datO))
}

run.1sim.tmlenet <- function(nsamp, psi0, Qform, f.gstar, trunc.const = 10, shift.const = 2, n_MCsims = 10) {
  datO <- get.iid.densityOdat(nsamp = nsamp, rndseed = NULL, trunc.const = trunc.const, shift.const = shift.const)$datO
  datO <- datO[,c("W1", "W2", "W3", "sA", "Y")]
  Kmax <- 1
  def_sW <- def.sW(W1 = "W1", W2 = "W2", W3 = "W3")
  def_sA <- def.sA(sA = "sA")
  tmlenet_res <- tmlenet(data = datO, Anode = "sA", Ynode = "Y",
                          Kmax = Kmax,
                          # nFnode = NULL,
                          f_gstar1 = f.gstar,
                          sW = def_sW, sA = def_sA,
                          Qform = Qform,
                          hform.g0 = "sA ~ W1 + W2 + W3",
                          hform.gstar = "sA ~ W1 + W2 + W3",
                          optPars = list(n_MCsims = n_MCsims))
                          # correct Q:
                          # Qform = "Y ~ W1 + W2 + W3 + sA",
                          # misspecified Q:
                          # Qform = "Y ~ W3 + sA",
  CIs <- tmlenet_res$EY_gstar1$CIs
  (tmle_B.CI <- CIs[rownames(CIs)%in%"tmle_B",])
  (h_iptw.CI <- CIs[rownames(CIs)%in%"h_iptw",])
  cover.tmle_B <- ((psi0 <= tmle_B.CI[2]) && (psi0 >= tmle_B.CI[1]))
  cover.h_iptw <- ((psi0 <= h_iptw.CI[2]) && (psi0 >= h_iptw.CI[1]))

  est_mat <- tmlenet_res$EY_gstar1$estimates
  est <- as.vector(est_mat)
  names(est) <- rownames(est_mat)
  cover <- apply(CIs, 1, function(row) ((psi0 <= row[2]) && (psi0 >= row[1])))
  cover2 <- c(tmle_B = cover.tmle_B, h_iptw = cover.h_iptw)
  # print("ests"); print(est)
  # print("CIs"); print(CIs)
  # print("cover"); print(cover)
  # print("cover2"); print(cover2)
  return(list(est = est, cover = cover))
}

# Function that returns a stochastic intervention function intervening on sA, for given shift
create_f.gstar <- function(shift, trunc.const) {
  shift.const <- shift
  trunc.const <- trunc.const
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    sA.mu <- 0.98 * data[,"W1"] + 0.58 * data[,"W2"] + 0.33 * data[,"W3"]
    untrunc.sA <- rnorm(n = nrow(data), mean = sA.mu + shift.const, sd = 1)
    r.new.sA <- exp(shift.const * (untrunc.sA - sA.mu - shift.const / 2))
    trunc.sA <- ifelse(r.new.sA > trunc.const, untrunc.sA - shift.const, untrunc.sA)
    return(trunc.sA)
  }
  return(f.gstar)
}

# NOTE:ADD THIS TO AN EXAMPLE OF STOCHASTIC INTERVENTION:
test.onesim.iid.tmlefit <- function() {
  # ---------------------------------------------------------------------------------------------------------
  trunc.const <- 10
  shift.const <- 1
  nsamp <- 10000
  # nsamp <- 2000
  # ---------------------------------------------------------------------------------------------------------
  # # get true psi.0:
  # datO <- get.iid.densityOdat(nsamp = 100000, rndseed = 12345, trunc.const = trunc.const, shift.const = shift.const)
  # psi0 <- datO$psi0
  # print("psi0: " %+% psi0)
  # [1] "psi0: 0.239584" (shift=1)
  # [1] "psi0: 0.2241085" (shift=2)

  # ---------------------------------------------------------------------------------------------------------
  # Misspecified Q:
  # ---------------------------------------------------------------------------------------------------------
  set.seed(33556)
  Qform.mis <- "Y ~ W3 + sA" # misspecified Q:
  f.gstar <- create_f.gstar(shift = shift.const, trunc.const = trunc.const)
  res <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.mis,
                          f.gstar = f.gstar,  trunc.const = trunc.const, shift.const = shift.const,
                          n_MCsims = 10)
  res$est
  #      tmle    h_iptw     gcomp
  # 0.2360305 0.2441800 0.2131320
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2360305
  # h_IPTW 0.2441800
  # MLE    0.2131320
  # test 1:
  checkTrue(abs(res$est["tmle"]-0.2360305) < 10^-4)
  # test 2:
  checkTrue(abs(res$est["h_iptw"]-0.2441800) < 10^-4)
  # test 3:
  checkTrue(abs(res$est["gcomp"]-0.2131320) < 10^-4)

  # ---------------------------------------------------------------------------------------------------------
  # Correct Q:
  # ---------------------------------------------------------------------------------------------------------
  set.seed(23466)
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA" # correct Q:
  res2 <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.corr,
                          f.gstar = f.gstar,  trunc.const = 10, shift.const = shift.const,
                          n_MCsims = 10)
  res2$est
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2391396
  # h_IPTW 0.2443753
  # MLE    0.2428439
  #      tmle    h_iptw     gcomp
  # 0.2391396 0.2443753 0.2428439
  # test 1:
  checkTrue(abs(res2$est["tmle"]-0.2391396) < 10^-4)
  # test 2:
  checkTrue(abs(res2$est["h_iptw"]-0.2443753) < 10^-4)
  # test 3:
  checkTrue(abs(res2$est["gcomp"]-0.2428439) < 10^-4)



  # ---------------------------------------------------------------------------------------------------------
  # Correct Q w/ glm.fit:
  # ---------------------------------------------------------------------------------------------------------
  old_opts <- tmlenet_options(useglm = FALSE)
  print_tmlenet_opts()
  set.seed(23466)
  Qform.corr <- "Y ~ W1 + W2 + W3 + sA" # correct Q:
  res3 <- run.1sim.tmlenet(nsamp = nsamp, psi0 = 0, Qform = Qform.corr,
                          f.gstar = f.gstar,  trunc.const = 10, shift.const = shift.const,
                          n_MCsims = 10)
  res3$est
  # [1] "new MC.ests mat: "
  #         estimate
  # TMLE   0.2391396
  # h_IPTW 0.2443753
  # MLE    0.2428439
  #      tmle    h_iptw     gcomp
  # 0.2391396 0.2443753 0.2428439
  # test 1:
  checkTrue(abs(res2$est["tmle"]-0.2391396) < 10^-4)
  # test 2:
  checkTrue(abs(res2$est["h_iptw"]-0.2443753) < 10^-4)
  # test 3:
  checkTrue(abs(res2$est["gcomp"]-0.2428439) < 10^-4)





}









