#-----------------------------------------------------------------------------
# Fit and Predict the IPTW (clever covariate) for summary measures (sW,sA): 
# P_{g^*}(sA | sW)/P_{g0}(sA | sW)
#-----------------------------------------------------------------------------

# @title Predict h weights under g_0 and g_star using existing m.h.fit model fit
# @name pred.hbars
# @export
# fit models for m_gAi
predict.hbars <- function(newdatnet = NULL, m.h.fit) {
# pred.hbars <- function(newdatnet = NULL, m.h.fit) {
    lbound <- m.h.fit$lbound
    # netA_names <- m.h.fit$m.gAi_vec_g$sA_nms
    # determ_cols_Friend <- m.h.fit$determ_cols_Friend # Should this be saved in m.gAi_vec_g and m.gAi_vec_gstar objects instead?
    # if (!is.null(newdatnet)) {
    #   NetInd_k <- newdatnet$netind_cl$NetInd_k
    #   determ.g_user <- newdatnet$determ.g
    #   determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
    #   determ_cols <- (determ_cols_user | determ_cols_Friend)
    # }
    # use original fitted data for prediction
    if (is.null(newdatnet)) {
      stop("newdatnet argument must be not null; this feature is not implemented")
      # newdatnet <- m.h.fit$cY_mtx_fitted
      # determ_cols <- m.h.fit$determ_cols_fitted
    }
    # PASS ENTIRE newdatnet which will get subset, rather than constructing cY_mtx...
    # if (h_user==FALSE) {
    h_gN <- m.h.fit$summeas.g0$predictAeqa(newdata = newdatnet)
    h_gstar <- m.h.fit$summeas.gstar$predictAeqa(newdata = newdatnet)
    # }
    h_gstar_h_gN <- h_gstar / h_gN
    h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
    h_gstar_h_gN <- bound(h_gstar_h_gN, c(0,1/lbound))
    return(h_gstar_h_gN)
}

# @title Defining and fitting the clever covariate h under g_0 and g_star, i.e. models P(sA[j] | sW,sA[j])
# @name fit.hbars
# @importFrom assertthat assert_that is.count
# @export
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars <- function(DatNet.ObsP0, est_params_list) {
  .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))  # defining the vector of c^A's that needs evaluation under h(c)
  #---------------------------------------------------------------------------------
  # PARAMETERS FOR LOGISTIC ESTIMATION OF h
  #---------------------------------------------------------------------------------
  O.datnetW <- DatNet.ObsP0$datnetW
  O.datnetA <- DatNet.ObsP0$datnetA
  lbound <- est_params_list$lbound
  max_npwt <- est_params_list$max_npwt # NOT IMPLEMENTED
  ng.MCsims <- est_params_list$ng.MCsims  # replace with p adaptive to k: p <- 100*(2^k)

  sW <- est_params_list$sW
  sA <- est_params_list$sA
  h.g0.sVars <- est_params_list$h.g0.sVars
  h.gstar.sVars <- est_params_list$h.gstar.sVars

  f.gstar <- est_params_list$f.gstar
  f.g0 <- est_params_list$f.g0

  h_g0_SummariesModel <- est_params_list$h_g0_SummariesModel
  if (!is.null(h_g0_SummariesModel)) {
    message("NOTE: Predictions for P(sA|sW) under g0 will be based on the model fit in h_g0_SummariesModel," %+%
            "all modeling settings will be ignored")
  }
  h_gstar_SummariesModel <- est_params_list$h_gstar_SummariesModel
  if (!is.null(h_gstar_SummariesModel)) {
    message("NOTE: Predictions for P(sA^*|sW^*) under f_gstar will be based on the model fit in h_g0_SummariesModel," %+%
      " all modeling settings will be ignored")
  }

  h_logit_sep_k <- est_params_list$h_logit_sep_k # NOT IMPLEMENTED
  # h_user=est_params_list$h_user; h_user_fcn=est_params_list$h_user_fcn; NOT IMPLEMENTED

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sW
  #---------------------------------------------------------------------------------
  # Summary measure names / expressions:
  sW.g0_nms <- h.g0.sVars$predvars
  sW.gstar_nms <- h.gstar.sVars$predvars

  # *****
  # Check that these summary measures exist in O.datnetW$names.sVar
  check.sW.g0.exist <- unlist(lapply(sW.g0_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  if (!all(check.sW.g0.exist)) stop("the following predictors from hform.g0 regression could not be located in sW summary measures: " %+%
                                    paste0(sW.g0_nms[!check.sW.g0.exist], collapse = ","))

  check.sW.gstar.exist <- unlist(lapply(sW.gstar_nms, function(sWname) sWname %in% O.datnetW$names.sVar))
  if (!all(check.sW.gstar.exist)) stop("the following predictors from hform.gstar regression could not be located in sW summary measures: " %+%
                                    paste0(sW.gstar_nms[!check.sW.gstar.exist], collapse = ","))

  #---------------------------------------------------------------------------------
  # Getting OBSERVED sA
  #---------------------------------------------------------------------------------
  # Summary measure names / expressions:
  sA_nms_g0 <- h.g0.sVars$outvars
  sA_nms_gstar <- h.gstar.sVars$outvars

  # ***********
  # Check that the outcome summary measures defined by h.g0.sVars$outvars and h.gstar.sVars$outvars are equivalent:
  # NOTE: might comment out in the future and allow different summary measures for sA_nms_g0 and sA_nms_gstar.
  # ***********
  if (!all(sA_nms_g0 == sA_nms_gstar)) stop("the outcome variable names defined by regressions hform.g0 & hform.gstar are not identical;" %+%
                                            " current implementation requires these to be the same.")

  # ***********
  # Check that these summary measures exist in O.datnetW$names.sVar
  check.sAg0.exist <- unlist(lapply(sA_nms_g0, function(sAname) sAname %in% O.datnetA$names.sVar))
  if (!all(check.sAg0.exist)) stop("the following outcomes from hform.g0 regression could not be located in sA summary measures: " %+%
                                    paste0(sA_nms_g0[!check.sAg0.exist], collapse = ","))

  check.sAgstar.exist <- unlist(lapply(sA_nms_gstar, function(sAname) sAname %in% O.datnetA$names.sVar))
  if (!all(check.sAgstar.exist)) stop("the following outcomes from hform.gstar regression could not be located in sA summary measures: " %+%
                                    paste0(sA_nms_gstar[!check.sAgstar.exist], collapse = ","))

  #-----------------------------------------------------------
  # DEFINING SUBSETING EXPRESSIONS (FOR DETERMINISTIC / DEGENERATE sA)
  #-----------------------------------------------------------
  # (1 subset expr per regression P(sA[j]|sA[j-1:0], sW))
  # Old examples of subsetting expressions:
  # based on the variable of gvars$misval (requires passing gvars envir for eval)
  # subset_exprs <- lapply(netvar("determ.g_Friend", c(0:Kmax)), function(var) {var%+%" != "%+%"misval"})
  # based on existing logical determ_g columns (TRUE = degenerate/determ):
  # subset_exprs <- lapply(netvar("determ.g_true", c(0:Kmax)), function(var) {var%+%" != "%+%TRUE})
  #-----------------------------------------------------------
  subsets_expr <- lapply(sA_nms_g0, function(var) {var})  # subsetting by !gvars$misval on sA:

  ##########################################
  # Summary class params:
  ##########################################
  sA_class <- O.datnetA$type.sVar[sA_nms_g0]

  if (gvars$verbose) {
    message("================================================================")
    message("fitting h_g0 with summary measures: ", "P(" %+% paste(sA_nms_g0, collapse = ",") %+% " | " %+% paste(sW.g0_nms, collapse = ",") %+% ")")
    message("================================================================")
  }

  p_h0 <- ifelse(is.null(f.g0), 1, ng.MCsims)
  if (!is.null(f.g0)) {
    if (gvars$verbose) message("generating DatNet.g0 under known g0")
    DatNet.g0 <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
    DatNet.g0$make.dat.sWsA(p = p_h0, f.g_fun = f.g0, sA.object = sA)
    print("head(DatNet.g0$dat.sWsA): "); print(head(DatNet.g0$dat.sWsA))
    # DatNet.g0$make.dat.sWsA(p = p_h0, f.g_fun = f.g0, f.g_args = f.g0_args, sA.object = sA)
  } else {
    DatNet.g0 <- DatNet.ObsP0
  }

  regclass.g0 <- RegressionClass$new(outvar.class = sA_class,
                                        outvar = sA_nms_g0,
                                        predvars = sW.g0_nms,
                                        subset = subsets_expr)

  summeas.g0 <- SummariesModel$new(reg = regclass.g0, DatNet.sWsA.g0 = DatNet.g0)
  if (!is.null(h_g0_SummariesModel)) {
    # 1) verify h_g0_SummariesModel is consistent with summeas.g0
    assert_that(inherits(h_g0_SummariesModel, "SummariesModel"))
    # 2) deep copy model fits in h_g0_SummariesModel to summeas.g0
    summeas.g0 <- h_g0_SummariesModel$clone(deep=TRUE)
  } else {
    summeas.g0$fit(data = DatNet.g0)
  }

  # *********
  # NEED TO PASS obsdat.sW.sA (observed data sWsA) to predict() funs.
  # If !is.null(f.g_fun) then DatNet.g0$dat.sWsA IS NOT THE OBSERVED data (sWsA), but rather sWsA data sampled under known g_0.
  # Option 1: Wipe out DatNet.g0$dat.sWsA with actually observed data - means that we can't use DatNet.g0$dat.sWsA in the future.
  # Option 2: Create a new class DatNet.Obs of DatNet.sWsA - pain in the ass...
  # Going with OPTION 1 for now:
  # Already generated DatNet.ObsP0 in tmlenet:
  h_gN <- summeas.g0$predictAeqa(newdata = DatNet.ObsP0)
  # *********

  if (gvars$verbose) {
    message("================================================================")
    message("fitting h_gstar based on summary measures: ", "P(" %+% paste(sA_nms_gstar, collapse = ",") %+% " | " %+% paste(sW.gstar_nms, collapse = ",") %+% ")")
    message("================================================================")
  }

  DatNet.gstar <- DatNet.sWsA$new(datnetW = O.datnetW, datnetA = O.datnetA)
  DatNet.gstar$make.dat.sWsA(p = ng.MCsims, f.g_fun = f.gstar, sA.object = sA)

  if (gvars$verbose) {
    print("Generated new summary measures by sampling A from f_gstar (DatNet.gstar): "); print(class(DatNet.gstar$dat.sWsA))
    print(dim(DatNet.gstar$dat.sWsA)); print(head(DatNet.gstar$dat.sWsA));
  }

  regclass.gstar <- RegressionClass$new(outvar.class = sA_class,
                                        outvar = sA_nms_gstar,
                                        predvars = sW.gstar_nms,
                                        subset = subsets_expr
                                        )
  # Define Intervals Under g_star to Be The Same as under g0:
  summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.g0)
  # Define Intervals Under g_star Based on Summary Measures Generated under g_star:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.gstar)
  # Define Intervals Under g_star Based on Union of Summary Measures under g_star and g0:
  # summeas.gstar <- SummariesModel$new(reg = regclass.gstar, DatNet.sWsA.g0 = DatNet.g0, datnet.gstar = DatNet.gstar)

  if (!is.null(h_gstar_SummariesModel)) {
    # 1) verify h_gstar_SummariesModel is consistent with summeas.gstar
    assert_that(inherits(h_gstar_SummariesModel, "SummariesModel"))
    # 2) deep copy the object with model fits to summeas.gstar
    summeas.gstar <- h_gstar_SummariesModel$clone(deep=TRUE)
  } else {
    summeas.gstar$fit(data = DatNet.gstar)
  }
  h_gstar <- summeas.gstar$predictAeqa(newdata = DatNet.ObsP0)

  ###########################################
  # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
  ##########################################
  h_gstar_h_gN <- h_gstar / h_gN
  h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
  h_gstar_h_gN <- bound(h_gstar_h_gN, c(0, 1/lbound))

  m.h.fit <- list(summeas.g0 = summeas.g0,
                  summeas.gstar = summeas.gstar,
                  lbound = lbound)

  return(list(h_gstar_h_gN = h_gstar_h_gN, m.h.fit = m.h.fit, DatNet.gstar = DatNet.gstar))
}
