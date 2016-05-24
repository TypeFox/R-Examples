# @title tmlenet-package
# @docType package
# @author Oleg Sofrygin, Mark J. van der Laan

#' @useDynLib tmlenet
#' @import R6
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices nclass.FD nclass.Sturges nclass.scott
#' @importFrom graphics axis barplot hist par text
#' @importFrom methods is
#' @importFrom stats approx binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm terms var predict glm.control
#' @importFrom utils data head str
#######################################################################
######### BETA VERSION - NOT FOR DISTRIBUTION #########################
#######################################################################

# NETWORK TMLE
# authors: Oleg Sofrygin <sofrygin@berkeley.edu> and Mark van der Laan <laan@berkeley.edu>

#------------------------------------
# TO DO LIST:
#------------------------------------
  # *) Implement data-adaptive weight truncation wrt minimization of MSE (taking max of 5% weight of total as truth)
  # *) Allow SL to fit Q_N, g_N and h (i.e P(A_j|A_{1},..,A_{j-1}, W_i\\inF_j))
#----------------------------------------------------------------------------------
# Input checks:
#----------------------------------------------------------------------------------
# todo 58 (tmlenet, Q.sVars, g.sVars) +0: Check that outvars & predvars in Q.sVars & g.sVars actually exist in sW, sA
# Currently, if sW,sA doesn't exist it will not be included, without any warning/message.
#todo 8 (tmlenet) +0: check all names exist in data (Anode, Ynode, etc...)
#----------------------------------------------------------------------------------
#todo 52 (tmlenet) +0: Accept sA & sW as character vectors / lists passed to tmlenet (in addition to current set-up)
  # When sW / sA are just lists of character vectors need to capture the calling env and call DefineSummariesClass constructor:
    # user.env <- parent.frame()
    # user.env_l <- list(user.env = user.env)
    # sW <- do.call(DefineSummariesClass$new, c(sW, list(type = "sW"), user.env_l))
    # sW.gstar <- do.call(DefineSummariesClass$new, c(sW.gstar, list(type = "sW.gstar"), user.env_l))
    # sA <- do.call(DefineSummariesClass$new, c(sA, list(type = "sA"), user.env_l))  
#todo 53 (tmlenet) +0: If no sVars were defined (default), use netW (Wnode[[0:Kmax]]) and netA for sVars (Anode[[0:Kmax]])
#todo 54 (tmlenet) +0: Check all sVar names are unique
#---------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DatNet.sWsA <- function(DatNet.sWsA) "DatNet.sWsA"%in%class(DatNet.sWsA)
is.DatNet <- function(DatNet) "DatNet"%in%class(DatNet)
#-----------------------------------------------------------------------------
# ALL NETWORK VARIABLE NAMES MUST BE CONSTRUCTED BY CALLING THIS FUNCTION.
# In the future might return the network variable (column vector) itself.
# Helper function that for given variable name (varnm) and friend index (fidx) 
# returns the characeter name of that network variable varnm[fidx], 
# for fidx = 0 (var itself), ..., kmax. fidx can be a vector, in which case a 
# character vector of network names is returned. If varnm is also a vector, a 
# character vector for all possible combinations of (varnm x fidx) is returned.
# OUTPUT format: Varnm_net.j:
#-----------------------------------------------------------------------------
netvar <- function(varnm, fidx) {
  cstr <- function(varnm, fidx) {
    slen <- length(fidx)
    rstr <- vector(mode = "character", length = slen)
    netidxstr <- ! (fidx %in% 0L)
    rstr[netidxstr] <- stringr::str_c('_netF', fidx[netidxstr])  # vs. 1
    # rstr[netidxstr] <- str_c('.net.', fidx[netidxstr])  # vs. 2
    return(stringr::str_c(varnm, rstr))
  }
  if (length(varnm) > 1) {
    return(unlist(lapply(varnm, cstr, fidx)))
  } else {
    return(cstr(varnm, fidx))
  }
}
# Examples:
# netvar("A", (0:5))
# netvar("A", c(0:5, 0, 3))
# netvar(c("A", "W"), c(0:5, 0, 3))
# netvar(c("A", "W"), c(0:5, 0, 3))

#-----------------------------------------------------------------------------
# General utilities / Global Vars
#-----------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)

checkpkgs <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg %+% " package needed for this function to work. Please install it.", call. = FALSE)
    }   
  }
}

# Bound g(A|W) probability within supplied bounds
bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds)
  x[x>max(bounds)] <- max(bounds)
  return(x)
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                            "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

CheckInputs <- function(data, inputparams) {    # Error checking for inputs
  # ....
}

# returns NULL if no factors exist, otherwise return the name of the factor variable(s)
CheckExistFactors <- function(data) {
  testvec <- unlist(lapply(data, is.factor))
  if (any(testvec)) {
    return(names(data)[which(testvec)])
  } else {
    return(NULL)
  }
}

# throw exception if 1) varname doesn't exist; 2) more than one varname is matched
CheckVarNameExists <- function(data, varname) {
  idvar <- names(data) %in% varname
  if (sum(idvar) < 1) stop("variable name " %+% varname %+% " not found in data input")
  if (sum(idvar) > 1) stop("more than one column in the input data has been matched to name " 
                            %+% varname %+% ". Consider renaming some of the columns: " %+% 
                            paste0(names(data)[idvar], collapse=","))
  return(invisible(NULL))
}

# THIS NEEDS MORE EVALUATION, DOESN'T SEEM TO WORK AS INTENDED DURING MC EVALUTION
# (GET HUGE ABSOLUTE VALUE WEIGHTS, THAT HAVE tiny % CONTRIBUTION)
scale_bound <- function(weights, max_npwt, n) {
  # scale weights
  weight_prop <- (weights/sum(weights))   # % contribution to the total weight, scaled by n
  weight_prop_byn <- (weights/sum(weights))*n   # % contribution to the total weight, scaled by n
  # print("weight summary before trunc"); print(summary(weights))
  # print("weight_prop before trunc, %"); print(summary(weight_prop))
  # print("weight_prop before trunc, scaled by n"); print(summary(weight_prop_byn))

  while (any(weight_prop_byn > (max_npwt+5))) {
    weights[which(weight_prop_byn >= max_npwt)] <- (max_npwt / n) * sum(weights)
    weight_prop_byn <- (weights/sum(weights)) * n   # % contribution to the total weight, sclaed by n
    # print("weight summary after trunc"); print(summary(weights))
    # print("weight_prop after trunc, %"); print(summary(weights/sum(weights)))
    # print("weight_prop after trunc, scaled by n"); print(summary(weight_prop_byn))
  }
  return(weights)
}
# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}
# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's data is an DatNet.sWsA object
#---------------------------------------------------------------------------------
get_all_ests <- function(estnames, DatNet.ObsP0, est_params_list) {
  #---------------------------------------------------------------------------------
  # unified estimator naming used throughout the package:
  # c("TMLE", "h_IPTW", "MLE")
  #---------------------------------------------------------------------------------
  # DatNet.ObsP0$det.Y             # TRUE/FALSE for deterministic Y's
  # DatNet.ObsP0$noNA.Ynodevals    # actual observed Y's
  # m.Q.init$getoutvarnm           # reg outvar name (Ynode)
  # DatNet.ObsP0$YnodeVals         # visible Y's with NA for det.Y
  # m.Q.init$getoutvarval          # Yvals used in prediction (with det.Y obs set to NA)
  # m.Q.init$getprobA1             # predictions (for non-DET Y)
  # m.Q.init$getsubset             # valid subset (!det.Y)
  # m.Q.init$reg                   # regression class (Qreg)

  nodes <- DatNet.ObsP0$nodes
  Y <- DatNet.ObsP0$noNA.Ynodevals # actual observed Y's
  determ.Q <- DatNet.ObsP0$det.Y
  m.Q.init <- est_params_list$m.Q.init

  QY.init <- DatNet.ObsP0$noNA.Ynodevals # getting all node vals, inc. deterministic  
  QY.init[!DatNet.ObsP0$det.Y] <- m.Q.init$predict(newdata = DatNet.ObsP0)$getprobA1[!DatNet.ObsP0$det.Y] # getting predictions P(Y=1) for non-DET Y
  off <- qlogis(QY.init)  # offset

  #************************************************
  # h^*/h_N clever covariate:
  #************************************************
  fit.hbars_t <- system.time(fit.hbars.res <- fit.hbars(DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_params_list)) # fit the clever covariate

  DatNet.gstar <- fit.hbars.res$DatNet.gstar
  m.h.fit <- fit.hbars.res$m.h.fit
  h_wts <- fit.hbars.res$h_gstar_h_gN

  #************************************************
  # IPTW_h estimator:
  #************************************************  
  h_IPTW <- Y
  h_IPTW[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]
  h_IPTW <- mean(h_IPTW)

  #************************************************
  # TMLEs
  #************************************************  
  if ("TMLE_A" %in% estnames) {
    #************************************************
    # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
    # #todo 19 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star
    #************************************************
    ctrl <- glm.control(trace = FALSE, maxit = 1000)
    SuppressGivenWarnings(m.Q.star <- glm(Y ~ -1 + h_wts + offset(off), data = data.frame(Y = Y, off = off, h_wts = h_wts),
                                              subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
    QY.star <- Y
    if (!is.na(coef(m.Q.star))) QY.star <- plogis(off + coef(m.Q.star) * h_wts)

  } else if ("TMLE_B" %in% estnames) {
    #************************************************
    # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
    # #todo 20 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star
    #************************************************
    ctrl <- glm.control(trace = FALSE, maxit = 1000)
    SuppressGivenWarnings(m.Q.star <- glm(Y ~ offset(off), data = data.frame(Y = Y, off = off), weights = h_wts,
                                              subset = !determ.Q, family = "quasibinomial", control = ctrl), GetWarningsToSuppress(TRUE))
    QY.star <- Y
    if (!is.na(coef(m.Q.star))) QY.star <- plogis(off + coef(m.Q.star))
  }

  #************************************************
  # (DISABLED) g_IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \\in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	# g_wts <- iptw_est(k = est_params_list$Kmax, data = data, node_l = nodes, m.gN = est_params_list$m.g0N,
  #                      f.gstar = est_params_list$f.gstar, f.g_args = est_params_list$f.g_args, family = "binomial",
  #                      NetInd_k = est_params_list$NetInd_k, lbound = est_params_list$lbound, max_npwt = est_params_list$max_npwt,
  #                      f.g0 = est_params_list$f.g0, f.g0_args = est_params_list$f.g0_args)
  # Y_IPTW_g <- Y
  # Y_IPTW_g[!determ.Q] <- Y[!determ.Q] * g_wts[!determ.Q]
  #************************************************
  # (DISABLED) g_IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  # #todo 21 (get_all_ests) +0: use glm.fit or speedglm.Wfit for m.Q.star_giptw
  #************************************************
	# SuppressGivenWarnings(m.Q.star_giptw <- glm(Y ~ -1 + g_wts + offset(off),
  #                                						data = data.frame(Y = Y, off = off, g_wts = g_wts),
  #                                						subset = !determ.Q, family = "quasibinomial", control = ctrl),
  #                                						GetWarningsToSuppress(TRUE))

  #************************************************
  # Run Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
	#************************************************
  MC_fit_params <- append(est_params_list, list(m.Q.star = m.Q.star))

  syst1 <- system.time(MCS_res <- get.MCS_ests(DatNet.ObsP0 = DatNet.ObsP0, 
                                                DatNet.gstar = DatNet.gstar, 
                                                MC_fit_params = MC_fit_params, 
                                                m.h.fit = m.h.fit))

  ests <- c(TMLE = MCS_res[["TMLE"]],
            h_IPTW = h_IPTW, # IPTW estimator based on h - clever covariate:
            MLE = MCS_res[["MLE"]])

  ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
  ests_mat[, 1] <- ests
  rownames(ests_mat) <- names(ests); colnames(ests_mat) <- "estimate"

  wts_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 1)
  colnames(wts_mat) <- c("h_wts")
  wts_mat[, "h_wts"] <- h_wts

  # Components of asymptotic variance for tmle_net (E_g^*[\bar{Q}_0(A^s,W^s|W^s)]-\psi_0):
  # SUBSTRACT overall estimate of psi_0 from fW_i i-specific components
  fWi_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 1)
  colnames(fWi_mat) <- c("fWi_Qinit")
  fWi_mat[,"fWi_Qinit"] <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]

  QY_mat <- matrix(0L, nrow = DatNet.ObsP0$nobs, ncol = 2)
  colnames(QY_mat) <- c("QY.init", "QY.star")
  QY_mat[,] <- cbind(QY.init, QY.star)


  if (gvars$verbose)  {
    print("time spent fitting new fit.hbars.res:"); print(fit.hbars_t)
    if ("TMLE_A" %in% estnames) {
      parsubmodel_fits <- rbind(coef(m.Q.star))
      rownames(parsubmodel_fits) <- c("epsilon (clever covariate coefficient)")
    } else if ("TMLE_B" %in% estnames) {
      parsubmodel_fits <- rbind(coef(m.Q.star))
      rownames(parsubmodel_fits) <- c("alpha (intercept)")
    }
    print("new parsubmodel_fits: "); print(parsubmodel_fits)
    print("time to run Monte Carlo target param evaluation: "); print(syst1);
    print(
          c(
          fWi_init = mean(fWi_mat[,"fWi_Qinit"] - ests["TMLE"])
          ));
    print("new MC.ests mat: "); print(ests_mat)
  }

  return(list( ests_mat = ests_mat,
               wts_mat = wts_mat,
               fWi_mat = fWi_mat,
               QY_mat = QY_mat,
               h_g0_SummariesModel = m.h.fit$summeas.g0,
               h_gstar_SummariesModel = m.h.fit$summeas.gstar
              ))
}


#---------------------------------------------------------------------------------
# (NEW INTERFACE FOR SPECIFYING regressions for hform.g0, hform.gstar & Qform)
#---------------------------------------------------------------------------------
get_vars_fromlist <- function(varname, sVar.map) {
  if (varname %in% names(sVar.map)) {
    as.vector(sVar.map[[varname]])
  } else {
    varname
  }
}
# Parse the formulas for summary measure names and create a map to actual covariate names in sA & sW
process_regform <- function(regform, sW.map = NULL, sA.map = NULL, NETIDnode = NULL, sep = ' ', NETIDmat = NULL) {
  if (length(regform)==0L) {
    return(list(outvars =  as.vector(unlist(sA.map)), predvars = as.vector(unlist(sW.map))))
  } else {
    # Getting predictors (sW names):
    regformterms <- terms(regform)
    sW.names <- attributes(regformterms)$term.labels 
    sW.names.alt <- colnames(attributes(regformterms)$factors)
    assert_that(all(sW.names == sW.names.alt))

    # Getting outcomes (sA names):
    out.var <- rownames(attributes(regformterms)$factors)[1] # character string
    out.vars.form <- as.formula(". ~ " %+% out.var)
    out.vars.terms <- terms(out.vars.form)
    sA.names <- attributes(out.vars.terms)$term.labels

    outvars <- unlist(lapply(sA.names, get_vars_fromlist, sA.map))
    predvars <- unlist(lapply(sW.names, get_vars_fromlist, sW.map))
    return(list(outvars = outvars, predvars = predvars))
  }
}

#' Evaluate Summary Measures sA and sW
#'
#' Take input data, create a network matrix (when input network matrix not provided) and evaluate the summary measures 
#'  previously defined with functions \code{\link{def.sW}} and \code{\link{def.sA}}. 
#'  This function is called internally by \code{tmlenet} for the evaluation of the summary measures.
#'  The R6 class object named \code{DatNet.ObsP0} that is returned by this function can be supplied as an input to the
#'  \code{tmlenet} function.
#'  When \code{DatNet.ObsP0} is used as an input to \code{tmlenet}, the rest of the input arguments already provided to 
#'  this function can be omitted from the \code{tmlenet} function call.
#' @param data Same as \code{\link{tmlenet}} input argument.
#' @param Kmax Same as \code{\link{tmlenet}} input argument.
#' @param sW Same as \code{\link{tmlenet}} input argument.
#' @param sA Same as \code{\link{tmlenet}} input argument.
#' @param IDnode (Optional) Same as \code{\link{tmlenet}} input argument.
#' @param NETIDnode (Optional) Same as \code{\link{tmlenet}} input argument.
#' @param sep Optional friend ID character separator for friends listed in \code{NETIDnode} column of \code{data}, default is \code{' '}; 
#'  same as \code{\link{tmlenet}} input argument \code{optPars$sep}.
#' @param NETIDmat (Optional) Same as \code{\link{tmlenet}} input argument.
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. 
#'  Turn this on by default using \code{options(tmlenet.verbose=TRUE)}.
#' @return A named list that contains:
#'  \itemize{
#'  \item \code{sW.matrix} - Matrix of evaluated summary measures for \code{sW}.
#'  \item \code{sA.matrix} - Matrix of evaluated summary measures for \code{sA}.
#'  \item \code{NETIDmat} - Network ID matrix that can be used for \code{NETIDmat} input argument to \code{tmlenet}.
#'  \item \code{DatNet.ObsP0} - R6 object of class \code{\link{DatNet.sWsA}} that stores all the summary measures and the network information.
#'    This object be passed to \code{\link{tmlenet}} as an argument, in which case the arguments already provided to \code{eval.summaries} no
#'    longer need to be specified to \code{tmlenet}.
#'  }
#' @seealso \code{\link{tmlenet}} for estimation of network effects and \code{\link{def.sW}} for defining the summary measures.
#' @example tests/examples/3_eval.summaries_examples.R
#' @export
eval.summaries <- function(data, Kmax, sW, sA, IDnode = NULL, NETIDnode = NULL, sep = ' ', NETIDmat = NULL, 
                            verbose = getOption("tmlenet.verbose")) {
  iid_data_flag <- FALSE  # set to true if no network is provided (will run iid TMLE)
  nFnode = "nF"
  #----------------------------------------------------------------------------------
  # SOME INPUT CHECKS
  #----------------------------------------------------------------------------------  
  assert_that(is.data.frame(data))
  assert_that(is.integerish(Kmax))
  Kmax <- as.integer(Kmax)
  assert_that(is.DefineSummariesClass(sW))
  assert_that(is.DefineSummariesClass(sA))
  nobs <- nrow(data)

  # Check no factors exist in the input data:
  check1 <- CheckExistFactors(data)
  if (!is.null(check1)) stop("found factor column(s) in the input data, consider removing or recoding such column(s) as strings: " 
                            %+% paste0(check1, collapse=','))

  if (is.null(NETIDnode) && is.null(NETIDmat)) {
    message("No network (friends) specified by NETIDnode or NETIDmat args, assuming the input data is i.i.d.")
    nFnode <- NULL
    iid_data_flag <- TRUE
    if (missing(Kmax)) Kmax <- 1 # need to allow Kmax = 0
  }

  #----------------------------------------------------------------------------------
  # Create an object with model estimates, data & network information that is passed on to estimation algorithm(s)
  #----------------------------------------------------------------------------------
  netind_cl <- NetIndClass$new(nobs = nobs, Kmax = Kmax)
  if (!is.null(NETIDnode)) {
    assert_that(is.character(NETIDnode))
    Net_str <- as.character(data[, NETIDnode])
    if (!is.null(IDnode)) {
      assert_that(is.character(IDnode))
      IDs_str <- as.character(data[, IDnode])
    } else {
      IDs_str <- NULL
    }
    netind_cl$makeNetInd.fromIDs(Net_str = Net_str, IDs_str = IDs_str, sep = sep)
  } else if (!is.null(NETIDmat)) {
    assert_that(is.matrix(NETIDmat))
    netind_cl$NetInd <- NETIDmat
    netind_cl$make.nF()
  }

  if (verbose) {
    message("evaluated the network ID matrix: "); print(head(netind_cl$NetInd))
    message("evaluated and added to summary measures the number of friends for each observation (nF): "); print(head(netind_cl$nF))
  }

  #----------------------------------------------------------------------------------
  # Test parsing and evaluating the summary measures (in class DefineSummariesClass):
  #----------------------------------------------------------------------------------
  # Testing the evaluation of summary measures:
  sW.matrix <- sW$eval.nodeforms(data.df = data, netind_cl = netind_cl)
  # sW.matrix <- sW$get.mat.sVar(data.df = data, netind_cl = netind_cl, addnFnode = nFnode)
  sA.matrix <- sA$eval.nodeforms(data.df = data, netind_cl = netind_cl)
  # sA.matrix <- sA$get.mat.sVar(data.df = data, netind_cl = netind_cl)
  if (verbose) {
    print("sample matrix of sW summary measurs: : "); print(head(sW.matrix))
    print("sample matrix of sA summary measurs: "); print(head(sA.matrix))
    print("map of sW names to column names: "); print(sW$sVar.names.map)    
    print("map of sA names to column names: "); print(sA$sVar.names.map)
  }

  #---------------------------------------------------------------------------------
  # BUILDING OBSERVED sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
  #---------------------------------------------------------------------------------
  datnetW <- DatNet$new(netind_cl = netind_cl, nFnode = nFnode)
  # datnetW <- DatNet$new(netind_cl = netind_cl, nFnode = nFnode, addnFnode = TRUE)
  datnetW$make.sVar(Odata = data, sVar.object = sW)
  datnetW$fixmiss_sVar() # permanently replace NA values in sW with 0
  datnetA <- DatNet$new(netind_cl = netind_cl)
  datnetA$make.sVar(Odata = data, sVar.object = sA)
  DatNet.ObsP0 <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA)$make.dat.sWsA()
  return(list(sW.matrix = sW.matrix, sA.matrix = sA.matrix, NETIDmat = netind_cl$NetInd, DatNet.ObsP0 = DatNet.ObsP0))
}

#---------------------------------------------------------------------------------
# MAIN TMLENET FUNCTION
#---------------------------------------------------------------------------------
# layers of tmlenet input spec's:
# 1) spec sW, sA, Qform => assumes hform.g0 = "sA ~ sW", hform.gstar = "sA ~ sW"
# 2) spec sW, sA, Qform and hform.g0 => assumes hform.gstar = hform.g0
# 3) 2) + spec separate hform.gstar (outcome have to be the same sA for both hform.g0 & hform.gstar)

#------------------------------------
#' Estimate Average Network Effects For Arbitrary (Stochastic) Interventions
#'
#' Estimate the average network effect among dependent units with known network structure (in presence of
#'  interference and/or spillover) using \strong{TMLE} (targeted maximum likelihood estimation), \strong{IPTW}
#'  (Horvitz-Thompson or the inverse-probability-of-treatment) and \strong{GCOMP} (parametric G-computation formula).
#' @param DatNet.ObsP0 Instance of class \code{\link{DatNet.sWsA}} returned under the same name by the \code{\link{eval.summaries}} function. 
#'  Stores the evaluated summary measures (\code{sW},\code{sA}) for the observed data, along with the network information.
#'  When this argument is specified, the
#'  following arguments no longer need to be provided: \code{data}, \code{Kmax}, \code{sW}, \code{sA},
#'  \code{IDnode}, \code{NETIDnode}, \code{optPars$sep}, \code{NETIDmat}.
#' @param data Observed data, a \code{data.frame} with named columns, containing the baseline covariates, 
#'  exposures (\code{Anode}), the outcomes (\code{Ynode}) and possibly the network column (\code{NETIDnode}), where 
#'  network is specified by a vector of strings of friend IDs, each string using \code{optPars$sep} character to separate different friend IDs 
#'  (default is \code{optPars$sep=' '}).
# @param estimators (NOT IMPLEMENTED) Character vector with estimator names.
#' @param Kmax Integer constant specifying the maximal number of friends for any observation in the input \code{data} data.frame.
#' @param sW Summary measures constructed from baseline covariates alone. This must be an object of class
#'  \code{DefineSummariesClass} that is returned by calling the function \code{\link{def.sW}}.
#' @param sA Summary measures constructed from exposures \code{Anode} and baseline covariates. This must be an object of class
#'  \code{DefineSummariesClass} that is returned by calling the function \code{\link{def.sW}}.
#' @param Anode Exposure (treatment) variable name (column name in \code{data}); exposures can be either binary, categorical or continuous.
#  This variable can be instead specified with argument \code{sA} by adding a call \code{+def.sA(Anode="ExposureVarName")} to \code{sA}.
# @param AnodeDET Optional column name for indicators of deterministic values of exposures in \code{Anode}, 
#  should be coded as (\code{TRUE}/\code{FALSE}) or (\code{1}/\code{0});
#  observations with \code{AnodeDET}=\code{TRUE}/\code{1} are assumed to have deterministically assigned exposures
#' @param Ynode  Outcome variable name (column name in \code{data}), assumed normalized between 0 and 1. This can instead be specified
#'  on the left-side of the regression formula in argument \code{Qform}.
#' @param NETIDmat Network specification via matrix of friend IDs (\code{ncol=Kmax}, \code{nrow=nrow(data)}),
#'  where each row \code{i} is a vector of \code{i}'s friends IDs or \code{i}'s friends row
#'  numbers in \code{data} if \code{IDnode=NULL}. See Details.
#' @param IDnode Subject identifier variable in the input data, if not supplied the network string in
#'  \code{NETIDnode} is assumed to be indexing the row numbers in the input \code{data}
#' @param NETIDnode Network specification by a column name in input \code{data} consisting of strings that identify the unit's friends
#'  by their IDs or their row numbers (two friends are separated by space, e.g., \code{"1 2"}; unit with no friends should have
#'  an empty \code{""} string). See Details.
#' @param f_gstar1 Either a function or a vector of counterfactual exposures. If a function, must return 
#'  a vector of counterfactual exposures evaluated based on the summary measures matrix (\code{sW,sA}) passed as a named 
#'  argument \code{"data"}, therefore, the function in \code{f_gstar1} must have a named argument \code{"data"} in its signature.
#'  The interventions defined by \code{f_gstar1} can be static, dynamic or stochastic. If \code{f_gstar1} is specified as a 
#'  vector, it must be of length \code{nrow(data)} or 1 (constant treatment assigned to all observations).
#'  See Details below and Examples in "EQUIVALENT WAYS OF SPECIFYING INTERVENTION \code{f_gstar1}" for demonstration.
# @param nFnode (Optional) Name of the variable for the number of friends each unit has, this name can then be used
#  inside the summary measures and regression formulas \code{sW}, \code{sA}, \code{Qform}, \code{hform.g0},
#  \code{hform.gstar} and \code{gform}. See Details.
#' @param Qform Regression formula for outcome in \code{Ynode}, when omitted (\code{Ynode}=\code{NULL}), the outcome variable is 
#'  regressed on all variables defined in \code{sW} and \code{sA}. See Details.
#' @param hform.g0 Regression formula for estimating the conditional density of P(\code{sA} | \code{sW}) under \code{g0}
#'  (the observed exposure mechanism), when omitted, this regression is defined by \code{sA~sW} where \code{sA} 
#'  are all summary measures defined by argument \code{sA} and \code{sW} are all baseline summary measures defined by argument \code{sW}.
#' @param hform.gstar Regression formula for estimating the conditional density P(\code{sA} | \code{sW}) under interventions 
#'  \code{f_gstar1} or \code{f_gstar2}.
#'  When omitted, the same regression formula as in \code{hform.g0} will be used for \code{hform.gstar}. See Details.
# @param gform  (Optional) Regression formula for the joint treatment mechanism, g, that includes the product of all
# friends treatments, P(A_i, A_{F_i} | W). See Details.
# @param args_f_g1star (Optional) Additional arguments to be passed to \code{f_gstar1} intervention function
# @param args_f_g2star (Optional) Additional arguments to be passed to \code{f_gstar2} intervention function
#' @param verbose Set to \code{TRUE} to print messages on status and information to the console. 
#'  Turn this on by default using \code{options(tmlenet.verbose=TRUE)}.
#' @param optPars A named list of additional optional parameters to be passed to \code{tmlenet}, such as
#'  \code{alpha}, \code{lbound}, \code{family}, \code{n_MCsims}, \code{runTMLE}, \code{YnodeDET}, \code{f_gstar2}, \code{sep}, 
#'  \code{f_g0}, \code{h_g0_SummariesModel} and
#'  \code{h_gstar_SummariesModel}. See Details below for the description of each parameter.
#((NOT IMPLEMENTED)) @param Q.SL.library SuperLearner libraries for outcome, Q
#((NOT IMPLEMENTED)) @param g.SL.library SuperLearner libraries for treatment mechanism, g
#((NOT IMPLEMENTED)) @param h_f.g0_args Additional arguments to be passed to f_g0
#((NOT IMPLEMENTED)) @param h_user_fcn User supplied function to calculate the clever covariate, h
#((NOT IMPLEMENTED)) @param h_logit_sep_k Flag for fitting a separate logistic regression for each strata of nFnode,
# used during estimation of the clever covariate, h
#' 
#' @section Details:
#' 
# (NOT IMPLEMENTED) When \code{sW} is missing, by default \code{sW} is constructed as follows. 
#  For each \code{"W"} in \code{Wnodes}, add a vector \code{data[,"W"]} as well as all friends covariate values of
#  \code{"W"} to \code{sW} by running \code{netW = def.sW(W[[0:Kmax]], noname = TRUE)}.
#  For each \code{"W"} in \code{Wnode}, a vector of \code{"W"} values of the first friend of observations
#  \code{i} = 1,...,\code{nrow(data)} will be created and named \code{"W_netF1"} 
#  (i.e., variable "W_netF1" is constructed as def.sW(W_netF1 = W[[1]])).
#  Similarly, the vector of "W" values of the jth friend for observations \code{i} = 1, ..., \code{nrow(data)}
#  will be created and named "W_netFj", for j from 1 to \code{Kmax}, with all \code{W_netFj} then being added to
#  \code{sW}.
# 
# (NOT IMPLEMENTED) Similarly, when \code{sA} is missing, it is constructed by running
#  \code{def.sW(netA = A[[0:Kmax]], noname = TRUE)} (assuming \code{"A"} is the value of \code{Anode} in \code{data}), 
#  which combines the column \code{data[,"A"]} with all the friends treatment assignments of variable "A".
# 
#'
#' Note that in case when both arguments \code{NETIDnode} and \code{NETIDmat} are left unspecified the input data are 
#'  assumed independent, i.e., no network dependency between the observations. 
#'  All inference will be performed based on the i.i.d. efficient influence curve for the target
#'  parameter \eqn{(E_{g^*_1}[Y])}. 
#'
#' Also note that the ordering of the friends IDs in \code{NETIDnode} or \code{NETIDmat} is unimportant.
#'
#' A special non-negative-integer-valued variable \code{nF} is automatically calculated each time
#'  \code{tmlenet} or \code{eval.summaries} functions are called. \code{nF} contains the total number of friends for
#'  each observation and it is always added as an additional column to the matrix of the baseline-covariates-based
#'  summary measures \code{sWmat}. The variable \code{nF} can be used in the same ways as any of the column names in
#'  the input data frame \code{data}. In particular, the name \code{nF} can be used inside the summary measure
#'  expressions (calls to functions \code{def.sW} and \code{def.sA}) and inside any of the regression formulas 
#'  (\code{Qform}, \code{hform.g0}, \code{hform.gstar}).
#' 
# However, if the column \code{data[,nFnode]}
# already exists in the input data it will be compared to the automatically calculated values, with an error produced
# if the two variables do not exactly match.
#' The regression formalas in \code{Qform}, \code{hform.g0} and \code{hform.gstar} can include any summary measures names defined in
#'  \code{sW} and \code{sA}, referenced by their individual variable names or by their aggregate summary measure names.
#'  For example, \code{hform.g0 = "netA ~ netW"} is equivalent to
#'  \code{hform.g0 = "A + A_netF1 + A_netF2 ~ W + W_netF1 + W_netF2"} for \code{sW,sA} summary measures defined by
#'  \code{def.sW(netW=W[[0:2]])} and \code{def.sA(netA=A[[0:2]])}.
#'
#' @section Additional parameters:
#' 
#' Some of the parameters that control the estimation in \code{tmlenet} can be set by calling the function 
#'  \code{\link{tmlenet_options}}. 
#'  
#' Additional parameters can be also specified as a named list \code{optPars} argument of the 
#' \code{tmlenet} function. The items that can be specified in \code{optPars} are:
#' \itemize{
#'
#' \item \code{alpha} - alpha-level for CI calculation (0.05 for 95% CIs); 
#'
#' \item \code{lbound} - One value for symmetrical bounds on P(sW | sW).
#'
#' \item \code{family} - Family specification for regression models, defaults to binomial (CURRENTLY ONLY BINOMIAL
#'  FAMILY IS IMPLEMENTED).
#'
#' \item \code{n_MCsims} - Number of Monte-Carlo simulations performed, each of sample size \code{nrow(data)}, 
#'    for generating new exposures under \code{f_gstar1} or \code{f_gstar2} (if specified) or \code{f_g0} (if specified).
#'    These newly generated exposures are utilized when fitting the conditional densities P(\code{sA}|\code{sW})
#'    and when evaluating the substitution estimators \strong{GCOMP} and \strong{TMLE}
#'    under stochastic interventions \code{f_gstar1} or \code{f_gstar2}.
#'
#' \item \code{runTMLE} - Choose which of the two TMLEs to run, "tmle.intercept" or "tmle.covariate". The default is "tmle.intercept".
#'
#' \item \code{YnodeDET} - Optional column name for indicators of deterministic values of outcome \code{Ynode}, 
#'    coded as (\code{TRUE}/\code{FALSE}) or 
#'    (\code{1}/\code{0}); observations with \code{YnodeDET}=\code{TRUE}/\code{1} are assumed to have constant value for their \code{Ynode}.
#'
#' \item \code{f_gstar2} - Either a function or a vector of counterfactual exposure assignments.
#'    Used for estimating contrasts (average treatment effect) for two interventions, if omitted, only the average 
#'    counterfactual outcome under intervention \code{f_gstar1} is estimated. The requirements for \code{f_gstar2}
#'    are identical to those for \code{f_gstar1}.
#'
#' \item \code{sep} - A character separating friend indices for the same observation in \code{NETIDnode}.
#'
#' \item \code{f_g0} - A function for generating true treatment mechanism under observed \code{Anode}, if known (for example in a
#'    randomized trial). This is used for estimating P(\code{sA}|\code{sW}) under \code{g0} by sampling large vector of \code{Anode}
#'    (of length \code{nrow(data)*n_MCsims}) from \code{f_g0} function;
#'
#' \item \code{h_g0_SummariesModel} - Previously fitted model for P(\code{sA}|\code{sW}) under observed exposure mechanism \code{g0}, 
#'    returned by the previous runs of the \code{tmlenet} function. 
#'    This has to be an object of \code{SummariesModel} \pkg{R6} class. When this argument is specified, all predictions 
#'    P(\code{sA}=\code{sa}|\code{sW}=\code{sw}) under \code{g0} will be based on the model fits provided by this argument.
#'
#' \item \code{h_gstar_SummariesModel} - Previously fitted model for P(\code{sA}|\code{sW}) under (stochastic) intervention
#'    specified by \code{f_gstar1} or \code{f_gstar2}. Also an object of \code{SummariesModel} \pkg{R6} class. 
#'    When this argument is specified, the predictions P(\code{sA}=\code{sa}|\code{sW}=\code{sw})
#'    under \code{f_gstar1} or \code{f_gstar2} will be based on the model fits provided by this argument.
#' }
#'
#' @section Specifying the counterfactual intervention function (\code{f_gstar1} and \code{optPars$f_gstar2}):
#' 
#' The functions \code{f_gstar1} and \code{f_gstar2} can only depend on variables specified by the combined matrix
#'  of summary measures (\code{sW},\code{sA}), which is passed using the argument \code{data}. The functions should
#'  return a vector of length \code{nrow(data)} of counterfactual treatments for observations in the input data.
#'
#' @section Specifying the Network of Friends:
#' 
#' The network of friends (connections) for observations in the input \code{data} can be specified in two
#'  alternative ways, using either \code{NETIDnode} or \code{NETIDmat} input arguments.
#' 
#' \code{NETIDnode} - The first (slower) method uses a vector of strings in \code{data[, NETIDnode]}, where each
#'  string \code{i} must contain the space separated IDs or row numbers of all units in \code{data} thought to be
#'  connected to observation i (friends of unit i);
#' 
#' \code{NETIDmat} - An alternative (and faster) method is to pass a matrix with \code{Kmax} columns and nrow(data)
#'  rows, where each row \code{NETIDmat[i,]} is a vector of observation \code{i}'s friends' IDs or \code{i}'s friends'
#'  row numbers in \code{data} if \code{IDnode=NULL}. If observation \code{i} has fewer than \code{Kmax} friends, the
#'  remainder of \code{NETIDmat[i,]} must be filled with \code{NA}s. Note that the ordering of friend indices is
#'  irrelevant.
#' 
#' @section IPTW estimator:
#' **********************************************************************
#' 
#' \itemize{
#' \item As described in the following section, the first step is to construct an estimator \eqn{P_{g_N}(sA | sW)}
#'    for the common (in \code{i}) conditional density \eqn{P_{g_0}(sA | sW)} for common (in \code{i}) unit-level summary
#'    measures (\code{sA},\code{sW}).
#'
#' \item The same fitting algorithm is applied to construct an estimator \eqn{P_{g^*_N}(sA^* | sW^*)} of the common (in \code{i})
#'    conditional density \eqn{P_{g^*}(sA^* | sW^*)} for common (in \code{i}) unit-level summary measures (\code{sA^*},\code{sW^*}) 
#'    implied by the user-supplied stochastic intervention \code{f_gstar1} or \code{f_gstar2} and the observed distribution of \code{W}.
#'
#' \item These two density estimators form the basis of the IPTW estimator,
#'    which is evaluated at the observed N data points \eqn{O_i=(sW_i, sA_i, Y_i), i=1,...,N} and is given by
#'    \deqn{\psi^{IPTW}_n = \sum_{i=1,...,N}{Y_i \frac{P_{g^*_N}(sA^*=sA_i | sW=sW_i)}{P_{g_N}(sA=sA_i | sW=sW_i)}}.}
#' }
#'
#' @section GCOMP estimator:
#' **********************************************************************
#'
#' @section TMLE estimator:
#' **********************************************************************
#' 
#' @section Modeling \code{P(sA|sW)} for summary measures \code{(sA,sW)}:
#' **********************************************************************
#'
#' Non-parametric 
#'  estimation of the common \strong{unit-level} multivariate joint conditional probability model \code{P_g0(sA|sW)}, 
#'  for unit-level summary measures \code{(sA,sW)} generated from the observed exposures and baseline covariates 
#'  \eqn{(A,W)=(A_i,W_i : i=1,...,N)} (their joint density given by \eqn{g_0(A|W)Q(W)}), is performed by first 
#'  constructing the dataset of N summary measures, \eqn{(sA_i,sW_i : i=1,...,N)}, and then fitting the usual i.i.d. MLE 
#'  for the common density \code{P_g0(sA|sW)} based on the pooled N sample of these summary measures.
#'  
#'  Note that \code{sA} can be multivariate and any of its components \code{sA[j]} can be either binary, categorical
#'  or continuous.
#'  The joint probability model for \code{P(sA|sA)} = \code{P(sA[1],...,sA[k]|sA)} can be factorized as
#'  \code{P(sA[1]|sA)} * \code{P(sA[2]|sA, sA[1])} * ... * \code{P(sA[k]|sA, sA[1],...,sA[k-1])},
#'  where each of these conditional probability models is fit separately, depending on the type of the outcome variable
#'  \code{sA[j]}.
#'  
#'  If \code{sA[j]} is binary, the conditional probability \code{P(sA[j]|sW,sA[1],...,sA[j-1])} is evaluated via logistic
#'  regression model.
#'  When \code{sA[j]} is continuous (or categorical), its range will be fist partitioned into \code{K} bins and the
#'  corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), where each bin indicator \code{B_j} is then used as an outcome in a 
#'  separate logistic regression model with predictors given by \code{sW, sA[1],...,sA[k-1]}.
#'  Thus, the joint probability \code{P(sA|sW)} is defined by such a tree of binary logistic regressions.
#' 
#' For simplicity, we now suppose \code{sA} is continuous and univariate and we describe here an algorithm for fitting
#'  \eqn{P_{g_0}(sA | sW)} (the algorithm 
#'  for fitting \eqn{P_{g^*}(sA^* | sW^*)} is equivalent, except that exposure \code{A} is replaced with exposure \code{A^*}
#'  generated under \code{f_gstar1} or \code{f_gstar2} and 
#'  the predictors \code{sW} from the regression formula \code{hform.g0} are replaced with predictors \code{sW^*}
#'  specified by the regression formula \code{hform.gstar}).
#'
#' \enumerate{
#' \item Generate a dataset of N observed continuous summary measures (\code{sa_i}:i=1,...,N) from observed
#'  ((\code{a_i},\code{w_i}):i=1,...,N). Let \code{sa}\\in{\code{sa_i}:i=1,...,M}.
#'
#' \item Divide the range of \code{sA} values into intervals S=(i_1,...,i_M,i_{M+1}) so that any observed data point
#'    \code{sa_i} belongs to one interval in S, namely, 
#'    for each possible value sa of \code{sA} there is k\\in{1,...,M}, such that, i_k < \code{sa} <= i_{k+1}.
#'    Let the mapping B(sa)\\in{1,...,M} denote a unique interval in S for sa, such that, i_{B(sa)} < sa <= i_{B(sa)+1}.
#'    Let bw_{B(sa)}:=i_{B(sa)+1}-i_{B(sa)} be the length of the interval (bandwidth) (i_{B(sa)},i_{B(sa)+1}).
#'    Also define the binary indicators b_1,...,b_M, where b_j:=I(B(sa)=j), for all j <= B(sa) and b_j:=NA for all j>B(sa).
#'    That is we set b_j to missing ones the indicator I(B(sa)=j) jumps from 0 to 1.
#'    Now let \code{sA} denote the random variable for the observed summary measure for one unit
#'    and denote by (B_1,...,B_M) the corresponding random indicators for \code{sA} defined as B_j := I(B(\code{sA}) = j) 
#'    for all j <= B(\code{sA}) and B_j:=NA for all j>B(\code{sA}).
#'
#' \item For each j=1,...,M, fit the logistic regression model for the conditional probability P(B_j = 1 | B_{j-1}=0, sW), i.e., 
#'    at each j this is defined as the conditional probability of B_j jumping from 0 to 1 at bin j, given that B_{j-1}=0 and 
#'    each of these logistic regression models is fit only among the observations that are still at risk of having B_j=1 with B_{j-1}=0.
#'
#' \item Normalize the above conditional probability of B_j jumping from 0 to 1 by its corresponding interval length (bandwidth) bw_j to 
#'    obtain the discrete conditional hazards h_j(sW):=P(B_j = 1 | (B_{j-1}=0, sW) / bw_j, for each j.
#'    For the summary measure \code{sA}, the above conditional hazard h_j(sW) is equal to P(\code{sA} \\in (i_j,i_{j+1}) | \code{sA}>=i_j, sW), 
#'    i.e., this is the probability that \code{sA} falls in the interval (i_j,i_{j+1}), conditional on sW and conditional on the fact that
#'    \code{sA} does not belong to any intervals before j.
#'
#' \item  Finally, for any given data-point \code{(sa,sw)}, evaluate the discretized conditional density for P(\code{sA}=sa|sW=sw) by first 
#'    evaluating the interval number k=B(sa)\\in{1,...,M} for \code{sa} and then computing \\prod{j=1,...,k-1}{1-h_j(sW))*h_k(sW)}
#'    which is equivalent to the joint conditional probability that \code{sa} belongs to the interval (i_k,i_{k+1}) and does not belong
#'    to any of the intervals 1 to k-1, conditional on sW. 
#'  }
#'
#' The evaluation above utilizes a discretization of the fact that any continuous density f of random variable X can be written as f_X(x)=S_X(x)*h_X(x), 
#'  for a continuous density f of X where S_X(x):=P(X>x) is the survival function for X, h_X=P(X>x|X>=x) is the hazard function for X; as well as the fact that
#'  the discretized survival function S_X(x) can be written as a of the hazards for s<x: S_X(x)=\\prod{s<x}h_X(x).
#'
#' @section Three methods for defining bin (interval) cuttoffs for a continuous one-dimenstional summary measure \code{sA[j]}:
#' **********************************************************************
#'
#' There are 3 alternative methods to defining the bin cutoffs S=(i_1,...,i_M,i_{M+1}) for a continuous summary measure
#'  \code{sA}. The choice of which method is used along with other discretization parameters (e.g., total number of
#'  bins) is controlled via the tmlenet_options() function. See \code{?tmlenet_options} argument \code{bin.method} for
#'  additional details.
#'
#' Approach 1 (\code{equal.len}): equal length, default.
#'
#' *********************
#'
#' The bins are defined by splitting the range of observed \code{sA} (sa_1,...,sa_n) into equal length intervals. 
#'  This is the dafault discretization method, set by passing an argument \code{bin.method="equal.len"} to 
#'  \code{tmlenet_options} function prior to calling \code{tmlenet()}. The intervals will be defined by splitting the
#'  range of (sa_1,...,sa_N) into \code{nbins} number of equal length intervals, where \code{nbins} is another argument
#'  of \code{tmlenet_options()} function. When \code{nbins=NA} (the default setting) the actual value of \code{nbins}
#'  is computed at run time by taking the integer value (floor) of \code{n/maxNperBin},
#'  for \code{n} - the total observed sample size and \code{maxNperBin=1000} - another argument of
#'  \code{tmlenet_options()} with the default value 1,000.
#'
#' Approach 2 (\code{equal.mass}): data-adaptive equal mass intervals.
#'
#' *********************
#'
#' The intervals are defined by splitting the range of \code{sA} into non-equal length data-adaptive intervals that 
#'  ensures that each interval contains around 
#'  \code{maxNperBin} observations from (sa_j:j=1,...,N).
#'  This interval definition approach can be selected by passing an argument \code{bin.method="equal.mass"} to 
#'  \code{tmlenet_options()} prior to calling \code{tmlenet()}.
#'  The method ensures that an approximately equal number of observations will belong to each interval, where that number
#'  of observations for each interval
#'  is controlled by setting \code{maxNperBin}. The default setting is \code{maxNperBin=1000} observations per interval.
#'
#' Approach 3 (\code{dhist}): combination of 1 & 2.
#'
#' *********************
#'
#' The data-adaptive approach dhist is a mix of Approaches 1 & 2. See Denby and Mallows "Variations on the Histogram" 
#'  (2009)). This interval definition method is selected by passing an argument \code{bin.method="dhist"} to 
#'  \code{tmlenet_options()}  prior to calling \code{tmlenet()}.
#' 
#' @return A named list with 3 items containing the estimation results for:
#'  \itemize{
#'  \item \code{EY_gstar1} - estimates of the mean counterfactual outcome under (stochastic) intervention function \code{f_gstar1} \eqn{(E_{g^*_1}[Y])}.
#'  \item \code{EY_gstar2} - estimates of the mean counterfactual outcome under (stochastic) intervention function \code{f_gstar2} \eqn{(E_{g^*_2}[Y])}, 
#'    \code{NULL} if \code{f_gstar2} not specified.
#'  \item \code{ATE} - additive treatment effect (\eqn{E_{g^*_1}[Y]} - \eqn{E_{g^*_2}[Y]}) under interventions \code{f_gstar1}
#'    vs. in \code{f_gstar2}, \code{NULL} if \code{f_gstar2} not specified.
#' }
#'  
#' Each list item above is itself a list containing the items:
#'  \itemize{
#'  \item \code{estimates} - various estimates of the target parameter (network population counterfactual mean under 
#'    (stochastic) intervention).
#'  \item \code{vars} - the asymptotic variance estimates, for \strong{IPTW} and \strong{TMLE}.
#'  \item \code{CIs} - CI estimates at \code{alpha} level, for \strong{IPTW} and \strong{TMLE}.
#'  \item \code{other.vars} - Placeholder for future versions.
#'  \item \code{h_g0_SummariesModel} - The model fits for P(\code{sA}|\code{sW}) under observed exposure mechanism 
#'    \code{g0}. This is an object of \code{SummariesModel} \pkg{R6} class.
#'  \item \code{h_gstar_SummariesModel} - The model fits for P(\code{sA}|\code{sW}) under intervention \code{f_gstar1}
#'    or \code{f_gstar2}. This is an object of \code{SummariesModel} \pkg{R6} class.
#' }
#' 
#' Currently implemented estimators are:
#'  \itemize{
#'  \item \code{tmle} - Either weighted regression intercept-based TMLE (\code{tmle.intercept} - the default)
#'    with weights defined by the IPTW weights \code{h_gstar/h_gN} or 
#'    covariate-based unweighted TMLE (\code{tmle.covariate}) that uses the IPTW weights as a covariate 
#'    \code{h_gstar/h_gN}.
#'  \item \code{h_iptw} - Efficient IPTW based on weights h_gstar/h_gN.
#'  \item \code{gcomp} - Parametric G-computation formula substitution estimator.
#' }
#' @seealso \code{\link{tmlenet-package}} for the general overview of the package,
#'  \code{\link{def.sW}} for defining the summary measures, \code{\link{eval.summaries}} for
#'  evaluation and validation of the summary measures,
#'  and \code{\link{df_netKmax2}}/\code{\link{df_netKmax6}}/\code{\link{NetInd_mat_Kmax6}}
#'  for examples of network datasets.
#' @example tests/examples/1_tmlenet_example.R
#' @export
tmlenet <- function(DatNet.ObsP0, data, Kmax, sW, sA, Anode, Ynode, f_gstar1,
                    Qform = NULL, hform.g0 = NULL, hform.gstar = NULL,
                    # estimators = c("tmle", "iptw", "gcomp"),
                    # AnodeDET = NULL, 
                    # YnodeDET = NULL,
                    NETIDmat = NULL,
                    IDnode = NULL, NETIDnode = NULL, 
                    verbose = getOption("tmlenet.verbose"),
                    optPars = list(
                      alpha = 0.05,
                      lbound = 0.005,
                      family = "binomial", # NOT YET IMPLEMENTED
                      n_MCsims = 10,
                      # n_MCsims = ifelse(!missing(data),ceiling(sqrt(nrow(data))),10),
                      runTMLE = c("tmle.intercept", "tmle.covariate"),
                      YnodeDET = NULL,
                      f_gstar2 = NULL,
                      sep = ' ',
                      f_g0 = NULL,
                      h_g0_SummariesModel = NULL,
                      h_gstar_SummariesModel = NULL)
                    ) {

  oldverboseopt <- getOption("tmlenet.verbose")
  options(tmlenet.verbose = verbose)
  gvars$verbose <- verbose
  #----------------------------------------------------------------------------------
  # ADDITIONAL ARGUMENTS (removed from input args of tmlenet())
  #----------------------------------------------------------------------------------
  AnodeDET <- NULL # removed from the input, since this is not implemented
  iid_data_flag <- FALSE  # set to true if no network is provided (usual iid TMLE)
  Q.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction") # NOT USED
  g.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction") # NOT USED
  max_npwt <- 50 # NOT USED YET
  h_logit_sep_k <- FALSE # NOT USED YET
  alpha <- ifelse(is.null(optPars$alpha), 0.05, optPars$alpha)
  lbound <- ifelse(is.null(optPars$lbound), 0.005, optPars$lbound)
  family <- ifelse(is.null(optPars$family), "binomial", optPars$family)
  sep <- ifelse(is.null(optPars$sep), ' ', optPars$sep)
  assert_that(is.character(sep) && length(sep)==1L)
  n_MCsims <- ifelse(is.null(optPars$n_MCsims), 10, optPars$n_MCsims)
  f_g0 <- if(is.null(optPars$f_g0)) {NULL} else {optPars$f_g0}
  if (!is.null(f_g0)) assert_that(is.function(f_g0))
  YnodeDET <- if(is.null(optPars$YnodeDET)) {NULL} else {optPars$YnodeDET}
  f_gstar2 <- if(is.null(optPars$f_gstar2)) {NULL} else {optPars$f_gstar2}
  if (!is.null(f_gstar2)) assert_that(is.function(f_gstar2) || is.vector(f_gstar2))
  # if TRUE, only evaluate the intercept-based TMLE (TMLE_B), if FALSE, evaluate only the covariate-based TMLE (TMLE_A)
  runTMLE <- optPars$runTMLE[1]
  if (is.null(runTMLE) || (runTMLE[1] %in% "tmle.intercept")) {
    onlyTMLE_B <- TRUE
  } else if (runTMLE %in% "tmle.covariate") {
    onlyTMLE_B <- FALSE
  } else {
    stop("optPars[['runTMLE']] argument must be either 'tmle.intercept' or 'tmle.covariate'")
  }

  message("Running tmlenet with the following settings from tmlenet_options(): "); str(gvars$opts)
  message("Running tmlenet with the following settings from optPars arg of tmlenet(): "); str(optPars)

  #----------------------------------------------------------------------------------
  # DETERMINING INTERNAL / EXTERNAL ESTIMATOR NAMES THAT WILL BE EVALUATED
  #----------------------------------------------------------------------------------
  # onlyTMLE_B <- TRUE
  assert_that(assertthat::is.flag(onlyTMLE_B))
  estnames.internal <- c("TMLE_A", "TMLE_B", "h_IPTW", "MLE")
  # estnames.internal <- c("TMLE_A", "TMLE_B", "TMLE_g_IPTW", "h_IPTW", "g_IPTW", "MLE")
  names(estnames.internal) <- estnames.internal
  estnames.out <- c("tmle", "h_iptw", "gcomp")
  if (onlyTMLE_B) {
    estnames.internal <- estnames.internal[-which(estnames.internal%in%"TMLE_A")]
  } else {
    estnames.internal <- estnames.internal[-which(estnames.internal%in%"TMLE_B")]
  }
  names(estnames.out) <- estnames.internal
  estnames.internal <- as.list(estnames.internal)

  #----------------------------------------------------------------------------------
  # MONTE-CARLO SIMULATION PARAMETERS
  #----------------------------------------------------------------------------------
  nQ.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for Q (each of size n)
  ng.MCsims <- as.integer(n_MCsims)  # number of times to sample MC sim for h (each of size n)
  assert_that(is.count(nQ.MCsims))
  assert_that(is.count(ng.MCsims))
  max.err_est <- 0.1    # maximum percent error for MCS estimators

  #----------------------------------------------------------------------------------
  # PARAMETERS FOR ESTIMATING h under g0 & gstar
  #----------------------------------------------------------------------------------
  f.g0 <- f_g0
  f.g0_args <- NULL
  h_user_fcn <- NULL
  h_user <- !(is.null(h_user_fcn))

  #----------------------------------------------------------------------------------
  # Perform some input checks
  # Create an object with data & network information that is passed on to estimation algorithm(s)
  # Build observed sW & sA: (obsdat.sW - a dataset (matrix) of n observed summary measures sW)
  #----------------------------------------------------------------------------------
  if (missing(Ynode)) Ynode <- NULL
  if (is.null(Ynode) && is.null(Qform)) stop("Either Qform or Ynode must be specified")
  if (!is.null(Qform) && !is.null(Ynode)) {
    message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+%
      Ynode)
  }
  if (!is.null(Qform) && is.null(Ynode)) {
    Ynode <- LhsVars(Qform)[1]
    message("Setting the Ynode to: " %+% Ynode)
  }

  if (missing(DatNet.ObsP0)) {
    DatNet.ObsP0 <- eval.summaries(data = data, Kmax = Kmax, sW = sW, sA = sA, 
                                    IDnode = IDnode, NETIDnode = NETIDnode, 
                                    sep = sep, NETIDmat = NETIDmat, verbose = FALSE)$DatNet.ObsP0
  } else {
    data <- DatNet.ObsP0$datnetW$Odata
    Kmax <- DatNet.ObsP0$Kmax
    sW <- DatNet.ObsP0$datnetW$sVar.object
    sA <- DatNet.ObsP0$datnetA$sVar.object
  }

  # variables that will only be available to eval.summaries:
  # node_l <- list(IDnode = IDnode, NETIDnode = NETIDnode, nFnode = nFnode)
  # new version of nodes:
  node_l <- list(nFnode = DatNet.ObsP0$datnetW$nFnode, Anode = Anode, AnodeDET = AnodeDET,
                  Ynode = Ynode, YnodeDET = YnodeDET)
  nobs <- DatNet.ObsP0$nobs
  DatNet.ObsP0$nodes <- node_l
  #----------------------------------------------------------------------------------
  # Defining (and checking) Deterministic Y and A node flags:
  #----------------------------------------------------------------------------------
  if (is.null(AnodeDET)) {
    determ.g <- rep_len(FALSE, nobs)
  } else {
    CheckVarNameExists(data, AnodeDET)
    determ.g <- (data[, AnodeDET] == 1)
  }
  if (is.null(YnodeDET)) {
    determ.Q <- rep_len(FALSE, nobs)
  } else {
    CheckVarNameExists(data, YnodeDET)
    determ.Q <- (data[, YnodeDET] == 1)
  }
  CheckVarNameExists(data, node_l$Anode)
  CheckVarNameExists(data, node_l$Ynode)
  #----------------------------------------------------------------------------------
  # NOTE: YnodeVals = obsYvals, det.Y = determ.Q need to be added to DatNet.ObsP0 after returned its eval.summaries()
  #----------------------------------------------------------------------------------
  obsYvals <- data[,node_l$Ynode]
  DatNet.ObsP0$addYnode(YnodeVals = obsYvals, det.Y = determ.Q)
  #----------------------------------------------------------------------------------
  # (OPTIONAL) ADDING DETERMINISTIC/DEGENERATE Anode FLAG COLUMNS TO sA:
  #----------------------------------------------------------------------------------
    # cancelled adding DET nodes to sVar since all sVar are automatically get added to A ~ predictors + DETnodes...
    # obsdat.sW <- O.datnetW$add_deterministic(Odata = data, userDETcol = "determ.g")$dat.sVar
  # Testing NA for visible det.Y and true observed Y as protected:
  # determ.Q <- c(FALSE, FALSE, FALSE, rep(TRUE, length(determ.Q)-3))
  # length(determ.Q) == length(obsYvals)
  # DatNet.ObsP0 <- DatNet.sWsA$new(datnetW = datnetW, datnetA = datnetA, YnodeVals = obsYvals, det.Y = determ.Q)$make.dat.sWsA()
  # print("DatNet.ObsP0: "); print(DatNet.ObsP0)
  # print(head(cbind(DatNet.ObsP0$noNA.Ynodevals, DatNet.ObsP0$YnodeVals, data[,node_l$Ynode]), 100))

  #----------------------------------------------------------------------------------
  # Optional regressions specs:
  #----------------------------------------------------------------------------------
  Q.sVars <- process_regform(as.formula(Qform), sW.map = c(sW$sVar.names.map, sA$sVar.names.map), sA.map = node_l$Ynode)
  h.g0.sVars <- process_regform(as.formula(hform.g0), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  if (!is.null(hform.gstar)) {
    h.gstar.sVars <- process_regform(as.formula(hform.gstar), sW.map = sW$sVar.names.map, sA.map = sA$sVar.names.map)
  } else {
    h.gstar.sVars <- h.g0.sVars
  }

  if (verbose) {
    print("Input regression Qform (E(Y|sA,sW)): " %+% Qform)
    print("Derived regression Qform (E(Y|sA,sW)):"); str(Q.sVars)
    print("Input regression hform.g0 (P(sA|sW) under g0): " %+% hform.g0)
    print("Derived regression hform.g0 (P(sA|sW) under g0): "); str(h.g0.sVars)
    print("Input regression hform.gstar (P(sA|sW) under g.star): " %+% hform.gstar)
    print("Derived regression hform.gstar (P(sA|sW) under g.star): "); str(h.gstar.sVars)
  }

  #-----------------------------------------------------------
  # Defining and fitting regression for Y ~ sW + sA:
  # todo 45 (m.Q.init) +0: In the future add option to fit separate m.Q.init models for each nF value
  # todo 45 (m.Q.init) +0: Move fitting of Q inside get_all_ests?
  #-----------------------------------------------------------
  check.Qpreds.exist <- unlist(lapply(Q.sVars$predvars, function(PredName) PredName %in% DatNet.ObsP0$names.sVar))
  if (!all(check.Qpreds.exist)) stop("the following predictors in Qform regression could not be located among the summary measures: " %+%
                                    paste0(Q.sVars$predvars[!check.Qpreds.exist], collapse = ","))

  if (verbose) {
    message("================================================================")
    message("fitting E(Y|sA,sW):= ", "P(" %+% node_l$Ynode %+% "=1 | " %+% paste(Q.sVars$predvars, collapse = ",") %+% ")")
    message("================================================================")
  }
  Qreg <- RegressionClass$new(outvar = node_l$Ynode,
                              predvars = Q.sVars$predvars,
                              subset = !determ.Q, ReplMisVal0 = TRUE)
  m.Q.init <- BinOutModel$new(glm = FALSE, reg = Qreg)$fit(data = DatNet.ObsP0)
  # m.Q.init <- BinOutModel$new(glm = FALSE, reg = Qreg)$fit(data = DatNet.ObsP0)$predict(newdata = DatNet.ObsP0)
  # if (verbose) {
  #   print("fit for E(Y|sA,sW) succeeded:")
  #   print("coef(m.Q.init): "); print(coef(m.Q.init))
  # }
  

  # DatNet.ObsP0$YnodeVals       # visible Y's with NA for det.Y
  # DatNet.ObsP0$det.Y           # TRUE/FALSE for deterministic Y's
  # DatNet.ObsP0$noNA.Ynodevals  # actual observed Y's
  # m.Q.init$getoutvarnm         # reg outvar name (Ynode)
  # m.Q.init$getoutvarval        # obsYvals after setting det.Y obs to NA
  # m.Q.init$getprobA1           # predictions (for non-DET Y)
  # m.Q.init$getsubset           # valid subset (!det.Y)
  # m.Q.init$reg                 # regression class (Qreg)
  # Not needed here, using only for comparison:
  # QY.init <- DatNet.ObsP0$noNA.Ynodevals
  # QY.init[!DatNet.ObsP0$det.Y] <- m.Q.init$getprobA1[!DatNet.ObsP0$det.Y] # getting predictions P(Y=1) for non-DET Y

  #----------------------------------------------------------------------------------
  # Create an object with model estimates, data & network information that is passed on to estimation procedure
  #----------------------------------------------------------------------------------
  # 1) define parameters for MC estimation of the substitution estimators
  # 2) define parameters for estimation of the efficient weights h(A^s|W^s)
  est_obj <- list(
                  estnames = estnames.internal,
                  lbound = lbound[1],
                  max.err_eps = max.err_est,  # error tolerance for the mean/var M.C. estimate
                  m.Q.init = m.Q.init,
                  f.g0 = f.g0,
                  sW = sW,
                  sA = sA,
                  Q.sVars = Q.sVars,
                  h.g0.sVars = h.g0.sVars,
                  h.gstar.sVars = h.gstar.sVars,
                  nQ.MCsims = nQ.MCsims,
                  ng.MCsims = ng.MCsims,
                  h_g0_SummariesModel = optPars$h_g0_SummariesModel,
                  h_gstar_SummariesModel = optPars$h_gstar_SummariesModel,
                  # Cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000):
                  max_npwt = max_npwt # NOT IMPLEMENTED  
                  # h_logit_sep_k = h_logit_sep_k, # NOT IMPLEMENTED
                  # h_user = h_user, # NOT IMPLEMENTED
                  # h_user_fcn = h_user_fcn, # NOT IMPLEMENTED
                  )

  est_obj_g1 <- append(est_obj,
                      list(
                        f.gstar = f_gstar1
                        )
                      )

  if (!is.null(f_gstar2)) {
    est_obj_g2 <- append(est_obj,
                      list(
                        f.gstar = f_gstar2
                        )
                      )
  }

  #----------------------------------------------------------------------------------
  # Running MC evaluation for substitution TMLE ests
  #----------------------------------------------------------------------------------
  tmle_g1_out <- get_all_ests(estnames = estnames.internal, DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_obj_g1)
  if (!is.null(f_gstar2)) {
    tmle_g2_out <- get_all_ests(estnames = estnames.internal, DatNet.ObsP0 = DatNet.ObsP0, est_params_list = est_obj_g2)
  } else {
    tmle_g2_out <- NULL
  }

  #----------------------------------------------------------------------------------
  # Create output list (estimates, as. variances, CIs)
  #----------------------------------------------------------------------------------
  EY_gstar1 <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha, DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out = tmle_g1_out)
  EY_gstar2 <- NULL
  ATE <- NULL	
  if (!is.null(f_gstar2)) {
    EY_gstar2 <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha, DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out=tmle_g2_out)
    ATE <- make_EYg_obj(estnames = estnames.internal, estoutnames = estnames.out, alpha = alpha, DatNet.ObsP0 = DatNet.ObsP0, tmle_g_out = tmle_g1_out, tmle_g2_out = tmle_g2_out)
	}

	tmlenet.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)
	class(tmlenet.res) <- c(class(tmlenet.res), "tmlenet")

  options(tmlenet.verbose = oldverboseopt)
  gvars$verbose <- oldverboseopt
	return(tmlenet.res)
}