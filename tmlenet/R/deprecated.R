# nocov start
#----------------------------------------------------------------------------------
# From tmlenet():
#----------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------
  # *** RETIRING g_iptw ESTIMATOR ***
  # Defining and fitting regression for A ~ sW:
  #----------------------------------------------------------------------------------
  # Anode.type <- DatNet.ObsP0$get.sVar.type(node_l$Anode)
  # if (verbose) print("Anode.type: " %+% Anode.type)
  # if (!(Anode.type %in% gvars$sVartypes$bin)) {
  #   message("Anode is not binary, full g_iptw cannot be estimated")
    # m.g0N <- NULL
  # } else {
  #   greg <- RegressionClass$new(outvar = node_l$Anode,
  #                               predvars = g.sVars$predvars,
  #                               subset = !determ.g)
  #   m.g0N <- BinOutModel$new(glm = FALSE, reg = greg)$fit(data = DatNet.ObsP0)
  #   if (verbose) {
  #     print("coef(m.g0N): "); print(coef(m.g0N))
  #   }
  # }
  #----------------------------------------------------------------------------------
  # DEPRECATED... MOVED TO all regs being defined via sW, sA summary measures...
  # Create net_d for fitting m.Q.init, m.g0N and m.h_g0, m.h_gstar
  #----------------------------------------------------------------------------------
  # When Qform.depr is provided, use the formula based fit for Q.init instead of Qreg (sW+sA) fit
  # if (!is.null(Qform.depr)) {
  #   net_d <- cbind(DatNet.ObsP0$dat.sWsA, subset(data, select = node_l$Ynode))
  #   net_d[gvars$misfun(net_d)] <- gvars$misXreplace
  #   m.Q.init.depr <- f_est(net_d[!determ.Q,], Qform.depr, family = family)
  #   QY.init.depr <- data[, node_l$Ynode] # setting deterministic node values
  #   QY.init.depr[!determ.Q] <- predict(m.Q.init.depr, newdata = net_d[!determ.Q,], type = "response") # predict p(Y) for non determ nodes    
  #   if (is.null(gform.depr)) {
  #     gform.depr <- node_l$Anode %+% " ~ " %+% paste0(DatNet.ObsP0$datnetW$names.sVar, collapse="+") # default to main terms in DatNet.ObsP0$datnetW
  #   }
  #   m.g0N.depr <- f_est(net_d[!determ.g,], gform.depr, family = family) # Set A = 0 when determ.g == 1
  #   d_sel <- cbind(d_sel, QY.init = QY.init.depr) # (DEPRECATED, TO BE REMOVED)
  #   # if (verbose) {
  #     print("head(net_d)"); print(head(net_d, 5))      
  #     print("new coef(m.Q.init): "); print(coef(m.Q.init))
  #     print("old coef(m.Q.init.depr): "); print(coef(m.Q.init.depr))
  #     print("coef(m.g0N.depr)"); print(coef(m.g0N.depr))
  #     print("head(d_sel) old: "); print(head(d_sel))
  #     message("Running tmlenet with... ");
  #     message("Qform.depr: " %+% Qform.depr)
  #     message("gform.depr: " %+% gform.depr)
  #     message("hform.depr: " %+% hform.depr)
  #   # }
  # }
  # dfcheck <- data.frame(QY.init = QY.init, QY.init.depr = QY.init.depr, diff = QY.init - QY.init.depr)
  # head(dfcheck, 50)
  # browser()
  # stop()  
  #----------------------------------------------------------------------------------
  # DEPRECATED: Run TMLE univariate fluctuations for each g.star and/or ATE:
  #----------------------------------------------------------------------------------
  # if (!is.null(Qform.depr)) {
  #   est_obj_g1$m.g0N <- m.g0N.depr
  #   est_obj_g1$m.Q.init <- m.Q.init.depr
  #   tmle_g1_out.depr <- get_all_ests.old(data = d_sel, est_obj = est_obj_g1)
  #   tmle_g2_out.depr <- NULL
  #   if (!is.null(f_gstar2)) {
  #     est_obj_g2$m.g0N <- m.g0N.depr
  #     est_obj_g2$m.Q.init <- m.Q.init.depr
  #     tmle_g2_out.depr <- get_all_ests.old(data = d_sel, est_obj = est_obj_g2)
  #   }
  # }

#-----------------------------------------------------------------------------
# Formula based glm model fit 
#-----------------------------------------------------------------------------
f_est <- function(d, form, family) {
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
    SuppressGivenWarnings({
              m <- glm(as.formula(form),
                  data = d,
                  family = family,
                  control = ctrl)
              },
              GetWarningsToSuppress())
    return(m)
}


#----------------------------------------------------------------------------------
# General S3 summary functions:
#----------------------------------------------------------------------------------

# Get summary measures for one or two tmlenet objects 
#(SEs, p-values, CIs)
# If two objects, include effect measures (additive effect, relative risk, odds ratio)
# summary.tmlenet <- function(object, control.object=NULL, 
            # estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  # #object is treatment, control.object is control
  # if (! is.null(control.object) && class(control.object) != "tmlenet") 
        # stop("the control.object argument to summary.tmlenet must be of class tmlenet")
  # if (! estimator %in% c("tmle", "iptw", "gcomp", "naive")) stop("estimator should be one of: tmle, iptw, gcomp, naive")
  # treatment.summary <- NULL
  # if (! is.null(control.object)) {
    # control.summary <- GetSummary(control.object$estimates[estimator], control.object$IC[[estimator]], loggedIC=FALSE)
    # effect.measures <- GetEffectMeasures(est0=control.object$estimates[estimator], 
                                          # IC0=control.object$IC[[estimator]], 
                                          # est1=object$estimates[estimator], 
                                          # IC1=object$IC[[estimator]])
    # effect.measures.summary <- lapply(effect.measures, function (x) GetSummary(x$est, x$IC, x$loggedIC))
  # } else {
    # control.summary <- effect.measures.summary <- NULL
  # }
  # ans <- list(treatment=treatment.summary, control=control.summary, 
              # effect.measures=effect.measures.summary, treatment.call=object$call, 
              # control.call=control.object$call, estimator=estimator)
  # class(ans) <- "summary.tmlenet"
  # return(ans)
# }

# # Get summary measures for tmlenet parameters (standard errors, p-values, confidence intervals)
# summary.tmlenet <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) 
# {
#   if (! estimator %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp")  
#   n <- nrow(IC)
#   v <- apply(IC, 2, var)
#   std.dev <- sqrt(v/n)
#   pval <- 2 * pnorm(-abs(estimate / std.dev))
#   CI <- GetCI(estimate, std.dev)
#   cmat <- cbind(estimate, std.dev, CI, pval)
#   dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", "CI 2.5%", "CI 97.5%", "p-value"))
#   ans <- list(cmat=cmat, estimator=estimator)
#   class(ans) <- "summary.tmlenet"
#   return(ans)
# }


# # Print method for summary.tmlenet
# print.summary.tmlenet <- function(x, ...) {
#   cat("Estimator: ", x$estimator, "\n")
#   if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
#   if (is.null(x$control)) {
#     PrintCall(x$treatment.call)
#     PrintSummary(x$treatment)
#   } else {
#     cat("\nControl ")
#     PrintCall(x$control.call)
#     cat("Treatment ")
#     PrintCall(x$treatment.call)
#     cat("Treatment Estimate:\n")
#     PrintSummary(x$treatment)
#     cat("\nControl Estimate:\n")
#     PrintSummary(x$control)
#     cat("\nAdditive Effect:\n")
#     PrintSummary(x$effect.measure$ATE)
#   }
#   invisible(x)
# }

# # Print method for tmlenet
# print.tmlenet <- function(x, ...) {
#   PrintCall(x$call)
#   cat("TMLE Estimate: ", x$estimates["tmle"], "\n")
#   invisible(x)
# }

# # Print a call
# PrintCall <- function(cl) {
#   cat("Call:  ", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
# }

# Print estimate, standard error, p-value, confidence interval
# PrintSummary <- function(x) {
#   cat("   Parameter Estimate: ", signif(x$estimate, 5))
#   cat("\n    Estimated Std Err: ", signif(x$std.dev^2, 5))
#   cat("\n              p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue, 5)))
#   cat("\n    95% Conf Interval:",paste("(", signif(x$CI[1], 5), ", ", signif(x$CI[2], 5), ")", sep=""),"\n")
#   invisible(x)
# }

# # Calculate estimate, SE, p-value, confidence interval
# GetSummary <- function(estimate, IC, loggedIC) {
#   if (is.null(IC)) {
#     std.dev <- NA
#   } else {
#     n <- length(IC)
#     std.dev <- sqrt(var(IC) / n)
#   }
#   if (loggedIC) {
#     pvalue <- 2 * pnorm(-abs(log(estimate) / std.dev))
#     CI <- exp(GetCI(log(estimate), std.dev))
#   } else {
#     pvalue <- 2 * pnorm(-abs(estimate / std.dev))
#     CI <- GetCI(estimate, std.dev)
#   }
  
#   return(list(estimate=estimate, std.dev=std.dev, pvalue=pvalue, CI=CI))
# }

# # Calculate 95% confidence interval
# GetCI <- function(estimate, std.dev) {
#   x <- qnorm(0.975) * std.dev
#   CI <- cbind("2.5%"=estimate - x, "97.5%"=estimate + x)
#   return(CI)
# }

# # Calculate Average Treatment Effect
# GetEffectMeasures <- function(est0, IC0, est1, IC1) {  
#   names(est0) <- names(est1) <- NULL
#   ATE <- est1 - est0
#   if (is.null(IC0)) {
#     ATE.IC <- NULL
#   } else {
#     ATE.IC <- -IC0 + IC1
#   }
#   return(list(ATE=list(est=ATE, IC=ATE.IC, loggedIC=FALSE)))
# }


# Run GLM or SuperLearner (in future vs) for fitting initial Q and g
# old call: (data, form, family)
# Estimate <- function(form, d, subs, family, newdata, SL.library) {
#   f <- as.formula(form)
#   if (any(is.na(d[subs, LhsVars(f)]))) stop("NA in Estimate")
#   # stats <- CalcRegressionStats(d, f, subs)
#   if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
#     #estimate using GLM
#     if (sum(subs) > 1) {
#       SuppressGivenWarnings({
#         m <- glm(f, data=d, subset=subs, family=family, control=glm.control(trace=FALSE, maxit=1000))
#         predicted.values <- predict(m, newdata=newdata, type="response")
#       }, GetWarningsToSuppress())
#     } else {
#       #glm breaks when sum(subs) == 1
#       predicted.values <- rep(d[subs, LhsVars(f)], nrow(newdata))
#     }
#   } else {
#     #estimate using SuperLearner
#     if (family == "quasibinomial") family <- "binomial"
#     #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
#     rhs <- setdiff(RhsVars(f), rownames(alias(f, data=d[subs,])$Complete))  
#     #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
#     new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))
    
#     m <- SuperLearner(Y=d[subs, LhsVars(f)], X=d[subs, rhs, drop=FALSE], SL.library=SL.library, 
#                       verbose=FALSE, family=family, newX=newdata[new.subs, rhs, drop=FALSE])
#     predicted.values <- rep(NA, nrow(newdata))
#     predicted.values[new.subs] <- m$SL.predict
#   }
#   return(list(values=predicted.values, model=m))
# }

# (DEPRECATED)
#-----------------------------------------------------------------------------
# USE glm.fit FUNCTION FOR FASTER FITTING of LOGISTIC REG TAKES DESIGN MAT AND Y VECTOR
#-----------------------------------------------------------------------------
.f.est_binom_fast <- function(X_mat, Y_vals) {
  ctrl <- glm.control(trace=FALSE, maxit=1000)          
    SuppressGivenWarnings({
              m.fit <- glm.fit(x=cbind(1,X_mat), y=Y_vals, family = binomial(), control=ctrl)
              }, GetWarningsToSuppress())
    return(m.fit)
}
.f_predict_fast <- function(glmfit, new_mtx) {
      new_mtx <- cbind(1,new_mtx)
      eta <- new_mtx[,!is.na(glmfit$coef), drop=FALSE] %*% glmfit$coef[!is.na(glmfit$coef)]
      return(glmfit$family$linkinv(eta))
}

# (DEPRECATED) (WAS USED FOR PLUG-IN NPMLE ESTIMATOR OF \bar{h}_gN)
# sample treatments with probA from above fcn
.f.gen.A_N <- function(df, deterministic, m.gN) {
	n <- nrow(df)
  	rbinom(n, 1, .f.gen.probA_N(df, deterministic, m.gN))
}


# (DEPRECATED)
#-----------------------------------------------------------------------------
# Calculate the joint probability of observing a=(a_1,..,a_k) from P(A^s=1)
# takes the matrix of predicted probabilities P(A^s=1): data.probA
# and a matrix of observed a^s (a_1,..,a_k): data.indA
#-----------------------------------------------------------------------------
.f.cumprod.matrix <- function(data.indA, data.probA) {
  y <- matrix(1, nrow=dim(data.probA)[1], ncol=dim(data.probA)[2])
  y[, 1] <- data.probA[,1]^as.integer(data.indA[,1]) *
              (1-data.probA[,1])^(1-as.integer(data.indA[,1]))
  if (dim(data.probA)[2] > 1) {
    for (i in 2:dim(data.probA)[2]) {
      y[,i] <- y[,i-1] * (data.probA[,i]^as.integer(data.indA[,i]) *
                (1-data.probA[,i])^(1-as.integer(data.indA[,i])))
    }
  }
  # return(round(y[,dim(data.probA)[2]],6))
  return(y[,dim(data.probA)[2]])
}

f.gen.probA.star.old <- function(k, df_AllW, fcn_name, f_args = NULL) {
  .f_g_wrapper <- function(k, df_AllW, fcn_name, ...) {
      args0 <- list(k = k, data = df_AllW)
      args <- c(args0, ...)
    do.call(fcn_name, args)
  }
  probA <- .f_g_wrapper(k, df_AllW, fcn_name, f_args)
  return(probA)
}


# (DEPRECATED) No longer called
f.gen.A.star.old <- function(k, df_AllW, fcn_name, f_args = NULL) {
  n <- nrow(df_AllW)
  rbinom(n, 1, f.gen.probA.star.old(k, df_AllW, fcn_name, f_args))
}

# (DEPRECATED)
# get the prob of A (under g_0) using fit g_N, estimated regression model for g0;
.f.gen.probA_N <- function(df, deterministic, m.gN) {
    SuppressGivenWarnings({
      g_N <- predict(m.gN, newdata=df, type="response")
          },
          GetWarningsToSuppress())
    #************************************************
    g_N[deterministic] <- 0
    #************************************************
    return(g_N)
}

#------------------------------------------------------------------------------
# (DISABLED, RETIRING)
# IPTW ESTIMATOR (est Y_g_star based on weights g_star(A|W)/g_N(A|W) )
#------------------------------------------------------------------------------
iptw_est <- function(k, data, node_l, m.gN, f.gstar, f.g_args, family="binomial", NetInd_k, lbound=0, max_npwt=50, f.g0=NULL, f.g0_args=NULL, iidIPTW=FALSE) {
  n <- nrow(data)
  netW <- NULL
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
  for (Wnode in node_l$Wnodes) {
    netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)))
  }
  cA.mtx <- cbind(netW, subset(data, select=nFnode))
  indA <- data.frame(.f.allCovars(k, NetInd_k, data[, Anode], Anode))

  determ.g <- data$determ.g
  determ.g_vals <- data[determ.g, Anode]
  # predict g*(A=1|W):
  pA_gstar <- f.gen.probA.star.old(k, cA.mtx, f.gstar, f.g_args)
  netpA_gstar <- .f.allCovars(k, NetInd_k, pA_gstar, Anode)

  # calculate likelihoods P_*(A=a|W):
  if (iidIPTW) { # for iid IPTW, use only A_i, no A_j, for j\in F_i:
    indA <- indA[, Anode, drop=FALSE]
    netpA_gstar <- netpA_gstar[, Anode, drop=FALSE]
  }

  gstar_A <- .f.cumprod.matrix(indA, netpA_gstar) # calculate likelihoods P_0(A=a|W) under g_star:
  #************************************************
  if (is.null(f.g0)) { #If g_0 is unknown, use logistic fit m.gN
    # print("RUNNING IPTW ON FITTED g0_N"); print(m.gN)
    pA_g0N <- .f.gen.probA_N(cA.mtx, determ.g, m.gN)
  }
  else {   # If g_0 is known get true P_g0
    # print("RUNNING IPTW ON TRUE g0")
    pA_g0N <- f.gen.probA.star.old(k, cA.mtx, f.g0, f.g0_args)
  }
  #************************************************
  netpA_g0 <- .f.allCovars(k, NetInd_k, pA_g0N, Anode)

  # for iid IPTW, use only A_i, no A_j, for j\in F_i:
  if (iidIPTW) {
    indA <- indA[,Anode,drop=FALSE]
    netpA_g0 <- netpA_g0[,Anode,drop=FALSE]
  }
  # print("head(indA)"); print(head(indA))
  # print("head(netpA_g0)"); print(head(netpA_g0))

  # calculate likelihoods P_0(A=a|W):
  g0_A <- .f.cumprod.matrix(indA, netpA_g0)
  ipweights <- gstar_A / g0_A
  ipweights[is.nan(ipweights)] <- 0     # 0/0 detection
  # ipweights <- bound(ipweights, c(0,10^6))
  # lower bound g0 by lbound
  ipweights <- bound(ipweights, c(0,1/lbound))

  # scale weights by total contribution (instead of bounding by consts):
  # cap the prop weights scaled at max_npwt (for =50 -> translates to max 10% of total weight for n=500 and 5% for n=1000)
  # ipweights <- scale_bound(ipweights, max_npwt, n)
  # print("iptw wts range after bound"); print(summary(ipweights))
  return(ipweights)
}


#-----------------------------------------------------------------------------
# (DEPRECATED)
# return entire network matrix from indiv. covar (Var) + covariate itself as first column
# NetInd_k is a matrix (N x k) of network friend indicies; 
# When obs i doesn't have j friend, NetInd[i, j] = NA
#-----------------------------------------------------------------------------
.f.allCovars <- function(k, NetInd_k, Var, VarNm, misval = gvars$misXreplace) {
# .f.allCovars <- function(k, NetInd_k, Var, VarNm, misval = 0L) {
  assertthat::assert_that(is.matrix(NetInd_k))
  n <- length(Var)
  netVar_names <- NULL
  netVar_full <- NULL
  d <- matrix(0L, nrow = n, ncol = k + 1)
  d[, 1] <- Var
  d[, c(2 : (k+1))] <- apply(NetInd_k, 2, function(k_indx) {
                                          netVar <- Var[k_indx]  # netVar values for non-existing friends are automatically NA
                                          netVar[is.na(netVar)] <- misval
                                          return(netVar)
                                          })
  if (k>0) netVar_names <- netvar(varnm = VarNm, fidx = c(1:k))
  Var_names <- c(VarNm, netVar_names)
  colnames(d) <- Var_names
  return(d)
}

# (DEPRECATED)
#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
#---------------------------------------------------------------------------------
# METHOD 1 (empirical distribution of \bar{h}=\sum{h_i}) - NO LONGER USED
# SEE OLDER git repo for this implementation
  # For given data, estimate g[A|cA] and calculate est. of h_i(c) for given value of c and each i. 
  # * Draw B samples from emp distribution of cY, for each i;
  # * Assume the same dimensionality for cY acorss all i;
  # * Assume W's are iid, use g_N(A|W) and keep N_i constant;
  # * Drawing from distribution of g_N or g* to get A's;
  # * Calculate h_bar from the joint distribution of all c_i^Y over all i;
  
# METHOD 2 (fit logit to each P(A_i|A's,W's)) 
# Predict P(A_1=a1|W1,..,Wk)*P(A_2=a2|A_1,W1,..,Wk)*...*P(A_k=a_k|A_1,...,A_k,W1,...,Wk)
pred.hbars.old <- function(new_data=NULL, fit_h_reg_obj, NetInd_k) {
   pred_h_fcn <- function(h_logit_sep_k, m.gAi_vec) {
      .f_pred_h_k <- function(k, sel_k_indx, m.gAi_vec) {
        k_sel <- sel_k_indx
        A_nms_arr <- colnames(indA)[c(1:(k+1))]
        indA <- indA[,c(1:(k+1)), drop=FALSE]

        if (is.null(hform)) {
          W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
          W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]         
        } else {
          W_nms_arr <- all.vars(as.formula(hform))[-1]
        }
        probAi_vec <- NULL
        for (k_i in (1:(k+1))) {
          Covars_nms <- c(A_nms_arr[-c(k_i:(k+1))], W_nms_arr)
          deterministic <- (determ_cols[,k_i] == 1L)
          predAs_idx <- (!deterministic & k_sel)
          if (!fit_fastmtx) {
            #************************************************                       
            probAi <- .f.gen.probA_N(data.frame(cY_mtx), !(predAs_idx), m.gAi_vec[[k_i]])
            #************************************************           
          }
          else {
            newdat_mtx <- cY_mtx[predAs_idx, Covars_nms, drop=FALSE]
            #************************************************           
            probAi <- rep_len(0,n)
            #************************************************
            if (sum(predAs_idx)>0) probAi[predAs_idx] <- .f_predict_fast(m.gAi_vec[[k_i]], newdat_mtx)
          }
          probAi_vec <- cbind(probAi_vec, probAi)
        }
        # return the likelihood for the entire vector of Ai's, based on cY_mtx (cY)
        return(list(h_vec = .f.cumprod.matrix(indA, probAi_vec)[k_sel]))
      }
      if (h_logit_sep_k) {
        #----------------------------------------------------------------------
        # VERSION 1: predict logit separately for each k in the network and combine results 
        h_vec <- rep_len(0,n)
        for (nFriends in (0:k)) {
          k_sel <- (new_data[,node_l$nFnode]==nFriends)
          # nothing to predict if no one observed with this size ntwk
          if (sum(k_sel)!=0) {
            h_vec[k_sel] <- .f_pred_h_k(k=nFriends, sel_k_indx=k_sel, m.gAi_vec[[nFriends+1]])$h_vec
          }
        }
        return(list(h_vec=h_vec))
      } else {
        #---------------------------------------------------------------------
        # VERSION 2: predict logistic to entire dataset with nFriends as an extra covariate
        return(.f_pred_h_k(k=k, sel_k_indx=rep_len(TRUE,n), m.gAi_vec))
      }
    }
    # ASSIGN VARIABLE NAMES BASED ON ARGUMENT fit_h_reg_obj
    k <- fit_h_reg_obj$k
    h_user <- fit_h_reg_obj$h_user
    h_user_fcn <- fit_h_reg_obj$h_user_fcn
    lbound <- fit_h_reg_obj$lbound
    max_npwt <- fit_h_reg_obj$max_npwt
    h_logit_sep_k <- fit_h_reg_obj$h_logit_sep_k
    node_l <- fit_h_reg_obj$node_l
    gform <- fit_h_reg_obj$gform
    hform <- fit_h_reg_obj$hform
    fit_fastmtx <- fit_h_reg_obj$fit_fastmtx
    netW_namesl <- fit_h_reg_obj$netW_namesl
    netA_names <-  fit_h_reg_obj$netA_names
    determ_cols_Friend <- fit_h_reg_obj$determ_cols_Friend
    # select only baseline covariates (W's) that are included in the original gform:

    if (!is.null(new_data)) {
      determ.g_user <- new_data$determ.g
      determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
      determ_cols <- (determ_cols_user | determ_cols_Friend)
    }
    if (is.null(new_data)) {    # use original fitted data for prediction
      new_data <- fit_h_reg_obj$cY_mtx_fitted      
      determ_cols <- fit_h_reg_obj$determ_cols_fitted
    }
    n <- nrow(new_data)
    indA <- as.matrix(new_data[, netA_names])
    if (is.null(hform)) {
      W_nms <- unlist(netW_namesl)
    } else {
      W_nms <- all.vars(as.formula(hform))[-1]
    }
    if (!(node_l$nFnode%in%W_nms)) { W_nms <- c(W_nms, node_l$nFnode) }
    cY_mtx <- cbind(determ_cols, as.matrix(new_data[, c(node_l$nFnode, W_nms, netA_names)]))
    #---------------------------------------------------------------------
    # MAIN BODY OF THE FUNCTION
    if (h_user==FALSE) {
      P.hbar.c <- pred_h_fcn(h_logit_sep_k, fit_h_reg_obj$m.gAi_vec_g)
      P.hbar.star.c <- pred_h_fcn(h_logit_sep_k, fit_h_reg_obj$m.gAi_vec_gstar)
    }
    else {
      h_bars_user <- h_user_fcn(k, new_data, node_l, NetInd_k)
      P.hbar.c <- h_bars_user$P.hbar.c
      P.hbar.star.c <- h_bars_user$P.hbar.star.c
    }

    h_tilde <- P.hbar.star.c$h_vec / P.hbar.c$h_vec
    h_tilde[is.nan(h_tilde)] <- 0     # 0/0 detection
    h_tilde <- bound(h_tilde, c(0,1/lbound))
    df_h_bar_vals <- data.frame(cY.ID = 0, 
                                h.star_c = P.hbar.star.c$h_vec, 
                                h_c = P.hbar.c$h_vec, 
                                h = h_tilde
                                )
    return(list(df_h_bar_vals=df_h_bar_vals))
}


# (DEPRECATED)
# fit models for m_gAi
#---------------------------------------------------------------------------------
fit.hbars.old <- function(data, h_fit_params) {
    #---------------------------------------------------------------------------------
    # 1) return observed network data cY_mtx if is.null(f.g_name)
    # 2) pass only one netW, which can be separate for g_star and g_0
    # SAMPLE A LARGE DATASET of cY's for given functions g* or g0, of size p*N for some p
    # Get rid of the loop, by assigning .f.allCovars ALL AT ONCE
    # NOTE (OS 06/02/15): The only reason we need to pass netW and netW_full is because g_0 and g^*
    # are assumed to be based on the same sW! This shouldn't be the case, need to allow sW_g ad sW_gstar to be different
    #---------------------------------------------------------------------------------
    get_hfit_data <- function(cY_mtx, k, Anode, NetInd_k, netW, f.g_name, f.g_args, p, misval = 0L)  {
      samp_g_data <- function(df_sel) {
        resamp_A <- f.gen.A.star.old(k, df_sel, f.g_name, f.g_args)
        resamp_netA <- .f.allCovars(k, NetInd_k, resamp_A, Anode, misval = misval)
        fit.g_data <- cbind(netW, resamp_netA)
        return(fit.g_data)
      }
      if (is.null(f.g_name)) { return(cY_mtx) }
      n <- nrow(netW)
      df_samp_g <- samp_g_data(netW_full)  # big resampled matrix of c's (repeated data p times)
      fit.g_data_large <- matrix(nrow=(n*p), ncol=ncol(df_samp_g))
      colnames(fit.g_data_large) <- colnames(df_samp_g)
      fit.g_data_large[c(1:n),] <- as.matrix(df_samp_g)
      for (i in (2:p)) {
        fit.g_data_large[c(((i - 1) * n + 1):(n * i)), ] <- as.matrix(samp_g_data(netW_full))
      }
      return(fit.g_data_large)
    }
    # defining the vector of c^A's that needs evaluation under h(c) 
    .f.mkstrNet <- function(Net) apply(Net, 1, function(Net_i) paste(Net_i, collapse=" "))
    #---------------------------------------------------------------------
    # Estimate h with logistic loss fcn (regression-based)
    # Estimate P(A_1|W1,..,Wk)*P(A_2|A_1,W1,..,Wk)*..*P(A_k|A_1,...,A_k,W1,...,Wk)
    # fitting each component P(A_i) as a logistic regression, for g_0 and g*
    # Does not rely on g_N - fitted model for g_0!!!!
    # Saves the estimated models for \bar{h}
    #---------------------------------------------------------------------  
    fit_h_fcn <- function(fit.g_data, p, h_logit_sep_k, netW_namesl, indA) {
      #----------------------------------------------------------------------
      # VERSION 1: Fit logit separately for each k in the network and combine together 
      # call logistic loss est. of h for each nFriends=0,..,k separately
      # defined k and k_sel values
      # .... REMOVED ...
      #---------------------------------------------------------------------
      # VERSION 2 for logistic loss h:
      # add extra smoothing: fit logistic to entire dataset with nFriends as an extra covariate
      # save the estimated model for \bar{h}
      k_sel <- rep_len(TRUE, n)
      #---------------------------------------------------------------------
      A_nms_arr <- colnames(indA)[c(1:(k+1))]
      indA <- indA[,c(1:(k+1)), drop=FALSE]
      if (is.null(hform)) {
        stop("NOT IMPLEMENTED, SUPPLY hform")
        # W_nms_arr <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
        # W_nms_arr <- W_nms_arr[!is.na(W_nms_arr)]
      } else {
        W_nms_arr <- all.vars(as.formula(hform))[-1]
      }
      probAi_vec <- NULL
      m.gAi_vec <- NULL
      for (k_i in (1:(k+1))) {  # factorization by nFriends (k) value
        A_i_nm <- A_nms_arr[k_i] # A variable we are predicting
        Covars_nms <- c(A_nms_arr[-c(k_i:(k+1))], W_nms_arr) # dependent covars
        deterministic <- (determ_cols[,k_i] == 1L)   # fit only to non-deterministic trt nodes
        fitAs_idx <- (!deterministic & k_sel)
        # USE glm.fit FUNCTION FOR FASTER FITTING (using design matrix format)
        Y_vals <- fit.g_data[rep.int(fitAs_idx, p), A_i_nm]
        X_mat <- as.matrix(fit.g_data[rep.int(fitAs_idx, p), Covars_nms, drop=FALSE], byrow=FALSE)
        m.gAi <- .f.est_binom_fast(X_mat, Y_vals)
        newdat_mtx <- cY_mtx[fitAs_idx, Covars_nms, drop=FALSE]   # original data to predict A

        # ************************************************
        probAi <- rep_len(0,n)
        # ************************************************

        if (sum(fitAs_idx)>0) probAi[fitAs_idx] <- .f_predict_fast(m.gAi, newdat_mtx)
        m.gAi_vec <- c(m.gAi_vec, list(m.gAi))
        probAi_vec <- cbind(probAi_vec, probAi)
      }
      # likelihood for the entire vector of Ai's, based on cY_mtx (cY)
      h_vec <- .f.cumprod.matrix(indA, probAi_vec)[k_sel]
      return(list(h_vec=h_vec, m.gAi_vec=m.gAi_vec))
    }

    #---------------------------------------------------------------------------------
    # PARAMETERS FOR LOGISTIC ESTIMATION OF h
    #---------------------------------------------------------------------------------
    # replace with p adaptive to k: p <- 100*(2^k)
    n <- nrow(data)
    n_samp_g0gstar <- h_fit_params$n_samp_g0gstar
    family <- h_fit_params$family
    k <- h_fit_params$Kmax
    node_l <- h_fit_params$node_l
    NetInd_k <- h_fit_params$NetInd_k
    lbound <- h_fit_params$lbound
    max_npwt <- h_fit_params$max_npwt

    h_logit_sep_k=h_fit_params$h_logit_sep_k
    h_user=h_fit_params$h_user; h_user_fcn=h_fit_params$h_user_fcn; hform=h_fit_params$hform
    # m.gN=h_fit_params$m.g0N # not used
    f.gstar=h_fit_params$f.gstar; f.g_args=h_fit_params$f.g_args
    f.g0=h_fit_params$f.g0; f.g0_args=h_fit_params$f.g0_args
    fit_fastmtx <- TRUE
    gform <- h_fit_params$gform
    #---------------------------------------------------------------------------------

    # select only baseline covariates (W and W_net) that are included in hform
    nFnode <- node_l$nFnode
    Anode <- node_l$Anode

    # W's aren't resampled (assumed independent but not iid) => get all network W's from the original data
    netW <- NULL
    netW_namesl <- list()
    hform_covars <- all.vars(as.formula(hform))[-1]
    for (Wnode in node_l$Wnodes) {
      netW_names <- NULL
      # New condition (02/16/14):
      # only add netW_i covars for bsl covariate W if W doesn't already have "net" in its name
      netWnm_true <- (length(agrep("netF", Wnode, max.distance=list(insertions=0, substitutions=0)))==1)
      if (k>0 & !netWnm_true) netW_names <- netvar(Wnode, c(1:k))
      allW_names <- c(Wnode, netW_names)
      # only include bsl covariates (W, netW) if model for g_N depends on any of them
      if (!netWnm_true) { # only get netW's for non-network covars:
        netW_i <- .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)
        netW <- cbind(netW, netW_i)
        if (any(allW_names %in% hform_covars)) netW_namesl <- c(netW_namesl, list(colnames(netW_i)))
      } else {
        netW <- cbind(netW, as.matrix(subset(data, select=Wnode)))
        if (any(allW_names %in% hform_covars)) netW_namesl <- c(netW_namesl, list(Wnode))
      }
    }

    # OS 06/01/15: moved from fit_h_fcn():
    if (is.null(hform)) {
      W_nms <- unlist(lapply(netW_namesl, function(netWi) netWi[c(1:(k+1))]))
      W_nms <- W_nms[!is.na(W_nms)]
    } else {
      W_nms <- all.vars(as.formula(hform))[-1]
    }
    if (!(nFnode%in%W_nms)) { W_nms <- c(W_nms,nFnode) }

    # 02/10/14: was a bug: get all W's for generating g_star (otherwise can get an error)
    netW_full <- NULL
    for (Wnode in node_l$Wnodes) {
        netW_i <- .f.allCovars(k, NetInd_k, data[, Wnode], Wnode)
        netW_full <- cbind(netW_full, netW_i)
    }
    netW_full <- cbind(netW_full, as.matrix(subset(data, select = nFnode)))

    #-----------------------------------------------------------
    # ADDING DETERMINISTIC FLAG COLUMNS TO netW
    #-----------------------------------------------------------
    # get deterministic As for the entire network of each unit (set by user)
    determ.g_user <- data$determ.g
    # determ.gvals_user <- data[,Anode] # add the values for deterministic nodes as well
    determ_cols_user <- .f.allCovars(k, NetInd_k, determ.g_user, "determ.g_true")
    determ_cols_Friend <- 1L - .f.allCovars(k, NetInd_k, rep_len(1L, n), "determ.g_Friend")
    determ_cols <- (determ_cols_user | determ_cols_Friend)
    print("head(determ_cols_Friend)"); print(head(determ_cols_Friend))
    print("head(determ_cols)"); print(head(determ_cols))

    #-----------------------------------------------------------
    # DEFINING LOGISTIC REGs AND THE SUBSETING EXPRESSIONS
    # (1 expression per 1 regression P(sA[j]|sA[j-1:0], sW))
    #-----------------------------------------------------------
    A_nms <- netvar(Anode, c(0:k))
    regs_idx <- seq_along(A_nms)-1
    indA <- .f.allCovars(k = k, NetInd_k = NetInd_k, Var = data[,Anode], VarNm = Anode, misval = 0L)
    netW <- cbind(determ_cols, netW, as.matrix(subset(data, select=nFnode)))
    cY_mtx <- cbind(netW, indA)
    #-----------------------------------------------------------
    # BODY OF MAIN FUNCTION:
    #-----------------------------------------------------------
    ##########################################
    message("fitting h under g_0...")
    ##########################################
    p_h0 <- ifelse(is.null(f.g0), 1, n_samp_g0gstar)
    fit.g0_dat <- get_hfit_data(cY_mtx = cY_mtx, k = k, Anode = Anode, NetInd_k = NetInd_k,
                                netW = netW, f.g_name = f.g0, f.g_args = f.g0_args, p = p_h0)
    fit.g0_dat <- data.frame(fit.g0_dat)
    P.hbar.c <- fit_h_fcn(fit.g_data = fit.g0_dat, p = p_h0, h_logit_sep_k = h_logit_sep_k, netW_namesl = netW_namesl, indA = indA)
    ##########################################
    message("fitting h under g_star...")
    ##########################################
    # produces a matrix (not df, so needs to be converted for faster glm.fit)
    t1 <- system.time(
      fit.gstar_dat <- get_hfit_data(k = k, Anode = Anode, NetInd_k = NetInd_k, 
                                      netW = netW, f.g_name = f.gstar, f.g_args = f.g_args,
                                      p = n_samp_g0gstar)
    )
    # determ_cols may change under gstar, e.g., when P_g^*(A=1|W)=1 for some W
    fit.gstar_dat <- data.frame(fit.gstar_dat)
    P.hbar.star.c <- fit_h_fcn(fit.g_data = fit.gstar_dat, p = n_samp_g0gstar, h_logit_sep_k = h_logit_sep_k, netW_namesl = netW_namesl, indA = indA)
    # ##########################################
    # 3) Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
    ##########################################
    h_tilde <- P.hbar.star.c$h_vec / P.hbar.c$h_vec
    h_tilde[is.nan(h_tilde)] <- 0     # 0/0 detection
    h_tilde <- bound(h_tilde, c(0,1/lbound))
    df_h_bar_vals <- data.frame(cY.ID = .f.mkstrNet(cY_mtx),
                                h.star_c = P.hbar.star.c$h_vec,
                                h_c = P.hbar.c$h_vec,
                                h = h_tilde
                                )

    fit_h_reg_obj <- list(m.gAi_vec_g = P.hbar.c$m.gAi_vec,
                          m.gAi_vec_gstar = P.hbar.star.c$m.gAi_vec,
                          k = k,
                          lbound = lbound,
                          h_logit_sep_k = h_logit_sep_k,
                          h_user = h_user,
                          h_user_fcn = h_user_fcn,
                          fit_fastmtx = fit_fastmtx,
                          node_l = node_l,
                          netW_namesl=netW_namesl,
                          gform = gform,
                          hform = hform,
                          netA_names=colnames(indA),
                          determ_cols_Friend = determ_cols_Friend,
                          determ_cols_fitted = determ_cols,
                          cY_mtx_fitted = cY_mtx)
    return(list(df_h_bar_vals = df_h_bar_vals, fit_h_reg_obj = fit_h_reg_obj))
}


#---------------------------------------------------------------------------------
# G-Comp & TMLEs: Use Monte-Carlo to estimate psi under stochastic g^* 
#---------------------------------------------------------------------------------
  # For given data, take Q[Y|cY]=m.Q.init and calcualte est. of psi under g*=f.gstar using Monte-Carlo integration:
  # * W_i can be iid or not (in latter case W are not resampled);
  # * Draw from the distributions of W and g*(A|W), keeping N_i constant, recalculate cY and cA each time;
  # * Recalculate Y^c under g^*
  # * Repeat nrep times and average.
#---------------------------------------------------------------------------------
get.MCS_ests.old <- function(data, MC_fit_params, fit_h_reg_obj) {
  family <- MC_fit_params$family

  onlyTMLE_B <- MC_fit_params$onlyTMLE_B
  n <- nrow(data)
  k <- MC_fit_params$Kmax
  n_MCsims <- MC_fit_params$n_MCsims
  lbound <- MC_fit_params$lbound
  max_npwt <- MC_fit_params$max_npwt
  max.err_eps <- MC_fit_params$max.err_eps

  node_l <- MC_fit_params$node_l
  nFnode <- node_l$nFnode
  Anode <- node_l$Anode
  Ynode <- node_l$Ynode
  iidW_flag <- MC_fit_params$iidW_flag
  NetInd_k <- MC_fit_params$NetInd_k

  m.Q.init <- MC_fit_params$m.Q.init
  m.Q.star_h_A <- MC_fit_params$m.Q.star_h_A
  m.Q.star_h_B <- MC_fit_params$m.Q.star_h_B
  m.Q.star_iptw <- MC_fit_params$m.Q.star_iptw
  m.g0N  <- MC_fit_params$m.g0N
  f.gstar <- MC_fit_params$f.gstar
  f.g_args <- MC_fit_params$f.g_args
  f.g0 <- MC_fit_params$f.g0
  f.g0_args <- MC_fit_params$f.g0_args

  # 02/07/12 Eliminated the need to make two passes through the data for monte-carlo # Creates one matrix of all Q(Y|W_i)
  .f_ests_all <- function(NetInd_k, emp_netW) {
    .f.gen.reps <- function(nrep, NetInd_k, emp_netW)   {
      # Generate full sample   c=(A,W) under g* (and iid or not for Ws)
       .f.gen.sample <- function(NetInd_k, n) {
        if (iidW_flag) {  # Version 1: Resample W's with replacement
          resamp_idx <- sample(c(1:n), n, replace=TRUE)
          netW <- NULL # get all network W's from the original data under resampled IDs
          for (Wnode in node_l$Wnodes) {
            netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[resamp_idx,Wnode], Wnode)))
          }
          full_dfW <- netW
        } else {
        # Version 2: No resampling of W's (W's are indep. but not iid, using NPMLE that puts mass 1 on all obs i=1,...,N)
          resamp_idx <- c(1:n)
          full_dfW <- emp_netW
        }
        nFriend <- subset(data, select=nFnode) # nFriends (network) is never resampled
        # print("full_dfW"); print(head(full_dfW,10))
        resamp_A <- f.gen.A.star.old(k, data.frame(full_dfW, subset(data, select=nFnode)), f.gstar, f.g_args)
        full_dfA <- .f.allCovars(k, NetInd_k, resamp_A, Anode)
        determ_df <-  data.frame(determ.g=data$determ.g[resamp_idx], determ.Q=data$determ.Q[resamp_idx]) # get deterministic nodes, also resampled, since fcn of W's
        Y_resamp <- subset(data, select=Ynode) # use observed Y's - INCORRECT, BASED ON NEW W (iidW=TRUE), deterministic Y's might change values
        
        resamp_data <- data.frame(full_dfW, full_dfA, Y_resamp, nFriend, determ_df)
        # print("resamp_data"); print(head(resamp_data))
        return(resamp_data)
      }

      # MLE - Predict E[Y|g_star] (QY.init) for each i, based on the initial model for E[Y|C^Y] (m.Q.init)
      .f_pred_barQ.init <- function(m.Q.init, samp_data) {
        #deterministic nodes for Q
        determ.Q <- samp_data$determ.Q
        # MLE Subs Estimator (est probY based on model for Q_Y)
        # predict only for those not infected at bsl, W2==0
        #************************************************
        # QY <- rep_len(1, nrow(samp_data))
        QY <- predict(m.Q.init, newdata=samp_data, type="response")
        QY[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
        #************************************************
        # QY[!determ.Q] <- predict(m.Q.init, newdata=samp_data[!determ.Q,], type="response")
        return(QY)
      }
      # TMLE - Predict E[Y|g_star] (QY.star) for each i, based on the coefficient epsilon update model for E[Y|C^Y] (m.Q.star_h_A)      
      .f_pred_barQ.star_A <- function(QY.init, samp_data) {
        h_bars <- pred.hbars.old(samp_data, fit_h_reg_obj, NetInd_k)
        h_iptw <- h_bars$df_h_bar_vals$h
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.Q.star_h_A))) {
          off <- qlogis(QY.init)
          QY.star <- plogis(off + coef(m.Q.star_h_A)*h_iptw)
          #************************************************
          # print("determ.Q"); print(sum(determ.Q))
          #************************************************
          # QY.star[determ.Q] <- 1
          QY.star[determ.Q] <- samp_data[determ.Q, Ynode] # will be incorrect when W's are resampled
          #************************************************
          return(QY.star)
        } else {
          return(QY.init)
        }
      }
      # TMLE B - Predict E[Y|g_star] (QY.star) for each i, based on the intercept epsilon update model for E[Y|C^Y] (m.Q.star_h_B)
      .f_pred_barQ.star_B <- function(QY.init, samp_data) {
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.Q.star_h_B))) {
          off <- qlogis(QY.init)
          QY.star_B <- plogis(off + coef(m.Q.star_h_B))
          #************************************************
          # QY.star[determ.Q] <- 1
          QY.star_B[determ.Q] <- samp_data[determ.Q, Ynode]  # will be incorrect when W's are resampled
          # print("QY.star_B"); print(QY.star_B)
          #************************************************
          return(QY.star_B)
        } else {
          return(QY.init)
        }
      }
      # get an estimate of fi_W (hold ALL W's fixed at once) - a component of TMLE Var
      .f.gen.fi_W <- function(NetInd_k, emp_netW) {
        determ_df <-  data.frame(determ.g=data$determ.g, determ.Q=data$determ.Q)
        resamp_A <- f.gen.A.star.old(k, data.frame(emp_netW,subset(data, select=nFnode)), f.gstar, f.g_args)
        samp_dataA <- .f.allCovars(k, NetInd_k, resamp_A, Anode)
        resamp_A_fixW <- data.frame(emp_netW, samp_dataA, subset(data, select=c(nFnode, Ynode)), determ_df)
        # *******fi_W based on Q,N.init model ******
        QY.init_fixW <- .f_pred_barQ.init(m.Q.init, resamp_A_fixW)
        fi_W_init <- QY.init_fixW   # vers 2
        # *******fi_W based on Q,N.star models (A & B) ******
        if (onlyTMLE_B) {
          QY.star_fixW_A <- rep_len(0, nrow(resamp_A_fixW))
          QY.star_fixW_B <- .f_pred_barQ.star_B(QY.init_fixW, resamp_A_fixW)
        } else {
          QY.star_fixW_A <- .f_pred_barQ.star_A(QY.init_fixW, resamp_A_fixW)
          QY.star_fixW_B <- .f_pred_barQ.star_B(QY.init_fixW, resamp_A_fixW)
        }
        return(list(fi_W_init=fi_W_init, fi_W_star_A=QY.star_fixW_A, fi_W_star_B=QY.star_fixW_B))
      }
      # IPTW NETWORK TMLE
      .f.gen_TMLEnetIPTW <- function(QY.init, samp_data, NetInd_k) {
        g_iptw <- iptw_est(k=k, data=samp_data, node_l=node_l, m.gN=m.g0N, f.gstar=f.gstar, f.g_args=f.g_args, 
                            family=family, NetInd_k=NetInd_k, lbound=lbound, max_npwt=max_npwt, f.g0=f.g0, f.g0_args=f.g0_args)
        determ.Q <- samp_data$determ.Q
        if (!is.na(coef(m.Q.star_iptw))) {
          off <- qlogis(QY.init)
          QY.star <- plogis(off + coef(m.Q.star_iptw)*g_iptw)
          #************************************************
          # QY.star[determ.Q] <- 1
          QY.star[determ.Q] <- samp_data[determ.Q, Ynode]
          #************************************************
          return(QY.star)
        }
      }
      #-------------------------------------------
      # Main body of .f.gen.reps()
      #-------------------------------------------
      resamp_d <- .f.gen.sample(NetInd_k, n) # Get a random sample of all A and W
      QY_gstar_mle <- .f_pred_barQ.init(m.Q.init, resamp_d) # QY.init (G-Comp estimator) - est probY based on model for Q_Y
      #-------------------------------------------
      if (onlyTMLE_B) {
        QY_gstar_TMLE_A <- rep_len(0, n) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- rep_len(0, n) # IPTW NETWORK TMLE
      } else {
        QY_gstar_TMLE_A <- .f_pred_barQ.star_A(QY_gstar_mle, resamp_d) # NETWORK TMLE A (adjusted by coefficient epsilon on h_bar ratio)
        QY_gstar_TMLE_B <- .f_pred_barQ.star_B(QY_gstar_mle, resamp_d) # NETWORK TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
        QY_gstar_TMLE_IPTW <- .f.gen_TMLEnetIPTW(QY_gstar_mle, resamp_d, NetInd_k) # IPTW NETWORK TMLE
      }
      fi_Ws_list <- .f.gen.fi_W(NetInd_k, emp_netW) # Get fi_W - hold W fixed to observed values

      # Put all estimators together and add names (defined in G_D_W_1_nms outside of this function):
      mean_psis_all <- c(mean(QY_gstar_mle), mean(QY_gstar_TMLE_A), mean(QY_gstar_TMLE_B), mean(QY_gstar_TMLE_IPTW), 
                        fi_Ws_list$fi_W_init, fi_Ws_list$fi_W_star_A, fi_Ws_list$fi_W_star_B)
      names(mean_psis_all) <- G_D_W_1_nms
      return(mean_psis_all)
    } # end of .f.gen.reps()

    #-------------------------------------------
    # Main body of .f_ests_all()
    #-------------------------------------------
    all_ests_reps <- t(sapply(seq(n_MCsims), .f.gen.reps, NetInd_k, emp_netW))
    return(all_ests_reps)
  }
  #---------------------------------------------------------------------------------
  # Main body of a fcn get.MCS_ests(): MC evalution of the estimators
  #---------------------------------------------------------------------------------
  # Names of all the estimators calculated during MC simulation:
  G_D_W_1_nms <- c("gcomp_mle", "tmle_A","tmle_B", "tmle_iptw",
                  paste("fWi_init_", c(1:n), sep = ""),
                  paste("fWi_star_A_", c(1:n), sep = ""),
                  paste("fWi_star_B_", c(1:n), sep = ""))
  
  #---------------------------------------------------------------------------------
  # Creating matrix of W's (fixed at observed Ws, for evalution of fi_W)
  netW <- NULL
  for (Wnode in node_l$Wnodes) {
    netW <- data.frame(cbind(netW, .f.allCovars(k, NetInd_k, data[,Wnode], Wnode)))
  }
  emp_netW <- netW

  #---------------------------------------------------------------------------------
  # Allow this part to loop, until desired MCS prob_epsilon for all estimators is reached:
  nrepeat <- 1
  psis_reps <- NULL
  G_comp_D_star_W_reps <- NULL
  repeat {
    G_comp_D_star_W_reps <- rbind(G_comp_D_star_W_reps, .f_ests_all(NetInd_k, emp_netW))
    # G_comp_D_star_W_reps <- rbind(G_comp_D_star_W_reps, .f_ests_all(NetInd_k, NetVec_k_D_Wi, NetVec_D_fullsamp))
    psi_est_mean <- apply(G_comp_D_star_W_reps, 2, mean, na.rm = T)
    psi_est_var <- apply(G_comp_D_star_W_reps, 2, var, na.rm = T)
    psi_percerr <- 2 * abs(psi_est_mean * max.err_eps) # estimate the maximum allowed epsilon for each estimator, based pre-defined % error:  
    # prob_epsilon <- psi_est_var / ((n_MCsims*nrepeat) * (max.err_eps)^2)
    prob_percerr <- psi_est_var / ((n_MCsims*nrepeat) * (psi_percerr)^2)
    prob_percerr[psi_est_var < 0.0001] <- 0.0001
    fin_ests_sel <- c(1:3) # final vec of estimators for which error is measured
    if ( (all(prob_percerr[fin_ests_sel] < 0.05)) | (nrepeat >= 100)) {
      break
    }
    nrepeat <- nrepeat + 1
  }   
  # print("nrepeat"); print(nrepeat)
  return(psi_est_mean)
}



# (DEPRECATED)
#---------------------------------------------------------------------------------
# Estimate h_bar under g_0 and g* given observed data and vector of c^Y's
#---------------------------------------------------------------------------------
get_all_ests.old <- function(data, est_obj) {
  # n <- nrow(data)
  # node_l <- data$nodes
  node_l <- est_obj$node_l
  # nFnode <- node_l$nFnode
  # Anode <- node_l$Anode
	Ynode <- node_l$Ynode

	Y <- data[, Ynode]
	determ.Q <- data[, "determ.Q"]
  # determ.g <- data[, "determ.g"]
  QY.init <- data[, "QY.init"] # initial Q fit
  off <- qlogis(QY.init)  # offset

  # new way of getting params (when data is an DatNet.sWsA object)
  # node_l <- data$nodes

  #************************************************
  # ESTIMATORS
  #************************************************

  #************************************************
  # IPTW_h estimator (based on h^*/h_N clever covariate):
  #************************************************
  fit.hbars_t <- system.time(h_bars <- fit.hbars.old(data = data, h_fit_params = est_obj)) # fit the clever covariat
  # fit.hbars_t.new <- system.time(h_bars.new <- fit.hbars.new(data = data, h_fit_params = est_obj)) # fit the clever covariat
  df_h_bar_vals <- h_bars$df_h_bar_vals
  fit_h_reg_obj <- h_bars$fit_h_reg_obj
  h_iptw <- df_h_bar_vals$h

  print("time to fit old h_bars:"); print(fit.hbars_t)
  print("old h est:"); print(head(df_h_bar_vals))

  Y_IPTW_h <- Y
  Y_IPTW_h[!determ.Q] <- Y[!determ.Q] * h_iptw[!determ.Q]
  print("IPW Est (h)"); print(mean(Y_IPTW_h))

  # print("time to fit h_bars new"); print(fit.hbars_t.new)
  # print("h est new"); print(head(h_bars.new$df_h_bar_vals))

  #************************************************
  # TMLE A: estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_A <- glm(Y ~ -1 + h_iptw + offset(off), data = data.frame(Y = data[, Ynode], off = off, h_iptw = h_iptw),
                                						subset = !determ.Q, family = est_obj$family, control = ctrl), GetWarningsToSuppress(TRUE))
	# QY.star <- Y
  QY.star_A <- Y
	if (!is.na(coef(m.Q.star_reg_A))) QY.star_A <- plogis(off + coef(m.Q.star_reg_A) * h_iptw)

  #************************************************
  # TMLE B: estimate the TMLE update via weighted univariate ML (espsilon is intercept)
  #************************************************
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  SuppressGivenWarnings(m.Q.star_reg_B <- glm(Y ~ offset(off), data = data.frame(Y = data[, Ynode], off = off), weights = h_iptw,
                                            subset = !determ.Q, family = est_obj$family, control = ctrl), GetWarningsToSuppress(TRUE))
  QY.star_B <- Y
  if (!is.na(coef(m.Q.star_reg_B))) QY.star_B <- plogis(off + coef(m.Q.star_reg_B))

  #************************************************
  # IPTW estimator (based on full likelihood factorization, prod(g^*)/prod(g_N):
  #************************************************
	# 02/16/13: IPTW estimator (Y_i * prod_{j \in Fi} [g*(A_j|c^A)/g0_N(A_j|c^A)])
	g_iptw <- iptw_est(k = est_obj$Kmax, data = data, node_l = node_l, m.gN = est_obj$m.g0N,
                      f.gstar = est_obj$f.gstar, f.g_args = est_obj$f.g_args, family = est_obj$family,
                      NetInd_k = est_obj$NetInd_k, lbound = est_obj$lbound, max_npwt = est_obj$max_npwt,
                      f.g0 = est_obj$f.g0, f.g0_args = est_obj$f.g0_args)
  Y_IPTW_net <- Y
  Y_IPTW_net[!determ.Q] <- Y[!determ.Q] * g_iptw[!determ.Q]

  #************************************************
  # IPTW-based clever covariate TMLE (based on FULL likelihood factorization), covariate based fluctuation
  #************************************************
	SuppressGivenWarnings(m.Q.star_iptw <- glm(Y ~ -1 + g_iptw + offset(off),
                                						data = data.frame(Y = data[, Ynode], off = off, g_iptw = g_iptw),
                                						subset = !determ.Q, family = est_obj$family, control = ctrl),
                                						GetWarningsToSuppress(TRUE))

  parsubmodel_fits <- rbind(coef(m.Q.star_reg_A), coef(m.Q.star_reg_B), coef(m.Q.star_iptw))
  rownames(parsubmodel_fits) <- c("epsilon (covariate)", "alpha (intercept)", "iptw epsilon (covariate)")
  print("parsubmodel_fits"); print(parsubmodel_fits)


  # run M.C. evaluation estimating psi under g^*:
  MC_fit_params <- append(est_obj,
                      list(m.Q.star_h_A = m.Q.star_reg_A,
                          m.Q.star_h_B = m.Q.star_reg_B,
                          m.Q.star_iptw = m.Q.star_iptw,
                          hstar = df_h_bar_vals$h.star_c,
                          hgN = df_h_bar_vals$h_c,
                          h_tilde = df_h_bar_vals$h))

  syst1 <- system.time(MCS_res <- get.MCS_ests.old(data = data,  MC_fit_params = MC_fit_params, fit_h_reg_obj = fit_h_reg_obj))

  print("time to run old MCS: "); print(syst1);

  psi_tmle_A <- MCS_res[names(MCS_res) == "tmle_A"]
  psi_tmle_B <- MCS_res[names(MCS_res) == "tmle_B"]
  psi_tmle_iptw <- MCS_res[names(MCS_res) == "tmle_iptw"]
  psi_iptw_h <- mean(Y_IPTW_h)  # IPTW estimator based on h - clever covariate
  psi_iptw_g <- mean(Y_IPTW_net)  # IPTW estimator based on full g factorization (prod(g))
  psi_mle <- MCS_res[names(MCS_res) == "gcomp_mle"]

  ests <- as.vector(c(psi_tmle_A, psi_tmle_B, psi_tmle_iptw, psi_iptw_h, psi_iptw_g, psi_mle))
  estnames <- c(      "tmle_A",   "tmle_B",   "tmle_g_iptw", "h_iptw",   "g_iptw", "mle")  
  names(ests) <- estnames
  ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
  ests_mat[, 1] <- ests
  rownames(ests_mat) <- estnames; colnames(ests_mat) <- "estimate"

  wts_mat <- matrix(0L, nrow = length(h_iptw), ncol = 2)
  colnames(wts_mat) <- c("h_wts", "g_wts")
  wts_mat[, "h_wts"] <- h_iptw
  # wts_mat[, "g_wts"] <- g_iptw

  # COMPONENTS OF ASYMPTOTIC VARIANCE FOR TMLE_NET (E_g^*[\bar{Q}_0(A^s,W^s|W^s)]-\psi_0):
  # SUBSTRACT overall estimate of psi_0 from fW_i i-specific components
  fWi_mat <- matrix(0L, nrow = nrow(data), ncol = 3)
  colnames(fWi_mat) <- c("fWi_Qinit", "fWi_Qstar_A", "fWi_Qstar_B")
  fWi_mat[,"fWi_Qinit"] <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_mat[,"fWi_Qstar_A"] <- MCS_res[agrep("fWi_star_A_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  fWi_mat[,"fWi_Qstar_B"] <- MCS_res[agrep("fWi_star_B_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]

  QY_mat <- matrix(0L, nrow = nrow(data), ncol = 3)
  colnames(QY_mat) <- c("QY.init", "QY.star_A", "QY.star_B")
  QY_mat[,"QY.init"] <- QY.init
  QY_mat[,"QY.star_A"] <- QY.star_A
  QY_mat[,"QY.star_B"] <- QY.star_B

  print(c(fWi_init_A = mean(fWi_mat[,"fWi_Qinit"] - ests["tmle_A"]), fWi_init_B = mean(fWi_mat[,"fWi_Qinit"] - ests["tmle_B"])));
  print(c(fWi_star_A = mean(fWi_mat[,"fWi_Qstar_A"] - ests["tmle_A"]), fWi_star_B = mean(fWi_mat[,"fWi_Qstar_B"] - ests["tmle_B"])));

  print("old MC.ests vec: "); print(ests)
  print("old MC.ests mat: "); print(ests_mat)

  return(list( ests_mat = ests_mat,
               wts_mat = wts_mat,
               fWi_mat = fWi_mat,
               QY_mat = QY_mat
              ))
}

# nocov end