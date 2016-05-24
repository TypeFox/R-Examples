# Copyright (C) 2010-2016 Johannes Ranke
# Portions of this code are copyright (C) 2013 Eurofins Regulatory AG
# Contact: jranke@uni-bremen.de
# The summary function is an adapted and extended version of summary.modFit
# from the FME package, v 1.1 by Soetart and Petzoldt, which was in turn
# inspired by summary.nls.lm

# This file is part of the R package mkin

# mkin is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>
if(getRversion() >= '2.15.1') utils::globalVariables(c("name", "time", "value"))

mkinfit <- function(mkinmod, observed,
  parms.ini = "auto",
  state.ini = "auto", 
  fixed_parms = NULL,
  fixed_initials = names(mkinmod$diffs)[-1],
  from_max_mean = FALSE,
  solution_type = c("auto", "analytical", "eigen", "deSolve"),
  method.ode = "lsoda",
  use_compiled = "auto",
  method.modFit = c("Port", "Marq", "SANN", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B"),
  maxit.modFit = "auto",
  control.modFit = list(),
  transform_rates = TRUE,
  transform_fractions = TRUE,
  plot = FALSE, quiet = FALSE,
  err = NULL, weight = "none", scaleVar = FALSE,
  atol = 1e-8, rtol = 1e-10, n.outtimes = 100,
  reweight.method = NULL,
  reweight.tol = 1e-8, reweight.max.iter = 10,
  trace_parms = FALSE,
  ...)
{
  # Check mkinmod and generate a model for the variable whith the highest value
  # if a suitable string is given
  parent_models_available = c("SFO", "FOMC", "DFOP", "HS", "SFORB", "IORE") 
  if (class(mkinmod) != "mkinmod") {
    presumed_parent_name = observed[which.max(observed$value), "name"]
    if (mkinmod[[1]] %in% parent_models_available) {
      speclist <- list(list(type = mkinmod, sink = TRUE))
      names(speclist) <- presumed_parent_name
      mkinmod <- mkinmod(speclist = speclist)
    } else {
      stop("Argument mkinmod must be of class mkinmod or a string containing one of\n  ",
           paste(parent_models_available, collapse = ", "))
    } 
  }

  # Check optimisation method and set maximum number of iterations if specified
  method.modFit = match.arg(method.modFit)
  if (maxit.modFit != "auto") {
    if (method.modFit == "Marq") control.modFit$maxiter = maxit.modFit
    if (method.modFit == "Port") {
      control.modFit$iter.max = maxit.modFit
      control.modFit$eval.max = maxit.modFit
    }
    if (method.modFit %in% c("SANN", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B")) {
        control.modFit$maxit = maxit.modFit
    }
  }

  # Get the names of the state variables in the model
  mod_vars <- names(mkinmod$diffs)

  # Get the names of observed variables
  obs_vars <- names(mkinmod$spec)

  # Subset observed data with names of observed data in the model
  observed <- subset(observed, name %in% obs_vars)

  # Obtain data for decline from maximum mean value if requested
  if (from_max_mean) {
    # This is only used for simple decline models
    if (length(obs_vars) > 1) 
      stop("Decline from maximum is only implemented for models with a single observed variable")

    means <- aggregate(value ~ time, data = observed, mean, na.rm=TRUE)
    t_of_max <- means[which.max(means$value), "time"]
    observed <- subset(observed, time >= t_of_max)
    observed$time <- observed$time - t_of_max
  }

  # Define starting values for parameters where not specified by the user
  if (parms.ini[[1]] == "auto") parms.ini = vector()

  # Prevent inital parameter specifications that are not in the model
  wrongpar.names <- setdiff(names(parms.ini), mkinmod$parms)
  if (length(wrongpar.names) > 0) {
    stop("Initial parameter(s) ", paste(wrongpar.names, collapse = ", "),
         " not used in the model")
  }

  # Warn that the sum of formation fractions may exceed one if they are not
  # fitted in the transformed way
  if (mkinmod$use_of_ff == "max" & transform_fractions == FALSE) {
    warning("The sum of formation fractions may exceed one if you do not use ",
            "transform_fractions = TRUE." )
    for (box in mod_vars) {
      # Stop if formation fractions are not transformed and we have no sink
      if (mkinmod$spec[[box]]$sink == FALSE) {
        stop("If formation fractions are not transformed during the fitting, ",
             "it is not supported to turn off pathways to sink.\n ",
             "Consider turning on the transformation of formation fractions or ",
             "setting up a model with use_of_ff = 'min'.\n")
      }
    }
  }

  # Do not allow fixing formation fractions if we are using the ilr transformation,
  # this is not supported
  if (transform_fractions == TRUE && length(fixed_parms) > 0) {
    if (grepl("^f_", fixed_parms)) {
      stop("Fixing formation fractions is not supported when using the ilr ",
           "transformation.")
    }
  }

  # Set initial parameter values, including a small increment (salt)
  # to avoid linear dependencies (singular matrix) in Eigenvalue based solutions
  k_salt = 0
  defaultpar.names <- setdiff(mkinmod$parms, names(parms.ini))
  for (parmname in defaultpar.names) {
    # Default values for rate constants, depending on the parameterisation
    if (grepl("^k", parmname)) {
      parms.ini[parmname] = 0.1 + k_salt
      k_salt = k_salt + 1e-4
    }
    # Default values for rate constants for reversible binding
    if (grepl("free_bound$", parmname)) parms.ini[parmname] = 0.1 
    if (grepl("bound_free$", parmname)) parms.ini[parmname] = 0.02
    # Default values for IORE exponents
    if (grepl("^N", parmname)) parms.ini[parmname] = 1
    # Default values for the FOMC, DFOP and HS models
    if (parmname == "alpha") parms.ini[parmname] = 1
    if (parmname == "beta") parms.ini[parmname] = 10
    if (parmname == "k1") parms.ini[parmname] = 0.1
    if (parmname == "k2") parms.ini[parmname] = 0.01
    if (parmname == "tb") parms.ini[parmname] = 5
    if (parmname == "g") parms.ini[parmname] = 0.5
  }
  # Default values for formation fractions in case they are present
  for (box in mod_vars) {
    f_names <- mkinmod$parms[grep(paste0("^f_", box), mkinmod$parms)]
    if (length(f_names) > 0) {
      # We need to differentiate between default and specified fractions
      # and set the unspecified to 1 - sum(specified)/n_unspecified
      f_default_names <- intersect(f_names, defaultpar.names)
      f_specified_names <- setdiff(f_names, defaultpar.names)
      sum_f_specified = sum(parms.ini[f_specified_names])
      if (sum_f_specified > 1) {
        stop("Starting values for the formation fractions originating from ",
             box, " sum up to more than 1.")
      }
      if (mkinmod$spec[[box]]$sink) n_unspecified = length(f_default_names) + 1
      else {
        n_unspecified = length(f_default_names)
      }
      parms.ini[f_default_names] <- (1 - sum_f_specified) / n_unspecified
    }
  }

  # Set default for state.ini if appropriate
  parent_name = names(mkinmod$spec)[[1]]
  if (state.ini[1] == "auto") {
    parent_time_0 = subset(observed, time == 0 & name == parent_name)$value
    parent_time_0_mean = mean(parent_time_0, na.rm = TRUE)
    if (is.na(parent_time_0_mean)) {
      state.ini = c(100, rep(0, length(mkinmod$diffs) - 1))
    } else {
      state.ini = c(parent_time_0_mean, rep(0, length(mkinmod$diffs) - 1))
    }
  }

  # Name the inital state variable values if they are not named yet
  if(is.null(names(state.ini))) names(state.ini) <- mod_vars

  # Transform initial parameter values for fitting
  transparms.ini <- transform_odeparms(parms.ini, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Parameters to be optimised:
  # Kinetic parameters in parms.ini whose names are not in fixed_parms
  parms.fixed <- parms.ini[fixed_parms]
  parms.optim <- parms.ini[setdiff(names(parms.ini), fixed_parms)]

  transparms.fixed <- transform_odeparms(parms.fixed, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)
  transparms.optim <- transform_odeparms(parms.optim, mkinmod,
                                       transform_rates = transform_rates,
                                       transform_fractions = transform_fractions)

  # Inital state variables in state.ini whose names are not in fixed_initials
  state.ini.fixed <- state.ini[fixed_initials]
  state.ini.optim <- state.ini[setdiff(names(state.ini), fixed_initials)]

  # Preserve names of state variables before renaming initial state variable
  # parameters
  state.ini.optim.boxnames <- names(state.ini.optim)
  state.ini.fixed.boxnames <- names(state.ini.fixed)
  if(length(state.ini.optim) > 0) {
    names(state.ini.optim) <- paste(names(state.ini.optim), "0", sep="_")
  }
  if(length(state.ini.fixed) > 0) {
    names(state.ini.fixed) <- paste(names(state.ini.fixed), "0", sep="_")
  }

  # Decide if the solution of the model can be based on a simple analytical
  # formula, the spectral decomposition of the matrix (fundamental system)
  # or a numeric ode solver from the deSolve package
  # Prefer deSolve over eigen if a compiled model is present and use_compiled
  # is not set to FALSE
  solution_type = match.arg(solution_type)
  if (solution_type == "analytical" && length(mkinmod$spec) > 1)
     stop("Analytical solution not implemented for models with metabolites.")
  if (solution_type == "eigen" && !is.matrix(mkinmod$coefmat))
     stop("Eigenvalue based solution not possible, coefficient matrix not present.")
  if (solution_type == "auto") {
    if (length(mkinmod$spec) == 1) {
      solution_type = "analytical"
    } else {
      if (!is.null(mkinmod$cf) & use_compiled[1] != FALSE) {
        solution_type = "deSolve"
      } else {
        if (is.matrix(mkinmod$coefmat)) {
          solution_type = "eigen"
          if (max(observed$value, na.rm = TRUE) < 0.1) {
            stop("The combination of small observed values (all < 0.1) and solution_type = eigen is error-prone")
          }
        } else {
          solution_type = "deSolve"
        }
      }
    }
  }

  # Define outtimes for model solution.
  # Include time points at which observed data are available
  # Make sure we include time 0, so initial values for state variables are for time 0
  outtimes = sort(unique(c(observed$time, seq(min(observed$time),
                                              max(observed$time),
                                              length.out = n.outtimes))))


  cost.old <- 1e100 # The first model cost should be smaller than this value
  calls <- 0 # Counter for number of model solutions
  out_predicted <- NA

  # Define the model cost function for optimisation, including (back)transformations
  cost <- function(P)
  {
    assign("calls", calls+1, inherits=TRUE) # Increase the model solution counter

    # Trace parameter values if requested
    if(trace_parms) cat(P, "\n")

    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
      odeini <- state.ini.fixed
      names(odeini) <- state.ini.fixed.boxnames
    }

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], transparms.fixed)

    parms <- backtransform_odeparms(odeparms, mkinmod,
                                    transform_rates = transform_rates,
                                    transform_fractions = transform_fractions)

    # Solve the system with current transformed parameter values
    out <- mkinpredict(mkinmod, parms, 
                       odeini, outtimes, 
                       solution_type = solution_type, 
                       use_compiled = use_compiled,
                       method.ode = method.ode,
                       atol = atol, rtol = rtol, ...)

    assign("out_predicted", out, inherits=TRUE)

    mC <- modCost(out, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    # Report and/or plot if the model is improved
    if (mC$model < cost.old) {
      if(!quiet) cat("Model cost at call ", calls, ": ", mC$model, "\n")

      # Plot the data and current model output if requested
      if(plot) {
        outtimes_plot = seq(min(observed$time), max(observed$time), length.out=100)

        out_plot <- mkinpredict(mkinmod, parms,
                                odeini, outtimes_plot, 
                                solution_type = solution_type, 
                                use_compiled = use_compiled,
                                method.ode = method.ode,
                                atol = atol, rtol = rtol, ...)

        plot(0, type="n", 
          xlim = range(observed$time), ylim = c(0, max(observed$value, na.rm=TRUE)),
          xlab = "Time", ylab = "Observed")
        col_obs <- pch_obs <- 1:length(obs_vars)
        lty_obs <- rep(1, length(obs_vars))
        names(col_obs) <- names(pch_obs) <- names(lty_obs) <- obs_vars
        for (obs_var in obs_vars) {
          points(subset(observed, name == obs_var, c(time, value)), 
                 pch = pch_obs[obs_var], col = col_obs[obs_var])
        }
        matlines(out_plot$time, out_plot[-1], col = col_obs, lty = lty_obs)
        legend("topright", inset=c(0.05, 0.05), legend=obs_vars, 
          col=col_obs, pch=pch_obs, lty=1:length(pch_obs))
      }
    
      assign("cost.old", mC$model, inherits=TRUE)
    }
    return(mC)
  }

  # Define the model cost function for the t-test, without parameter transformation
  cost_notrans <- function(P)
  {
    if(length(state.ini.optim) > 0) {
      odeini <- c(P[1:length(state.ini.optim)], state.ini.fixed)
      names(odeini) <- c(state.ini.optim.boxnames, state.ini.fixed.boxnames)
    } else {
      odeini <- state.ini.fixed
      names(odeini) <- state.ini.fixed.boxnames
    }

    odeparms <- c(P[(length(state.ini.optim) + 1):length(P)], parms.fixed)

    # Solve the system with current transformed parameter values
    out <- mkinpredict(mkinmod, odeparms, 
                       odeini, outtimes, 
                       solution_type = solution_type, 
                       use_compiled = use_compiled,
                       method.ode = method.ode,
                       atol = atol, rtol = rtol, ...)

    mC <- modCost(out, observed, y = "value",
      err = err, weight = weight, scaleVar = scaleVar)

    return(mC)
  }

  # Define lower and upper bounds other than -Inf and Inf for parameters
  # for which no internal transformation is requested in the call to mkinfit.
  lower <- rep(-Inf, length(c(state.ini.optim, transparms.optim)))
  upper <- rep(Inf, length(c(state.ini.optim, transparms.optim)))
  names(lower) <- names(upper) <- c(names(state.ini.optim), names(transparms.optim))
  if (!transform_rates) {
    index_k <- grep("^k_", names(lower))
    lower[index_k] <- 0
    index_k__iore <- grep("^k__iore_", names(lower))
    lower[index_k__iore] <- 0
    other_rate_parms <- intersect(c("alpha", "beta", "k1", "k2", "tb"), names(lower))
    lower[other_rate_parms] <- 0
  }

  if (!transform_fractions) {
    index_f <- grep("^f_", names(upper))
    lower[index_f] <- 0
    upper[index_f] <- 1
    other_fraction_parms <- intersect(c("g"), names(upper))
    lower[other_fraction_parms] <- 0
    upper[other_fraction_parms] <- 1
  }

  # Do the fit and take the time
  fit_time <- system.time({
    fit <- modFit(cost, c(state.ini.optim, transparms.optim), 
                  method = method.modFit, control = control.modFit, 
                  lower = lower, upper = upper, ...)

    # Reiterate the fit until convergence of the variance components (IRLS)
    # if requested by the user
    weight.ini = weight
    if (!is.null(err)) weight.ini = "manual"

    if (!is.null(reweight.method)) {
      if (reweight.method != "obs") stop("Only reweighting method 'obs' is implemented")
      if(!quiet) {
        cat("IRLS based on variance estimates for each observed variable\n")
      }
      if (!quiet) {
        cat("Initial variance estimates are:\n")
        print(signif(fit$var_ms_unweighted, 8))
      }
      reweight.diff = 1
      n.iter <- 0
      if (!is.null(err)) observed$err.ini <- observed[[err]]
      err = "err.irls"
      while (reweight.diff > reweight.tol & n.iter < reweight.max.iter) {
        n.iter <- n.iter + 1
        sigma.old <- sqrt(fit$var_ms_unweighted)
        observed[err] <- sqrt(fit$var_ms_unweighted)[as.character(observed$name)]
        fit <- modFit(cost, fit$par, method = method.modFit,
                      control = control.modFit, lower = lower, upper = upper, ...)
        reweight.diff = sum((sqrt(fit$var_ms_unweighted) - sigma.old)^2)
        if (!quiet) {
          cat("Iteration", n.iter, "yields variance estimates:\n")
          print(signif(fit$var_ms_unweighted, 8))
          cat("Sum of squared differences to last variance estimates:",
              signif(reweight.diff, 2), "\n")
        }
      }
    }
  })

  # Check for convergence
  if (method.modFit == "Marq") {
    if (!fit$info %in% c(1, 2, 3)) {
      fit$warning = paste0("Optimisation by method ", method.modFit, 
                           " did not converge.\n",
                           "The message returned by nls.lm is:\n",
                                    fit$message)
      warning(fit$warning)
    }
    else {
      if(!quiet) cat("Optimisation by method", method.modFit, "successfully terminated.\n")
    }
  }
  if (method.modFit %in% c("Port", "SANN", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B")) {
    if (fit$convergence != 0) {
      fit$warning = paste0("Optimisation by method ", method.modFit, 
                           " did not converge.\n",
                           "Convergence code is ", fit$convergence,
                           ifelse(is.null(fit$message), "", 
                                  paste0("\nMessage is ", fit$message)))
      warning(fit$warning)
    } else {
      if(!quiet) cat("Optimisation by method", method.modFit, "successfully terminated.\n")
    }
  }

  # Return number of iterations for SANN method
  if (method.modFit == "SANN") {
    fit$iter = maxit.modFit
    if(!quiet) cat("Termination of the SANN algorithm does not imply convergence.\n")
  }

  # We need to return some more data for summary and plotting
  fit$solution_type <- solution_type
  fit$transform_rates <- transform_rates
  fit$transform_fractions <- transform_fractions
  fit$method.modFit <- method.modFit
  fit$maxit.modFit <- maxit.modFit
  fit$calls <- calls
  fit$time <- fit_time

  # We also need the model for summary and plotting
  fit$mkinmod <- mkinmod

  # We need data and predictions for summary and plotting
  fit$observed <- observed
  fit$obs_vars <- obs_vars
  fit$predicted <- mkin_wide_to_long(out_predicted, time = "time")

  # Backtransform parameters
  bparms.optim = backtransform_odeparms(fit$par, fit$mkinmod,
                                        transform_rates = transform_rates,
                                        transform_fractions = transform_fractions)
  bparms.fixed = c(state.ini.fixed, parms.fixed)
  bparms.all = c(bparms.optim, parms.fixed)

  # Attach the cost functions to the fit for post-hoc parameter uncertainty analysis
  fit$cost <- cost
  fit$cost_notrans <- cost_notrans

  # Estimate the Hessian for the model cost without parameter transformations
  # to make it possible to obtain the usual t-test
  # Code ported from FME::modFit
  Jac_notrans <- gradient(function(p, ...) cost_notrans(p)$residuals$res, 
                          bparms.optim, centered = TRUE)
  fit$hessian_notrans <- 2 * t(Jac_notrans) %*% Jac_notrans

  # Collect initial parameter values in three dataframes
  fit$start <- data.frame(value = c(state.ini.optim, 
                                    parms.optim))
  fit$start$type = c(rep("state", length(state.ini.optim)), 
                     rep("deparm", length(parms.optim)))

  fit$start_transformed = data.frame(
      value = c(state.ini.optim, transparms.optim),
      lower = lower,
      upper = upper)

  fit$fixed <- data.frame(value = c(state.ini.fixed, parms.fixed))
  fit$fixed$type = c(rep("state", length(state.ini.fixed)), 
                     rep("deparm", length(parms.fixed)))

  # Collect observed, predicted, residuals and weighting
  data <- merge(fit$observed, fit$predicted, by = c("time", "name"))
  data$name <- ordered(data$name, levels = obs_vars)
  data <- data[order(data$name, data$time), ]

  fit$data <- data.frame(time = data$time,
                         variable = data$name,
                         observed = data$value.x,
                         predicted = data$value.y)
  fit$data$residual <- fit$data$observed - fit$data$predicted
  if (!is.null(data$err.ini)) fit$data$err.ini <- data$err.ini
  if (!is.null(err)) fit$data[[err]] <- data[[err]]

  fit$atol <- atol
  fit$rtol <- rtol
  fit$weight.ini <- weight.ini
  fit$reweight.method <- reweight.method
  fit$reweight.tol <- reweight.tol
  fit$reweight.max.iter <- reweight.max.iter

  # Return different sets of backtransformed parameters for summary and plotting
  fit$bparms.optim <- bparms.optim 
  fit$bparms.fixed <- bparms.fixed

  # Return ode and state parameters for further fitting
  fit$bparms.ode <- bparms.all[mkinmod$parms] 
  fit$bparms.state <- c(bparms.all[setdiff(names(bparms.all), names(fit$bparms.ode))],
                        state.ini.fixed)
  names(fit$bparms.state) <- gsub("_0$", "", names(fit$bparms.state))

  fit$date <- date()

  class(fit) <- c("mkinfit", "modFit") 
  return(fit)
}

summary.mkinfit <- function(object, data = TRUE, distimes = TRUE, alpha = 0.05, ...) {
  param  <- object$par
  pnames <- names(param)
  bpnames <- names(object$bparms.optim)
  p      <- length(param)
  mod_vars <- names(object$mkinmod$diffs)
  covar  <- try(solve(0.5*object$hessian), silent = TRUE)   # unscaled covariance
  covar_notrans  <- try(solve(0.5*object$hessian_notrans), silent = TRUE)   # unscaled covariance
  rdf    <- object$df.residual
  resvar <- object$ssr / rdf
  if (!is.numeric(covar)) {
    covar <- NULL
    se <- lci <- uci <- rep(NA, p)
  } else {
    rownames(covar) <- colnames(covar) <- pnames
    se     <- sqrt(diag(covar) * resvar)
    lci    <- param + qt(alpha/2, rdf) * se
    uci    <- param + qt(1-alpha/2, rdf) * se
  }

  if (!is.numeric(covar_notrans)) {
    covar_notrans <- NULL
    se_notrans <- tval <- pval <- rep(NA, p)
  } else {
    rownames(covar_notrans) <- colnames(covar_notrans) <- bpnames
    se_notrans     <- sqrt(diag(covar_notrans) * resvar)
    tval   <- object$bparms.optim/se_notrans
    pval   <- pt(abs(tval), rdf, lower.tail = FALSE)
  }

  names(se) <- pnames
  modVariance <- object$ssr / length(object$residuals)

  param <- cbind(param, se, lci, uci)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "Lower", "Upper"))

  bparam <- cbind(Estimate = object$bparms.optim, se_notrans, "t value" = tval, "Pr(>t)" = pval, Lower = NA, Upper = NA)

  # Transform boundaries of CI for one parameter at a time,
  # with the exception of sets of formation fractions (single fractions are OK).
  f_names_skip <- character(0)
  for (box in mod_vars) { # Figure out sets of fractions to skip
    f_names <- grep(paste("^f", box, sep = "_"), pnames, value = TRUE)
    n_paths <- length(f_names)
    if (n_paths > 1) f_names_skip <- c(f_names_skip, f_names)
  }

  for (pname in pnames) {
    if (!pname %in% f_names_skip) {
      par.lower <- param[pname, "Lower"]
      par.upper <- param[pname, "Upper"]
      names(par.lower) <- names(par.upper) <- pname
      bpl <- backtransform_odeparms(par.lower, object$mkinmod, 
                                            object$transform_rates, 
                                            object$transform_fractions)
      bpu <- backtransform_odeparms(par.upper, object$mkinmod,
                                            object$transform_rates, 
                                            object$transform_fractions)
      bparam[names(bpl), "Lower"] <- bpl
      bparam[names(bpu), "Upper"] <- bpu
    } 
  }

  ans <- list(
    version = as.character(utils::packageVersion("mkin")),
    Rversion = paste(R.version$major, R.version$minor, sep="."),
	  date.fit = object$date,
	  date.summary = date(),
	  solution_type = object$solution_type,
	  method.modFit = object$method.modFit,
	  warning = object$warning,
	  use_of_ff = object$mkinmod$use_of_ff,
    weight.ini = object$weight.ini,
    reweight.method = object$reweight.method,
    residuals = object$residuals,
    residualVariance = resvar,
    sigma = sqrt(resvar),
    modVariance = modVariance,
    df = c(p, rdf), 
    cov.unscaled = covar,
    cov.scaled = covar * resvar,
    info = object$info, 
    niter = object$iterations,
    calls = object$calls,
    time = object$time,
    stopmess = message,
    par = param,
    bpar = bparam)

  ans$diffs <- object$mkinmod$diffs
  if(data) ans$data <- object$data
  ans$start <- object$start
  ans$start_transformed <- object$start_transformed

  ans$fixed <- object$fixed

  ans$errmin <- mkinerrmin(object, alpha = 0.05)

  ans$bparms.ode <- object$bparms.ode
  ep <- endpoints(object)
  if (length(ep$ff) != 0)
    ans$ff <- ep$ff
  if(distimes) ans$distimes <- ep$distimes
  if(length(ep$SFORB) != 0) ans$SFORB <- ep$SFORB
  class(ans) <- c("summary.mkinfit", "summary.modFit") 
  return(ans)  
}

# Expanded from print.summary.modFit
print.summary.mkinfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("mkin version:   ", x$version, "\n")
  cat("R version:      ", x$Rversion, "\n")
  cat("Date of fit:    ", x$date.fit, "\n")
  cat("Date of summary:", x$date.summary, "\n")

  if (!is.null(x$warning)) cat("\n\nWarning:", x$warning, "\n\n")

  cat("\nEquations:\n")
  writeLines(strwrap(x[["diffs"]], exdent = 11))
  df  <- x$df
  rdf <- df[2]

  cat("\nModel predictions using solution type", x$solution_type, "\n")

  cat("\nFitted with method", x$method.modFit, 
      "using", x$calls, "model solutions performed in", x$time[["elapsed"]],  "s\n")

  cat("\nWeighting:", x$weight.ini)
  if(!is.null(x$reweight.method)) cat(" then iterative reweighting method",
                                      x$reweight.method)
  cat("\n")

  cat("\nStarting values for parameters to be optimised:\n")
  print(x$start)

  cat("\nStarting values for the transformed parameters actually optimised:\n")
  print(x$start_transformed)

  cat("\nFixed parameter values:\n")
  if(length(x$fixed$value) == 0) cat("None\n")
  else print(x$fixed)
  
  cat("\nOptimised, transformed parameters with symmetric confidence intervals:\n")
  print(signif(x$par, digits = digits))

  if (x$calls > 0) {
    cat("\nParameter correlation:\n")
    if (!is.null(x$cov.unscaled)){
      Corr <- cov2cor(x$cov.unscaled)
      rownames(Corr) <- colnames(Corr) <- rownames(x$par)
      print(Corr, digits = digits, ...)
    } else {
      cat("Could not estimate covariance matrix; singular system:\n")
    }
  }

  cat("\nResidual standard error:",
      format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

  cat("\nBacktransformed parameters:\n")
  cat("Confidence intervals for internally transformed parameters are asymmetric.\n")
  if ((x$version) < "0.9-36") {
    cat("To get the usual (questionable) t-test, upgrade mkin and repeat the fit.\n")
    print(signif(x$bpar, digits = digits))
  } else {
    cat("t-test (unrealistically) based on the assumption of normal distribution\n")
    cat("for estimators of untransformed parameters.\n")
    print(signif(x$bpar[, c(1, 3, 4, 5, 6)], digits = digits))
  }

  cat("\nChi2 error levels in percent:\n")
  x$errmin$err.min <- 100 * x$errmin$err.min
  print(x$errmin, digits=digits,...)

  printSFORB <- !is.null(x$SFORB)
  if(printSFORB){
    cat("\nEstimated Eigenvalues of SFORB model(s):\n")
    print(x$SFORB, digits=digits,...)
  }

  printff <- !is.null(x$ff)
  if(printff){
    cat("\nResulting formation fractions:\n")
    print(data.frame(ff = x$ff), digits=digits,...)
  }

  printdistimes <- !is.null(x$distimes)
  if(printdistimes){
    cat("\nEstimated disappearance times:\n")
    print(x$distimes, digits=digits,...)
  }

  printdata <- !is.null(x$data)
  if (printdata){
    cat("\nData:\n")
    print(format(x$data, digits = digits, ...), row.names = FALSE)
  }

  invisible(x)
}
# vim: set ts=2 sw=2 expandtab:
