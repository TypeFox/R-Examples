endpoints <- function(fit) {
  # Calculate dissipation times DT50 and DT90 and formation
  # fractions as well as SFORB eigenvalues from optimised parameters
  # Additional DT50 values are calculated from the FOMC DT90 and k1 and k2 from
  # HS and DFOP, as well as from Eigenvalues b1 and b2 of any SFORB models
  ep <- list()
  obs_vars <- fit$obs_vars
  parms.all <- c(fit$bparms.optim, fit$bparms.fixed)
  ep$ff <- vector()
  ep$SFORB <- vector()
  ep$distimes <- data.frame(DT50 = rep(NA, length(obs_vars)), 
			    DT90 = rep(NA, length(obs_vars)), 
    row.names = obs_vars)
  for (obs_var in obs_vars) {
    type = names(fit$mkinmod$map[[obs_var]])[1]  

    # Get formation fractions if directly fitted, and calculate remaining fraction to sink
    f_names = grep(paste("f", obs_var, sep = "_"), names(parms.all), value=TRUE)
    if (length(f_names) > 0) {
      f_values = parms.all[f_names]
      f_to_sink = 1 - sum(f_values)
      names(f_to_sink) = ifelse(type == "SFORB", 
                              paste(obs_var, "free", "sink", sep = "_"), 
                              paste(obs_var, "sink", sep = "_"))
      for (f_name in f_names) {
        ep$ff[[sub("f_", "", sub("_to_", "_", f_name))]] = f_values[[f_name]]
      }
      ep$ff = append(ep$ff, f_to_sink)
    }

    # Get the rest
    if (type == "SFO") {
      k_names = grep(paste("k", obs_var, sep="_"), names(parms.all), value=TRUE)
      k_tot = sum(parms.all[k_names])
      DT50 = log(2)/k_tot
      DT90 = log(10)/k_tot
      if (fit$mkinmod$use_of_ff == "min") {
        for (k_name in k_names)
        {
          ep$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
        }
      }
    }
    if (type == "FOMC") {
      alpha = parms.all["alpha"]
      beta = parms.all["beta"]
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      ep$distimes[obs_var, c("DT50back")] = DT50_back
    }
    if (type == "IORE") {
      k_names = grep(paste("k__iore", obs_var, sep="_"), names(parms.all), value=TRUE)
      k_tot = sum(parms.all[k_names])
      # From the NAFTA kinetics guidance, p. 5
      n = parms.all[paste("N", obs_var, sep = "_")]
      k = k_tot
      # Use the initial concentration of the parent compound
      source_name = fit$mkinmod$map[[1]][[1]]
      c0 = parms.all[paste(source_name, "0", sep = "_")]
      alpha = 1 / (n - 1)
      beta = (c0^(1 - n))/(k * (n - 1))
      DT50 = beta * (2^(1/alpha) - 1)
      DT90 = beta * (10^(1/alpha) - 1)
      DT50_back = DT90 / (log(10)/log(2)) # Backcalculated DT50 as recommended in FOCUS 2011
      ep$distimes[obs_var, c("DT50back")] = DT50_back
      if (fit$mkinmod$use_of_ff == "min") {
        for (k_name in k_names)
        {
          ep$ff[[sub("k_", "", k_name)]] = parms.all[[k_name]] / k_tot
        }
      }
    }
    if (type == "DFOP") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      g = parms.all["g"]
      f <- function(log_t, x) {
        t <- exp(log_t)
        fraction <- g * exp( - k1 * t) + (1 - g) * exp( - k2 * t)
        (fraction - (1 - x/100))^2
      }
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      DT90_k1 = log(10)/k1
      DT90_k2 = log(10)/k2

      DT50 <- try(exp(optimize(f, c(log(DT50_k1), log(DT50_k2)), x=50)$minimum),
                  silent = TRUE)
      DT90 <- try(exp(optimize(f, c(log(DT90_k1), log(DT90_k2)), x=90)$minimum),
                  silent = TRUE)
      if (inherits(DT50, "try-error")) DT50 = NA
      if (inherits(DT90, "try-error")) DT90 = NA
      
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "HS") {
      k1 = parms.all["k1"]
      k2 = parms.all["k2"]
      tb = parms.all["tb"]
      DTx <- function(x) {
        DTx.a <- (log(100/(100 - x)))/k1
        DTx.b <- tb + (log(100/(100 - x)) - k1 * tb)/k2
        if (DTx.a < tb) DTx <- DTx.a
        else DTx <- DTx.b
        return(DTx)
      }
      DT50 <- DTx(50)
      DT90 <- DTx(90)
      DT50_k1 = log(2)/k1
      DT50_k2 = log(2)/k2
      ep$distimes[obs_var, c("DT50_k1")] = DT50_k1
      ep$distimes[obs_var, c("DT50_k2")] = DT50_k2
    }
    if (type == "SFORB") {
      # FOCUS kinetics (2006), p. 60 f
      k_out_names = grep(paste("k", obs_var, "free", sep="_"), names(parms.all), value=TRUE)
      k_out_names = setdiff(k_out_names, paste("k", obs_var, "free", "bound", sep="_"))
      k_1output = sum(parms.all[k_out_names])
      k_12 = parms.all[paste("k", obs_var, "free", "bound", sep="_")]
      k_21 = parms.all[paste("k", obs_var, "bound", "free", sep="_")]

      sqrt_exp = sqrt(1/4 * (k_12 + k_21 + k_1output)^2 + k_12 * k_21 - (k_12 + k_1output) * k_21)
      b1 = 0.5 * (k_12 + k_21 + k_1output) + sqrt_exp
      b2 = 0.5 * (k_12 + k_21 + k_1output) - sqrt_exp

      DT50_b1 = log(2)/b1
      DT50_b2 = log(2)/b2
      DT90_b1 = log(10)/b1
      DT90_b2 = log(10)/b2

      SFORB_fraction = function(t) {
        ((k_12 + k_21 - b1)/(b2 - b1)) * exp(-b1 * t) +
        ((k_12 + k_21 - b2)/(b1 - b2)) * exp(-b2 * t)
      }

      f_50 <- function(log_t) (SFORB_fraction(exp(log_t)) - 0.5)^2
      log_DT50 <- try(optimize(f_50, c(log(DT50_b1), log(DT50_b2)))$minimum,
                      silent = TRUE)
      f_90 <- function(log_t) (SFORB_fraction(exp(log_t)) - 0.1)^2
      log_DT90 <- try(optimize(f_90, c(log(DT90_b1), log(DT90_b2)))$minimum,
                      silent = TRUE)

      DT50 = if (inherits(log_DT50, "try-error")) NA
             else exp(log_DT50)
      DT90 = if (inherits(log_DT90, "try-error")) NA
             else exp(log_DT90)

      for (k_out_name in k_out_names)
      {
        ep$ff[[sub("k_", "", k_out_name)]] = parms.all[[k_out_name]] / k_1output
      }

      # Return the eigenvalues for comparison with DFOP rate constants
      ep$SFORB[[paste(obs_var, "b1", sep="_")]] = b1
      ep$SFORB[[paste(obs_var, "b2", sep="_")]] = b2

      ep$distimes[obs_var, c(paste("DT50", obs_var, "b1", sep = "_"))] = DT50_b1
      ep$distimes[obs_var, c(paste("DT50", obs_var, "b2", sep = "_"))] = DT50_b2
    }
    ep$distimes[obs_var, c("DT50", "DT90")] = c(DT50, DT90)
  }
  return(ep)
}
