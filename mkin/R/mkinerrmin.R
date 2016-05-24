# Copyright (C) 2010-2014 Johannes Ranke
# Contact: jranke@uni-bremen.de

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
if(getRversion() >= '2.15.1') utils::globalVariables(c("name", "value_mean"))

mkinerrmin <- function(fit, alpha = 0.05)
{
  parms.optim <- fit$par

  kinerrmin <- function(errdata, n.parms) {
    means.mean <- mean(errdata$value_mean, na.rm = TRUE)
    df = length(errdata$value_mean) - n.parms
  
    err.min <- sqrt((1 / qchisq(1 - alpha, df)) *
               sum((errdata$value_mean - errdata$value_pred)^2)/(means.mean^2))

    return(list(err.min = err.min, n.optim = n.parms, df = df))
  }

  means <- aggregate(value ~ time + name, data = fit$observed, mean, na.rm=TRUE)
  errdata <- merge(means, fit$predicted, by = c("time", "name"), 
    suffixes = c("_mean", "_pred"))
  errdata <- errdata[order(errdata$time, errdata$name), ]

  # Remove values at time zero for variables whose value for state.ini is fixed,
  # as these will not have any effect in the optimization and should therefore not 
  # be counted as degrees of freedom.
  fixed_initials = gsub("_0$", "", rownames(subset(fit$fixed, type = "state")))
  errdata <- subset(errdata, !(time == 0 & name %in% fixed_initials))

  n.optim.overall <- length(parms.optim)

  errmin.overall <- kinerrmin(errdata, n.optim.overall)
  errmin <- data.frame(err.min = errmin.overall$err.min,
    n.optim = errmin.overall$n.optim, df = errmin.overall$df)
  rownames(errmin) <- "All data"

  # The degrees of freedom are counted according to FOCUS kinetics (2011, p. 164)
  for (obs_var in fit$obs_vars)
  {
    errdata.var <- subset(errdata, name == obs_var)

    # Check if initial value is optimised
    n.initials.optim <- length(grep(paste(obs_var, ".*", "_0", sep=""), names(parms.optim)))

    # Rate constants and IORE exponents are attributed to the source variable
    n.k.optim <- length(grep(paste("^k", obs_var, sep="_"), names(parms.optim)))
    n.k.optim <- n.k.optim + length(grep(paste("^log_k", obs_var, sep="_"), 
                                         names(parms.optim)))
    n.k__iore.optim <- length(grep(paste("^k__iore", obs_var, sep="_"), names(parms.optim)))
    n.k__iore.optim <- n.k__iore.optim + length(grep(paste("^log_k__iore", obs_var, 
							 sep = "_"),
						   names(parms.optim)))

    n.N.optim <- length(grep(paste("^N", obs_var, sep="_"), names(parms.optim)))

    n.ff.optim <- 0
    # Formation fractions are attributed to the target variable, so look
    # for source compartments with formation fractions
    for (source_var in fit$obs_vars) {
      n.ff.source = length(grep(paste("^f", source_var, sep = "_"),
                                 names(parms.optim)))
      n.paths.source = length(fit$mkinmod$spec[[source_var]]$to)
      for (target_var in fit$mkinmod$spec[[source_var]]$to) {
        if (obs_var == target_var) {
          n.ff.optim <- n.ff.optim + n.ff.source/n.paths.source
        }
      }
    }

    n.optim <- sum(n.initials.optim, n.k.optim, n.k__iore.optim, n.N.optim, n.ff.optim)

    # FOMC, DFOP and HS parameters are only counted if we are looking at the
    # first variable in the model which is always the source variable
    if (obs_var == fit$obs_vars[[1]]) {
      special_parms = c("alpha", "log_alpha", "beta", "log_beta", 
                        "k1", "log_k1", "k2", "log_k2", 
                        "g", "g_ilr", "tb", "log_tb")
      n.optim <- n.optim + length(intersect(special_parms, names(parms.optim)))
    }

    # Calculate and add a line to the dataframe holding the results
    errmin.tmp <- kinerrmin(errdata.var, n.optim)
    errmin[obs_var, c("err.min", "n.optim", "df")] <- errmin.tmp
  }

  return(errmin)
}
