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

transform_odeparms <- function(parms, mkinmod,
                               transform_rates = TRUE, 
                               transform_fractions = TRUE) 
{
  # We need the model specification for the names of the model 
  # variables and the information on the sink
  spec = mkinmod$spec

  # Set up container for transformed parameters
  transparms <- numeric(0)

  # Do not transform initial values for state variables
  state.ini.optim <- parms[grep("_0$", names(parms))]
  transparms[names(state.ini.optim)] <- state.ini.optim

  # Log transformation for rate constants if requested
  k <- parms[grep("^k_", names(parms))]
  k__iore <- parms[grep("^k__iore_", names(parms))]
  k <- c(k, k__iore)
  if (length(k) > 0) {
    if(transform_rates) {
      transparms[paste0("log_", names(k))] <- log(k)
    } else transparms[names(k)] <- k
  }

  # Do not transform exponents in IORE models
  N <- parms[grep("^N", names(parms))]
  transparms[names(N)] <- N

  # Go through state variables and apply isometric logratio transformation to
  # formation fractions if requested
  mod_vars = names(spec)
  for (box in mod_vars) {
    f <- parms[grep(paste("^f", box, sep = "_"), names(parms))]

    if (length(f) > 0) {
      if(transform_fractions) {
        if (spec[[box]]$sink) {
          trans_f <- ilr(c(f, 1 - sum(f)))
          trans_f_names <- paste("f", box, "ilr", 1:length(trans_f), sep = "_")
          transparms[trans_f_names] <- trans_f
        } else {
          if (length(f) > 1) {
            trans_f <- ilr(f)
            trans_f_names <- paste("f", box, "ilr", 1:length(trans_f), sep = "_")
            transparms[trans_f_names] <- trans_f
          }
        }
      } else {
        transparms[names(f)] <- f
      }
    }
  }

  # Transform also FOMC parameters alpha and beta, DFOP and HS rates k1 and k2
  # and HS parameter tb if transformation of rates is requested
  for (pname in c("alpha", "beta", "k1", "k2", "tb")) {
    if (!is.na(parms[pname])) {
      if (transform_rates) {
        transparms[paste0("log_", pname)] <- log(parms[pname])
      } else {
        transparms[pname] <- parms[pname]
      }
    } 
  }

  # DFOP parameter g is treated as a fraction
  if (!is.na(parms["g"])) {
    g <- parms["g"]
    if (transform_fractions) {
      transparms["g_ilr"] <- ilr(c(g, 1 - g))
    } else {
      transparms["g"] <- g
    }
  }

  return(transparms)
}

backtransform_odeparms <- function(transparms, mkinmod,
                                   transform_rates = TRUE,
                                   transform_fractions = TRUE)
{
  # We need the model specification for the names of the model 
  # variables and the information on the sink
  spec = mkinmod$spec

  # Set up container for backtransformed parameters
  parms <- numeric(0)

  # Do not transform initial values for state variables
  state.ini.optim <- transparms[grep("_0$", names(transparms))]
  parms[names(state.ini.optim)] <- state.ini.optim

  # Exponential transformation for rate constants
  if(transform_rates) {
    trans_k <- transparms[grep("^log_k_", names(transparms))]
    trans_k__iore <- transparms[grep("^log_k__iore_", names(transparms))]
    trans_k = c(trans_k, trans_k__iore)
    if (length(trans_k) > 0) {
      k_names <- gsub("^log_k", "k", names(trans_k))
      parms[k_names] <- exp(trans_k)
    }
  } else {
    trans_k <- transparms[grep("^k_", names(transparms))]
    parms[names(trans_k)] <- trans_k
    trans_k__iore <- transparms[grep("^k__iore_", names(transparms))]
    parms[names(trans_k__iore)] <- trans_k__iore
  }

  # Do not transform exponents in IORE models
  N <- transparms[grep("^N", names(transparms))]
  parms[names(N)] <- N

  # Go through state variables and apply inverse isometric logratio transformation
  mod_vars = names(spec)
  for (box in mod_vars) {
    # Get the names as used in the model
    f_names = grep(paste("^f", box, sep = "_"), mkinmod$parms, value = TRUE)
    # Get the formation fraction parameters
    trans_f = transparms[grep(paste("^f", box, sep = "_"), names(transparms))]
    if (length(trans_f) > 0) {
      if(transform_fractions) {
        f <- invilr(trans_f)
        if (spec[[box]]$sink) {
          parms[f_names] <- f[1:length(f)-1]
        } else {
          parms[f_names] <- f
        }
      } else {
        parms[names(trans_f)] <- trans_f
      }
    }
  }

  # Transform parameters also for FOMC, DFOP and HS models
  for (pname in c("alpha", "beta", "k1", "k2", "tb")) {
    if (transform_rates) {
      pname_trans = paste0("log_", pname)
      if (!is.na(transparms[pname_trans])) {
        parms[pname] <- exp(transparms[pname_trans])
      }
    } else {
      if (!is.na(transparms[pname])) {
        parms[pname] <- transparms[pname]
      }
    }
  }

  # DFOP parameter g is treated as a fraction
  if (!is.na(transparms["g_ilr"])) {
    g_ilr <- transparms["g_ilr"]
    parms["g"] <- invilr(g_ilr)[1]
  }
  if (!is.na(transparms["g"])) {
    parms["g"] <- transparms["g"]
  }

  return(parms)
}
# vim: set ts=2 sw=2 expandtab:
