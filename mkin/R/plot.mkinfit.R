# Copyright (C) 2010-2015 Johannes Ranke
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
if(getRversion() >= '2.15.1') utils::globalVariables(c("type", "variable", "observed"))

plot.mkinfit <- function(x, fit = x,
  obs_vars = names(fit$mkinmod$map),
  xlab = "Time", ylab = "Observed",
  xlim = range(fit$data$time), 
  ylim = "default",
  col_obs = 1:length(fit$mkinmod$map),
  pch_obs = col_obs, 
  lty_obs = rep(1, length(fit$mkinmod$map)),
  add = FALSE, legend = !add, 
  show_residuals = FALSE, maxabs = "auto",
  lpos = "topright", inset = c(0.05, 0.05), ...)
{
  if (add && show_residuals) stop("If adding to an existing plot we can not show residuals")

  if (ylim[[1]] == "default") {
    ylim = c(0, max(subset(fit$data, variable %in% obs_vars)$observed, na.rm = TRUE))
  }

  solution_type = fit$solution_type
  parms.all <- c(fit$bparms.optim, fit$bparms.fixed)

  ininames <- c(
    rownames(subset(fit$start, type == "state")),
    rownames(subset(fit$fixed, type == "state")))
  odeini <- parms.all[ininames]

  # Order initial state variables
  names(odeini) <- sub("_0", "", names(odeini))
  odeini <- odeini[names(fit$mkinmod$diffs)]

  outtimes <- seq(xlim[1], xlim[2], length.out=100)

  odenames <- c(
    rownames(subset(fit$start, type == "deparm")),
    rownames(subset(fit$fixed, type == "deparm")))
  odeparms <- parms.all[odenames]

  out <- mkinpredict(fit$mkinmod, odeparms, odeini, outtimes, 
          solution_type = solution_type, atol = fit$atol, rtol = fit$rtol)

  # Set up the plot if not to be added to an existing plot
  if (add == FALSE) {
    if (show_residuals) {
      oldpar <- par(no.readonly = TRUE)
      layout(matrix(c(1, 2), 2, 1), heights = c(2, 1.3))
      par(mar = c(3, 4, 4, 2) + 0.1)
    }
    plot(0, type="n", 
      xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, ...)
  }
  # Plot the data and model output
  names(col_obs) <- names(pch_obs) <- names(lty_obs) <- names(fit$mkinmod$map)
  for (obs_var in obs_vars) {
    points(subset(fit$data, variable == obs_var, c(time, observed)), 
      pch = pch_obs[obs_var], col = col_obs[obs_var])
  }
  matlines(out$time, out[obs_vars], col = col_obs[obs_vars], lty = lty_obs[obs_vars])
  if (legend == TRUE) {
    legend_names = lapply(names(fit$mkinmod$spec), function(x) {
                          if (!is.null(fit$mkinmod$spec[[x]]$full_name))
                            if (is.na(fit$mkinmod$spec[[x]]$full_name)) x
                            else fit$mkinmod$spec[[x]]$full_name
                          else x
      })
    legend(lpos, inset= inset, legend = legend_names,
      col = col_obs[obs_vars], pch = pch_obs[obs_vars], lty = lty_obs[obs_vars])
  }
  # Show residuals if requested
  if (show_residuals) {
    par(mar = c(5, 4, 0, 2) + 0.1)
    residuals <- subset(fit$data, variable %in% obs_vars, residual)
    if (maxabs == "auto") maxabs = max(abs(residuals), na.rm = TRUE)
    plot(0, type="n", 
      xlim = xlim, 
      ylim = c(-1.2 * maxabs, 1.2 * maxabs),
      xlab = xlab, ylab = "Residuals")
    for(obs_var in obs_vars){
      residuals_plot <- subset(fit$data, variable == obs_var, c("time", "residual"))
      points(residuals_plot, pch = pch_obs[obs_var], col = col_obs[obs_var])
    }
    abline(h = 0, lty = 2)
    par(oldpar, no.readonly = TRUE)
  }
}
