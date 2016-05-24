# Copyright (C) 2008-2014 Johannes Ranke
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
if(getRversion() >= '2.15.1') utils::globalVariables(c("variable", "residual"))

mkinresplot <- function (object, 
  obs_vars = names(object$mkinmod$map),
  xlim = c(0, 1.1 * max(object$data$time)),
  xlab = "Time", ylab = "Residual",
  maxabs = "auto", legend= TRUE, lpos = "topright", ...) 
{
	obs_vars_all <- as.character(unique(object$data$variable))

  if (length(obs_vars) > 0){
      obs_vars <- intersect(obs_vars_all, obs_vars)	
  } else obs_vars <- obs_vars_all

  residuals <- subset(object$data, variable %in% obs_vars, residual)

  if (maxabs == "auto") maxabs = max(abs(residuals), na.rm = TRUE)

	col_obs <- pch_obs <- 1:length(obs_vars)
 	names(col_obs) <- names(pch_obs) <- obs_vars

  plot(0, type = "n",
       xlab = xlab, ylab = ylab, 
       xlim = xlim,
       ylim = c(-1.2 * maxabs, 1.2 * maxabs), ...)

	for(obs_var in obs_vars){
		residuals_plot <- subset(object$data, variable == obs_var, c("time", "residual"))
		points(residuals_plot, pch = pch_obs[obs_var], col = col_obs[obs_var])
	}

  abline(h = 0, lty = 2)

  if (legend == TRUE) {
    legend(lpos, inset = c(0.05, 0.05), legend = obs_vars, 
      col = col_obs[obs_vars], pch = pch_obs[obs_vars])
  }
}
