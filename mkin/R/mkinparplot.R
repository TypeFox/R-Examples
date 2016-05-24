# Copyright (C) 2014 Johannes Ranke
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
mkinparplot <- function(object) {
  state.optim = rownames(subset(object$start, type == "state"))
  deparms.optim = rownames(subset(object$start, type == "deparm"))
  fractions.optim = grep("^f_", deparms.optim, value = TRUE)
  N.optim = grep("^N_", deparms.optim, value = TRUE)
  if ("g" %in% deparms.optim) fractions.optim <- c("g", fractions.optim)
  rates.optim.unsorted = setdiff(deparms.optim, union(fractions.optim, N.optim))
  rates.optim <- rownames(object$start[rates.optim.unsorted, ])
  n.plot <- c(state.optim = length(state.optim), 
              rates.optim = length(rates.optim), 
	      N.optim = length(N.optim),
              fractions.optim = length(fractions.optim))
  n.plot <- n.plot[n.plot > 0]

  oldpars <- par(no.readonly = TRUE)
  layout(matrix(1:length(n.plot), ncol = 1), heights = n.plot + 1)

  s <- summary(object)
  bpar <- data.frame(t(s$bpar[, c("Estimate", "Lower", "Upper")]))
  par(mar = c(2.1, 2.1, 0.1, 2.1))
  par(cex = 1)
  for (type in names(n.plot)) {
    parnames <- get(type)
    values <- bpar[parnames]
    values_with_confints <- data.frame(t(subset(data.frame(t(values)), !is.na("Lower"))))
    xlim = switch(type,
                  state.optim = range(c(0, unlist(values)), 
                                      na.rm = TRUE, finite = TRUE),
                  rates.optim = range(c(0, unlist(values)), 
                                      na.rm = TRUE, finite = TRUE),
                  N.optim = range(c(0, 1, unlist(values)), 
                                      na.rm = TRUE, finite = TRUE),
                  fractions.optim = range(c(0, 1, unlist(values)), 
                                          na.rm = TRUE, finite = TRUE))
    parname_index <- length(parnames):1 # Reverse order for strip chart

    stripchart(values["Estimate", ][parname_index], 
               xlim = xlim,
               ylim = c(0.5, length(get(type)) + 0.5),
               yaxt = "n")
    if (type %in% c("rates.optim", "fractions.optim")) abline(v = 0, lty = 2)
    if (type %in% c("N.optim", "fractions.optim")) abline(v = 1, lty = 2)
    position <- ifelse(values["Estimate", ] < mean(xlim), "right", "left")
    text(ifelse(position == "left", min(xlim), max(xlim)), 
         parname_index, parnames, 
         pos = ifelse(position == "left", 4, 2))

    values.upper.nonInf <- ifelse(values["Upper", ] == Inf, 1.5 * xlim[[2]], values["Upper", ])
    # Suppress warnings for non-existing arrow lengths
    suppressWarnings(arrows(as.numeric(values["Lower", ]), parname_index, 
	   as.numeric(values.upper.nonInf), parname_index, 
	   code = 3, angle = 90, length = 0.05))
  }
  par(oldpars)
}
