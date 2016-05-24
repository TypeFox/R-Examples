# Copyright (C) 2015 Johannes Ranke
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

plot.mmkin <- function(x, main = "auto", legends = 1, errmin_var = "All data", errmin_digits = 2, 
                       cex = 0.7, rel.height.middle = 0.9, ...) {
  n.m <- nrow(x)
  n.d <- ncol(x)
  if (n.m > 1 & n.d > 1) stop("Please select fits either for one model or for one dataset")
  if (n.m == 1 & n.d == 1) loop_over = "none"
  if (n.m > 1) loop_over <- "models"
  if (n.d > 1) loop_over <- "datasets"
  n.fits <- length(x)

  if (main == "auto") {
    main = switch(loop_over,
                  none = paste(rownames(x), colnames(x)),
                  models = colnames(x),
                  datasets = rownames(x))
  }

  oldpar <- par(no.readonly = TRUE)
  rel.heights <- if (n.fits > 2) c(1, rep(rel.height.middle, n.fits - 2), 1)
                 else rep(1, n.fits)
  layout(matrix(1:(2 * n.fits), n.fits, 2, byrow = TRUE), heights = rel.heights)

  #par(mfrow = c(n.fits, 2))
  par(mar = c(3.0, 4.1, 4.1, 2.1)) # Reduce bottom margin by 2.1 - hides x axis legend
  par(cex = cex)

  for (i.fit in 1:n.fits) {
    if (i.fit == 2) {
      # Reduce top margin by 2 after the first plot as we have no main title, 
      # reduced plot height, therefore we need rel.height.middle in the layout
      par(mar = c(3.0, 4.1, 2.1, 2.1))
    }
    if (i.fit == n.fits) {
      # Reduce top margin by 2 after the first plot as we have no main title, 
      # plot height remains about constant
      par(mar = c(5.1, 4.1, 2.1, 2.1))

    }
    fit <- x[[i.fit]]
    plot(fit, legend = legends == i.fit, ...)

    title(main, outer = TRUE, line = -2)

    fit_name <- switch(loop_over,
                       models = rownames(x)[i.fit],
                       datasets = colnames(x)[i.fit],
                       none = "") 

    if (!is.null(errmin_var)) {
      chi2 <- paste0(round(100 * mkinerrmin(fit)[errmin_var, "err.min"], errmin_digits), "%")
      mtext(bquote(.(fit_name) ~ chi^2 ~ "error level" == .(chi2)), cex = cex, line = 0.4)
    }
    mkinresplot(fit, legend = FALSE, ...)
    mtext(paste(fit_name, "residuals"), cex = cex, line = 0.4)
  }

  par(oldpar, no.readonly = TRUE)
}
