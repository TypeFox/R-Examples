# $Id: kinplot.R 117 2011-06-14 08:52:14Z kati $

# Copyright (C) 2008-2013 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

# This file is part of the R package kinfit

# kinfit is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

if(getRversion() >= '2.15.1') utils::globalVariables("x")

kinplot <- function(kinobject, 
	main = "",
	xlab = "Time [days]", ylab = "Parent [% of applied radioactivity]",
        ylim = c("auto", "auto"),
	lpos = "topright")
{
	kindata <- na.omit(kinobject$data)
	kinfits <- kinobject$fits
  if (ylim[1] == "auto") ylim[1] <- 0
  if (ylim[2] == "auto") ylim[2] <- max(kindata$parent)
  ylim <- as.numeric(ylim)
        
	plot(kindata$t, kindata$parent,
	  main = main,
	  xlab = xlab,
	  ylab = ylab,
	  ylim = ylim
  )
	n.m <- length(kinfits)
	colors <- ltys <- 1:n.m
	names(colors) <- names(ltys) <- names(kinfits)
  ltext <- paste(kinobject$parent, "measured")
	for (kinmodel in names(kinfits))
	{
		m = kinfits[[kinmodel]]
		if(class(m) == "nls") {
      if (!"parent.0" %in% names(coef(m))) {
        switch(kinmodel,
          SFO = lines(
              t <- seq(min(kindata$t), max(kindata$t), length.out=500),
              predict(m, 
              newdata = data.frame(t)),
              col = colors[[kinmodel]],
              lty = ltys[[kinmodel]]),
          FOMC = lines(
              t <- seq(min(kindata$t), max(kindata$t), length.out=500),
              predict(m, 
              newdata = data.frame(t)),
              col = colors[[kinmodel]],
              lty = ltys[[kinmodel]]), 
 	    HS = lines(
              t <- seq(min(kindata$t), max(kindata$t), length.out=500),
              predict(m, 
              newdata = data.frame(t)),
              col = colors[[kinmodel]],
              lty = ltys[[kinmodel]]), 
          DFOP = lines(
              t <- seq(min(kindata$t), max(kindata$t), length.out=500),
              predict(m, 
              newdata = data.frame(t)),
              col = colors[[kinmodel]],
              lty = ltys[[kinmodel]])
		)
        ltext <- c(ltext, paste("Fitted", kinmodel, "model"))
      } else {
        switch(kinmodel,
          SFO = curve(SFO(x, 
            coef(m)[["parent.0"]], 
            coef(m)[["k"]]),
            from = min(kindata$t), to = max(kindata$t), add=TRUE,
            col = colors[[kinmodel]],
            lty = ltys[[kinmodel]]),
          FOMC = curve(FOMC(x, 
            coef(m)[["parent.0"]],
            coef(m)[["alpha"]],
            coef(m)[["beta"]]),
            from = min(kindata$t), to = max(kindata$t), add=TRUE,
            col = colors[[kinmodel]],
            lty = ltys[[kinmodel]]),
          HS = curve(HS(x, 
            coef(m)[["parent.0"]], 
            coef(m)[["k1"]],
            coef(m)[["k2"]],
            coef(m)[["tb"]]),
            from = min(kindata$t), to = max(kindata$t), add=TRUE,
            col = colors[[kinmodel]],
            lty = ltys[[kinmodel]]),
          DFOP = curve(DFOP(x, 
            coef(m)[["parent.0"]], 
            coef(m)[["k1"]],
            coef(m)[["k2"]],
            coef(m)[["g"]]),
            from = min(kindata$t), to = max(kindata$t), add=TRUE,
            col = colors[[kinmodel]],
            lty = ltys[[kinmodel]]))
        ltext <- c(ltext, paste("Fitted", kinmodel, "model"))
      }
		} else {
        ltext <- c(ltext, paste(kinmodel, "model failed"))
            ltys[[kinmodel]] <- NA
		} 
	}
	legend(lpos, bty="n", inset = 0.05, 
		legend = ltext,
		pch = c(1, rep(NA, n.m)),
		lty = c(NA, ltys),
		col = c(1, colors))
}
