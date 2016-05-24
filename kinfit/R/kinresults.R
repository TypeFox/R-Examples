# $Id: kinresults.R 124 2011-11-09 10:58:27Z jranke $

# Copyright (C) 2008-2011 Johannes Ranke
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

kinresults <- function(kinfits, alpha = 0.05, DTmax = 1000, SFORB=TRUE)
{
	kindata <- data.frame(
    		t = kinfits[[1]]$model$t, 
    		parent = kinfits[[1]]$model$parent)
 	kindata.means <- aggregate(kindata, list(kindata$t), mean)
	kindata.means.mean <- mean(kindata.means$parent, na.rm=TRUE)
	n.times <- length(kindata.means$parent)
	parms <- list()
	df <- err.min <- R2 <- RSS <- TSS <- RSS.means <- TSS.means <- vector()
	DT50 <- DT90 <- vector()
	f <- list()
	
	for (kinmodel in names(kinfits))
	{
		m = kinfits[[kinmodel]]
		if(class(m) == "nls") {
			kindata.means$est <- predict(m, kindata.means)
			
      		if (!"parent.0" %in% names(coef(m))) {
				
       			parms[[kinmodel]] <- switch(kinmodel,
          			SFO = list(k = coef(m)[["k"]]),
          			FOMC = list(alpha = coef(m)[["alpha"]],
              			beta = coef(m)[["beta"]]),
	    			HS = list(k1 = coef(m)[["k1"]],
              			k2 = coef(m)[["k2"]], 
              			tb = coef(m)[["tb"]]),
          			DFOP= list(k1 = coef(m)[["k1"]],
              			k2 = coef(m)[["k2"]], 
              			g = coef(m)[["g"]]))
      		} else {
        			parms[[kinmodel]] <- switch(kinmodel,
          				SFO = list(parent.0 = coef(m)[["parent.0"]], 
              				k = coef(m)[["k"]]),
          				FOMC = list(parent.0 = coef(m)[["parent.0"]],
              				alpha = coef(m)[["alpha"]],
             				 beta = coef(m)[["beta"]]),
          				HS = list(parent.0 = coef(m)[["parent.0"]], 
              				k1 = coef(m)[["k1"]],
             			 	k2 = coef(m)[["k2"]], 
              				tb = coef(m)[["tb"]]),
          				DFOP = list(parent.0 = coef(m)[["parent.0"]],
              				k1 = coef(m)[["k1"]],
              				k2 = coef(m)[["k2"]], 
              				g = coef(m)[["g"]]))
        			if(kinmodel == "DFOP" & SFORB) {
          				k1 = coef(m)[["k1"]]
          				k2 = coef(m)[["k2"]]
          				g = coef(m)[["g"]]
          				parms[["SFORB"]] = 
            				list(parent.0 = coef(m)[["parent.0"]],
              				k1out = g * k1 + (1 - g) * k2,
             			 	k21 = k1 * k2 / (g * k1 + (1 - g) * k2),
              				k12 = (g * (1 - g) * (k1 - k2)^2) / (g * k1 + (1 - g) * k2))
        			}
      		}

  			n.parms = length(coef(m))
			if (!"parent.0" %in% names(coef(m))) {

				f[[kinmodel]] = switch(kinmodel,
					HS = function(t, x) {
						(HS(t, kinfits[[kinmodel]]$model[["parent.0.user"]], 
							coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["tb"]]) - 
						(1 - x/100) * kinfits[[kinmodel]]$model[["parent.0.user"]])^2
					},
					DFOP = function(t, x) {
						(DFOP(t, kinfits[[kinmodel]]$model[["parent.0.user"]], 
							coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["g"]]) - 
						(1 - x/100) * kinfits[[kinmodel]]$model[["parent.0.user"]])^2
					}
				)
			}else{
				f[[kinmodel]] = switch(kinmodel,
					HS = function(t, x) {
						(HS(t, coef(m)[["parent.0"]], 
							coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["tb"]]) - 
						(1 - x/100) * coef(m)[["parent.0"]])^2
					},
					DFOP = function(t, x) {
						(DFOP(t, coef(m)[["parent.0"]], 
							coef(m)[["k1"]], coef(m)[["k2"]], coef(m)[["g"]]) - 
						(1 - x/100) * coef(m)[["parent.0"]])^2
					}
				)
			}
				
			coef(m)

			df[[kinmodel]] = n.times - n.parms
			RSS[[kinmodel]] = sum(summary(m)$residuals^2)
      TSS[[kinmodel]] = sum((m$model$parent - mean(m$model$parent))^2)
	
			DT50.o = switch(kinmodel,
        SFO = log(2)/coef(m)[["k"]],
				FOMC = coef(m)[["beta"]] * (2^(1/coef(m)[["alpha"]]) - 1),
				HS = optimize(f[[kinmodel]], c(0, DTmax), x=50)$minimum,
				DFOP = optimize(f[[kinmodel]], c(0, DTmax), x=50)$minimum)
			
			
      DT50[[kinmodel]] = ifelse(abs(DT50.o - DTmax) < 0.1, NA, DT50.o)
      DT90.o = switch(kinmodel,
          SFO = log(10)/coef(m)[["k"]],
          FOMC = coef(m)[["beta"]] * (10^(1/coef(m)[["alpha"]]) - 1),
          HS = optimize(f[[kinmodel]], c(0, DTmax), x=90)$minimum,
          DFOP = optimize(f[[kinmodel]], c(0, DTmax), x=90)$minimum)
      DT90[[kinmodel]] = ifelse(abs(DT90.o - DTmax) < 0.1, NA, DT90.o)

      # Chi2 error level as defined in FOCUS kinetics (see ?kinerrmin)
      err.min[[kinmodel]] <- kinerrmin(kinfits, kinmodel)

      # Coefficient of determination calculated from residual sum of squares and totals sum of squares
      # so this r2 is what is called model efficiency in FOCUS kinetics (2006), p. 99
      R2[[kinmodel]] = 1 - RSS[[kinmodel]]/TSS[[kinmodel]]

		}
	}

	stats <- data.frame(n.times = n.times, df = df, 
   	mean.means = kindata.means.mean, 
		RSS = RSS, err.min = err.min, R2 = R2)
	results <- data.frame(DT50 = DT50, DT90 = DT90)
	list(parms = parms, stats = stats, results = results)
}
