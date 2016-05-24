# $Id: kinfit.R 116 2011-06-14 08:46:47Z kati $

# Copyright (C) 2008-2010 Johannes Ranke
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

kinfit <- function(kindata, kinmodels = c("SFO"), 
	parent.0.user = NA, 
	parent.0.fixed = FALSE, 
	start.SFO = list(parent.0 = NA, k = NA), 
	start.FOMC = list(parent.0 = NA, alpha = NA, beta = NA), 
	start.DFOP = list(parent.0 = NA, k1 = NA, k2 = NA, g = NA),
	start.HS = list(parent.0 = NA, k1 = NA, k2 = NA, tb = NA),
        algorithm = "default")
{
	
	kindata <- subset(kindata, !is.na(kindata$parent))
	kinfits <- list()

	if (!is.na(parent.0.user)) {
		start.SFO$parent.0 = parent.0.user
		start.FOMC$parent.0 = parent.0.user
		start.DFOP$parent.0 = parent.0.user
		start.HS$parent.0 = parent.0.user
	}

	lmlogged = lm(log(parent) ~ t, data = kindata)
      k.est = -coef(lmlogged)[["t"]]

	for (kinmodel in kinmodels)
	{

		if (kinmodel == "SFO") {
			if (is.na(start.SFO$parent.0)) {
        			start.SFO$parent.0 = max(kindata$parent)
			}
			if (is.na(start.SFO$k)) {
				start.SFO$k = - coef(lmlogged)[["t"]]
			}
      if (parent.0.fixed)
        {
          start.SFO = start.SFO[-1]
          kinfits[[kinmodel]] = try(
            nls(parent ~ SFO(t, parent.0.user, k),
              data = kindata, model = TRUE,
              start = start.SFO,
              algorithm = algorithm), silent=TRUE)
        } else {
          kinfits[[kinmodel]] = try(
            nls(parent ~ SFO(t, parent.0, k),
              data = kindata, model = TRUE,
              start = start.SFO,
              algorithm = algorithm), silent=TRUE)
        }
      k.est = coef(kinfits$SFO)[["k"]]
		}	
		if (kinmodel == "FOMC") {
			if (is.na(start.FOMC$parent.0)) {
       			start.FOMC$parent.0 = max(kindata$parent)
			}
			if (is.na(start.FOMC$alpha)) {
				start.FOMC$alpha = 1
			}
			if (is.na(start.FOMC$beta)) {
				start.FOMC$beta = start.FOMC$alpha / k.est 
			}
			
      		if (parent.0.fixed)
      		{
        			start.FOMC = list(alpha = start.FOMC$alpha, beta = start.FOMC$beta)
				
        			kinfits[[kinmodel]] = try(
         			nls(parent ~ FOMC(t, parent.0.user, alpha, beta),
           				 data = kindata, model = TRUE,
          				  start = start.FOMC,
            			algorithm = algorithm), silent=TRUE)

      		} else {
        			kinfits[[kinmodel]] = try(
          			nls(parent ~ FOMC(t, parent.0, alpha, beta),
            		data = kindata, model = TRUE,
           			start = start.FOMC,
           			algorithm = algorithm), silent=TRUE)

		      }
		}	
		if (kinmodel == "DFOP") {
			if (is.na(start.DFOP$parent.0)) {
        			start.DFOP$parent.0 = max(kindata$parent)
			}
			if (is.na(start.DFOP$k1)) {
				start.DFOP$k1 = k.est * 2
			}
			if (is.na(start.DFOP$k2)) {
				start.DFOP$k2 = k.est / 2
			}
			if (is.na(start.DFOP$g)) {
				start.DFOP$g = 0.5
			}
			if (parent.0.fixed)
      		{
				start.DFOP = list(k1 = start.DFOP$k1, k2 = start.DFOP$k2, g = start.DFOP$g)

				kinfits[[kinmodel]] = try(
				nls(parent ~ DFOP(t, parent.0.user, k1, k2, g),
					data = kindata, model = TRUE,
					start = start.DFOP,
          				algorithm = algorithm), silent=TRUE)
			}else{
				kinfits[[kinmodel]] = try(
				nls(parent ~ DFOP(t, parent.0, k1, k2, g),
					data = kindata, model = TRUE,
					start = start.DFOP,
          				algorithm = algorithm), silent=TRUE)
			}
		}	
		if (kinmodel == "HS") {
			if (is.na(start.HS$parent.0)) {
        			start.HS$parent.0 = max(kindata$parent)
			}
			if (is.na(start.HS$k1)) {
				start.HS$k1 = k.est
			}
			if (is.na(start.HS$k2)) {
				start.HS$k2 = k.est / 10
			}
			if (is.na(start.HS$tb)) {
				start.HS$tb = 0.05 * max(kindata$t)
			}
			
			if (parent.0.fixed)
      		{		
				
				start.HS = list(k1 = start.HS$k1, k2 = start.HS$k2, tb = start.HS$tb)	

				kinfits[[kinmodel]] = try(
				nls(parent ~ HS(t, parent.0.user, k1, k2, tb),
					data = kindata, model = TRUE,
					start = start.HS,
          				algorithm = algorithm), silent=TRUE)
			}else{
				kinfits[[kinmodel]] = try(
				nls(parent ~ HS(t, parent.0, k1, k2, tb),
					data = kindata, model = TRUE,
					start = start.HS,
          				algorithm = algorithm), silent=TRUE)
			}
		}	
	}
	return(kinfits)		
}
