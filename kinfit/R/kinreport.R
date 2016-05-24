# $Id: kinreport.R 123 2011-11-01 12:26:41Z jranke $

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

kinreport <- function(kinobject, file = NA, data = TRUE, R2 = FALSE, vcov = FALSE, endpoint.digits = 1)
{
	if (!is.na(file)) {
		sink(file, split=TRUE)
	}

	cat("Parent compound: ", kinobject$parent, "\n")
  if (!is.null(kinobject$label)) {
    cat("Label position:\t\t", kinobject$label, "\n")
  }
	cat("Study type:      ", kinobject$type, "\n")
	cat("System:          ", kinobject$system, "\n")
	if (!is.null(kinobject$source)) {
    cat("Source:          ", kinobject$source, "\n")
  }
  cat("kinfit version:  ", as.character(packageVersion("kinfit")), "\n")
  cat("R version:       ", paste(R.version$major, R.version$minor, sep="."), "\n")
  cat("Report generated:", date(), "\n")
	cat("\n")
  if (data) {
    cat("Data:\n")
    print(kinobject$data)
    cat("\n")
  }
	fit.names <- names(kinobject$fits)
	for (kinmodel in fit.names)
	{
    m <- kinobject$fits[[kinmodel]]
    if (!(class(m) == "try-error")) {
      cat("\n\n---\n")
      cat("Nonlinear least squares fit of the", kinmodel, "model\n")
      if (!"parent.0" %in% names(coef(m))) {
        cat(paste("Initial value of parent fixed to ", m$model$parent.0.user, "\n", sep=""))
      }
      cat("\n")
      cat("Parameter estimation:\t")
      s <- summary(m)
      df <- s$df[2]
      p <- 1 - pt(s$parameters[,3], df = df)
      parms <- matrix(nrow = nrow(s$parameters), ncol=4)
      dimnames(parms) = list(rownames(s$parameters), 
        c("Estimate", "Std. Error", "t value", "Pr(>t)"))
      parms[, c(1,2,3)] <- s$parameters[,c(1,2,3)]
      parms[, 4] <- p
      cat("\n")
      print(parms, digits=3)
      cat("\n")
      if(vcov)
      {
        cat("Variance-covariance matrix:\n")
        print(vcov(m))
        cat("\n")
      }
      cat("Chi2 error estimation: ", 
        round(100 * kinobject$results$stats[kinmodel, "err.min"], digits=2), 
          " %\n", sep="")
      cat("\n")
      if(R2)
      {
        cat("Coefficient of determination R2: ",
          round(kinobject$results$stats[kinmodel, "R2"], digits=3), "\n")
      }
    }
	}
	cat("\n\n---\n")
	cat("Endpoint estimates\n\n")
	print(round(kinobject$results$results, digits=endpoint.digits))

	if (!is.na(file)) sink()
}
