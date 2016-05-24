####**********************************************************************
####**********************************************************************
####
####  SPIKE AND SLAB 1.1.5
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Written and Developed by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    J. Sunil Rao, Ph. D.
####    Professor and Director of the Division of Biostatistics, 
####    Department of Epidemiology & Public Health
####    Clinical Research Bldg, R-669
####    1120 NW 14th Street, Room 1056
####    Miami, FL 33136
####    email:  rao.jsunil@gmail.com
####    URL:    http://biostat.med.miami.edu/people/primary-faculty/sunil-rao
####  ----------------------------------------------------------------
####  Maintained by:
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

print.spikeslab <- function(x, ...)
{

  ### Check that object is compatible
  if (!inherits(x, "spikeslab"))
     stop("This function only works for objects of class `spikeslab'")
  
  ###check whether object inherits mixing type
  ###make suitable alterations to merge mixing results
  if (sum(inherits(x, c("spikeslab", "mixing"), TRUE) == c(1, 2)) == 2) {
     o.r <- order(abs(x$spikeslab.obj$bma), decreasing = TRUE)
     grr.mix <- x$grr.mix[o.r]
     grr.mix.scale <- x$grr.mix.scale[o.r]
     new.summary <- cbind(x$spikeslab.obj$summary, grr.mix, grr.mix.scale)
     new.summary <- new.summary[, c(1, 2, 5, 3, 4, 6)]
     x <- x$spikeslab.obj
     x$summary <- new.summary
  }

  ### extract summary data
  verbose <- x$verbose
  
  ### --------------------------------------------------------------
  ###	Terminal Output
  ### --------------------------------------------------------------	

  if (sum(inherits(x, c("spikeslab", "cv"), TRUE) == c(1, 2)) != 2) {

    ### spike slab basic output
    
    cat("-------------------------------------------------------------------","\n")
    cat("Variable selection method     :",verbose[[1]],"\n")
    cat("Big p small n                 :",verbose[[2]],"\n")
    cat("Screen variables              :",verbose[[3]],"\n")
    cat("Fast processing               :",verbose[[4]],"\n")
    cat("Sample size                   :",verbose[[5]],"\n")
    cat("No. predictors                :",verbose[[6]],"\n")
    cat("No. burn-in values            :",verbose[[7]],"\n")
    cat("No. sampled values            :",verbose[[8]],"\n")
    cat("Estimated mse                 :",verbose[[9]],"\n")
    cat("Model size                    :",verbose[[10]],"\n")
    cat("\n\n")
    cat("---> Top variables:\n")
    print(round(x$summary[x$summary[, 2] != 0, ], 3))
    
  }
  
  if (sum(inherits(x, c("spikeslab", "cv"), TRUE) == c(1, 2)) == 2) {
    
    ### cv output
    
    cat("-------------------------------------------------------------------","\n")
    cat("Variable selection method     :",verbose[[1]],"\n")
    cat("Big p small n                 :",verbose[[2]],"\n")
    cat("Screen variables              :",verbose[[3]],"\n")
    cat("Fast processing               :",verbose[[4]],"\n")
    cat("Sample size                   :",verbose[[5]],"\n")
    cat("No. predictors                :",verbose[[6]],"\n")
    cat("No. burn-in values            :",verbose[[7]],"\n")
    cat("No. sampled values            :",verbose[[8]],"\n")
    cat("K-fold                        :",verbose[[9]],"\n")
    cat("CV mean-squared error         :",verbose[[10]],"\n")
    cat("Model size                    :",verbose[[11]],"\n")
    cat("\n\nTop variables in terms of stability:\n")
    print(head(round(x$stability, 3), verbose[[12]]))
    
  }

  cat("-------------------------------------------------------------------","\n")

}
