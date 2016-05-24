# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"list.gam.settings"<-
  function(object)
{
### Displays the name of the current database.
  if(exists("object")) {
    cat(paste("\nThe current run number is ", object@Runno, ".\n\n", sep=""))

    if(!any(is.null(object@Prefs@Xvardef$parms)))
      cat("Parameters:",object@Prefs@Xvardef$parms,fill=60)
    
    if(!any(is.null(object@Prefs@Xvardef$covariates))) {
      cat("Covariates:",object@Prefs@Xvardef$covariates,fill=60)
      conts <- cats <- character(0)
      for(i in xvardef("covariates", object))
        if(!is.factor(object@Data[[i]])) {
          if(length(conts))
            conts <- c(conts,i)
          else
            conts <- i
        } else {
          if(length(cats))
            cats <- c(cats,i)
          else
            cats <- i
        }
      cat("  ( Continuous:",conts,")",fill=60)
      cat("  ( Categorical:",cats,")",fill=60)
    }

    if(!any(is.null(object@Prefs@Gam.prefs$disp))){
      cat("Use a dispersion factor (null/true): ",
          object@Prefs@Gam.prefs$disp,"\n")
    }

    if(!any(is.null(object@Prefs@Gam.prefs$steppit))){
      cat("Use stepwise search for covariates (true/false): ",
          object@Prefs@Gam.prefs$steppit,"\n")
    }
    

    if(!any(is.null(object@Prefs@Gam.prefs$onlyfirst))){
      cat("Use only the first value in each individual: ",
          object@Prefs@Gam.prefs$onlyfirst,"\n")
    }

    if(!any(is.null(object@Prefs@Subset))) 
      cat("Subset:",object@Prefs@Subset,"\n")

    if(!any(is.null(object@Prefs@Gam.prefs$start.mod))){
      cat("Starting model: ",
          object@Prefs@Gam.prefs$start.mod,"\n")
    }

    cat("Normalize to median: ",TRUE,"\n")
    
    if(!any(is.null(object@Prefs@Gam.prefs$plot.ids))){
      cat("Use ID labels in GAM plots: ",
          object@Prefs@Gam.prefs$plot.ids,"\n")
    }

  }
  else {
    cat("The current run number is", object@Runno, 
        "but no matching database was found.\n")
  }
  return()
}
