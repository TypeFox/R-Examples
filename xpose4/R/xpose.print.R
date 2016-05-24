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

"xpose.print" <- function(object ,long=TRUE)
{
  cat("The database contains the following observed items:\n")
  cat(names(object@Data),fill=60)
  
  if(!any(is.null(object@SData))) {
    cat("\nThe database contains the following simulated items:\n")
    cat(names(object@SData),fill=60)
  }

  cat("\nThe following variables are defined:\n\n")

  if(!any(is.null(object@Prefs@Xvardef$id)))
    cat("ID variable:",object@Prefs@Xvardef$id,"\n")
  
  if(!any(is.null(object@Prefs@Xvardef$idlab)))
    cat("Label variable:",object@Prefs@Xvardef$idlab,"\n")

#  if(!any(is.null(flag(data)))) {
#    cat("Flag variable:",vname(flag(data)),"\n")
#    if(!any(is.null(cur.flag(data)))) 
#      cat("Current value of flag:",cur.flag(data),"\n")
#  }

  if(!any(is.null(object@Prefs@Xvardef$idv)))
    cat("Independent variable:",object@Prefs@Xvardef$idv,"\n")
  
  if(!any(is.null(object@Prefs@Xvardef$occ)))
    cat("Occasion variable:",object@Prefs@Xvardef$occ,"\n")
  
  
 if(!any(is.null(object@Prefs@Xvardef$dv))) {

   if(is.factor(object@Prefs@Xvardef$dv)) {
     cat("Dependent variable (categorical):",object@Prefs@Xvardef$dv,"\n")
   } else {
     cat("Dependent variable:",object@Prefs@Xvardef$dv,"\n")
   }
 }
 
  if (long) {
    if(!any(is.null(object@Prefs@Xvardef$pred)))
      cat("Population prediction variable:",object@Prefs@Xvardef$pred,"\n")
    
    if(!any(is.null(object@Prefs@Xvardef$ipred)))
      cat("Individual prediction variable:",object@Prefs@Xvardef$ipred,"\n")
    
    if(!any(is.null(object@Prefs@Xvardef$wres)))
      cat("Weighted population residual variable:",object@Prefs@Xvardef$wres,"\n")
    
    if(!any(is.null(object@Prefs@Xvardef$iwres)))
      cat("Weighted individual residual variable:",object@Prefs@Xvardef$iwres,"\n")
    
    if(!any(is.null(object@Prefs@Xvardef$res)))
      cat("Population residual variable:",object@Prefs@Xvardef$res,"\n")
  }
    
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

  if(!any(is.null(object@Prefs@Xvardef$tvparms)))  {
    cat("Typical parameters:",object@Prefs@Xvardef$tvparms,fill=60)
  }

  if(!any(is.null(object@Prefs@Xvardef$ranpar)))  {
    cat("Variability parameters:",object@Prefs@Xvardef$ranpar,fill=60)
  }
  
  if(!any(is.null(object@Prefs@Miss))) 
    cat("Missing value label:",object@Prefs@Miss,"\n")
    
  if(!any(is.null(object@Prefs@Subset))) 
    cat("Subset:",object@Prefs@Subset,"\n")

  invisible()
  return(cat(""))
}
