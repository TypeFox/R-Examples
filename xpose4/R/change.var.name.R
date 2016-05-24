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

"change.var.name" <- function(object, classic = FALSE)
{
  
  data <- object
  cat("The database contains the following items:\n")
  cat(names(object@Data), fill = 60)

  cat("\nType the name of the variable you want to change the name of:\n")
  ans <- readline()
  
  if(is.na(index <- match(ans, names(data@Data)))) {
    cat("The current database does not contain a variable with that name.\n")
    return(cat(""))
  }
  
  if(is.na(sindex <- match(ans, names(data@SData)))) {
    cat("The simulated data does not contain a variable with that name.\n")
  }
  
  lindex <- match(ans, names(data@Prefs@Labels)) 
  
  cat("Type the new name:\n")
  ans1 <- readline()
  ans1 <- paste(ans1,collapse=" ")
  
  if(!is.na(match(ans1,names(data@Data)))) {
    cat("There is already a variable with that name in the current database!\n")
    return(cat(""))
  }

  names(data@Data)[index] <- ans1       # Data
  names(data@SData)[sindex] <- ans1       # SData
  names(data@Prefs@Labels)[lindex] <- ans1  # Labels
  
  data@Prefs@Xvardef$parms <- replace(data@Prefs@Xvardef$parms, grep(ans, data@Prefs@Xvardef$parms), ans1)  # replace in parms
  data@Prefs@Xvardef$covariates <- replace(data@Prefs@Xvardef$covariates, grep(ans, data@Prefs@Xvardef$covariates), ans1)  # replace in parms
  data@Prefs@Xvardef$id <- replace(data@Prefs@Xvardef$id, grep(ans, data@Prefs@Xvardef$id), ans1)
  data@Prefs@Xvardef$idv <- replace(data@Prefs@Xvardef$idv, grep(ans, data@Prefs@Xvardef$idv), ans1)
  data@Prefs@Xvardef$idlab <- replace(data@Prefs@Xvardef$idlab, grep(ans, data@Prefs@Xvardef$idlab), ans1)
  data@Prefs@Xvardef$occ <- replace(data@Prefs@Xvardef$occ, grep(ans, data@Prefs@Xvardef$occ), ans1)
  data@Prefs@Xvardef$dv <- replace(data@Prefs@Xvardef$dv, grep(ans, data@Prefs@Xvardef$dv), ans1)
  data@Prefs@Xvardef$res <- replace(data@Prefs@Xvardef$res, grep(ans, data@Prefs@Xvardef$res), ans1)
  data@Prefs@Xvardef$iwres <- replace(data@Prefs@Xvardef$iwres, grep(ans, data@Prefs@Xvardef$iwres), ans1)
  data@Prefs@Xvardef$wres <- replace(data@Prefs@Xvardef$wres, grep(ans, data@Prefs@Xvardef$wres), ans1)
  data@Prefs@Xvardef$pred <- replace(data@Prefs@Xvardef$pred, grep(ans, data@Prefs@Xvardef$pred), ans1)
  data@Prefs@Xvardef$ipred <- replace(data@Prefs@Xvardef$ipred, grep(ans, data@Prefs@Xvardef$ipred), ans1)
  
  for (i in 1:length(data@Prefs@Labels)) {
    data@Prefs@Labels[i] <- replace(data@Prefs@Labels[i], grep(ans, data@Prefs@Labels[i]), ans1)
  }
  
  #class(object@Data) <- "data.frame"
  #vname(data@Data[,index]) <- ans1
  #class(object@Data) <- "xpose"
            
        if (classic==TRUE) {
          c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
          eval(c1)
          c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
          eval(c2)
          return(cat(""))
          
        } else {
          return(data)
        }
}
