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

"check.gamobj"<-
  function() {
    ## another function definition
    getit <- function()
      {
        cat("\nYou have to specify the parameter name and the run number",
            "of the gam objects you want to plot. The following", 
            "gam objects are available:\n", fill = 60)
        if(.Platform$OS == "windows") {
          cat(objects(pattern="gam.xpose*",pos=1),fill=60)
        } else {
          cat(objects(pattern = "^gam.xpose",pos=1), fill = 60)
        }
        cat("\nParameter (0 to exit): ")
        ans <- readline()
        if(ans == 0){
          return(ans <- NULL)
        }
        cat("Run number (0 to exit):")
        ans1 <- readline()
        if(ans1 == 0){
          return(ans1 <- NULL)
        }
        gobjname <- paste("gam.xpose.", ans, ".", ans1, sep = "")
        if(!exists(gobjname, where = 1)) {
          cat("\n*There are no object that matches", gobjname, 
              "\n")
          gobjname <- Recall()
        }
        return(gobjname)
      }
    ##
    ## The real code starts here
    ##
    if(exists("current.gam", where = 1)) {
      
      cat("\nThe current GAM object is for",
          eval(parse(text=paste("current.gam","$pars",sep=""))),
          #current.gam$pars,
          "in run",
          #current.gam$runno,
          eval(parse(text=paste("current.gam","$runno",sep=""))),
          ".\n")
      cat("\nDo you want to proceed with this gam object? y(n) ")
      ans <- readline()
      if(ans != "y" && ans != "") {
        gobjname <- getit()
        if(!is.null(gobjname)){
          c1 <- call("assign",pos = 1, "current.gam", eval(as.name(gobjname)),immediate=T)
          eval(c1)
        }
      } else {
        gobjname <- T
      }
    }  else {
      gobjname <- getit()
      if(!is.null(gobjname)){
        c2 <- call("assign",pos=1, "current.gam", eval(as.name(gobjname)),immediate=T)
        eval(c2)
      }
    }
    return(gobjname)
  }






