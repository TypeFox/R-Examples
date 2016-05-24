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

"xp.gam"<-
  function(object,
           parnam=xvardef("parms", object),
           covnams = xvardef("covariates", object),
           wts.col=NULL,
           ask.for.input=TRUE,
           overwrite=TRUE, # should be false for classic
           ...){
    

    ## begin function definition    
    ask.for.par <- function(...){
      cat("\nEnter name of parameter to use for this GAM search (0 to exit):")
      ans <- readline()
      if(ans == 0) return(NULL)
      if(length(ans) > 1) {
        cat("\nYou have specified more than one parameter.\n")
        cat("The GAM can be run on only one parameter at a time.\n")
        ans <- Recall(...)
      } else {
        ans.exists <- check.vars(ans,object)
        #cat("The name you typed doesn't match any of\n")
        #cat("the names in the current database\n")
        if(is.null(ans.exists)) ans <- Recall(...)
      }
      return(ans)
    }

    get.par <- function(nams, get.input=FALSE,...){
      ans <- NULL
      if(length(nams)==0) {
        cat("\nNo parameter is defined for this GAM search\n")
        if(get.input){
          ans <- ask.for.par()
        } else {
          cat("\nType '?xp.gam' for more information.\n")
        }
      }
      if(length(nams)>1) {
        cat("\nThere is more than one parameter defined\n")
        cat("for this GAM search. The parameters defined are:\n\n")
        cat(nams, fill = 60)
        cat("\nThe GAM can be run on only one parameter at a time.\n")
        if(get.input) {
          ans <- ask.for.par()
        } else {
          cat("\nType '?xp.gam' for more information.\n")
        }
      }
      if(length(nams)==1) {
        ans <- nams
      }
      return(ans)
    }

    ask.for.wts <- function(...){
      cat("\nWeight column to use (0 to exit, NULL for no weight):")
      ans <- readline()
      if(ans == "NULL") return("NULL")
      if(ans == 0) return(NULL)
      if(length(ans) > 1) {
        cat("\nYou have specified more than one weight.\n")
        cat("Only one weight is allowed.\n")
        ans <- Recall(...)
      } else {
        if(is.na(pmatch(ans,names(object@Data.firstonly)))){
          cat(paste("\n","-----------Variable(s) not defined!-------------\n",
                    ans, "is not defined in the current database\n",
                    "and must be defined for this command to work!\n",
                    "------------------------------------------------\n"))
          ans <- Recall(...)
        }
        return(ans)
      }
    }

    get.wts <- function(nams, get.input=FALSE,...){
      ans <- NULL
      if(length(nams)==0) {
        cat("\nNo weights are defined for this GAM search\n")
        if(get.input){
          ans <- ask.for.wts()
        } else {
          cat("\nType '?xp.gam' and '?xpose.gam' for more information.\n")
        }
      }
      if(length(nams)>1) {
        cat("\nPlease specify a the weights for the parameter.\n")
        cat("The weights come from columns in the data contained\n")
        cat("in the Data.firstonly section of the xpose data object.\n")
        cat("These values usualy come from the .phi file of a NONMEM run.\n")
        cat("Possible weight values (column names) are:\n\n")
        cat(nams, fill = 60)
        cat("\nOnly one weight can be specified.\n")
        if(get.input) {
          ans <- ask.for.wts()
        } else {
          cat("\nType '?xp.gam' and '?xpose.gam' for more information.\n")
        }
      }
      if(length(nams)==1) {
        ans <- nams
      }
      return(ans)
    }
    ## end function definition

    pars <- get.par(parnam,get.input=ask.for.input,...)
    if(is.null(pars)) {
      return(invisible())
    }


    ## check for weighting
    if(object@Prefs@Gam.prefs$wts & ask.for.input){
      wts <- get.wts(names(object@Data.firstonly),get.input=ask.for.input,...)
      if(is.null(wts)) {
        return(invisible())
      }
      if(wts=="NULL") wts <- NULL
      wts.col <- wts
    }
    
    ##
    ## Check if we have an existing GAM objects
    ##
    if(exists(paste("gam.xpose.", pars, ".",
                    object@Runno, sep = ""),
              where = 1)
       & !overwrite) {
      if(ask.for.input){
        cat("\nThere is already a gam object associated with the current\n")
        cat("run number and parameter. It will be overwritten if you proceed.\n")
        cat("Proceed? n(y): ")
        ans <- readline()
        cat("\n")
        if(ans != "y") return()
      } else {
        cat("\nThere is already a gam object associated with the current\n")
        cat("run number and parameter. It will NOT be overwritten.\n")
        return()
      }
    }
    
    
    ##
    ## Run the GAM
    ##
    gamobj1 <- xpose.gam(object,parnam=pars,covnams=covnams,wts.col=wts.col,...)
    
    ## add things to gam object
    gamobj1$pars <- pars
    gamobj1$runno <- object@Runno
    
    
    ##
    ## Save the gam object
    ##
    c1 <- call("assign",pos = 1, paste("gam.xpose.", pars, ".", object@Runno,
             sep= ""), gamobj1, immediate = T)
    eval(c1)
    
    if(exists("current.gam",where=1)){
      remove(pos=1,"current.gam")
    }
    c2 <- call("assign",pos = 1, "current.gam", gamobj1,immediate=T)
    eval(c2)
    
    ##
    ## Return
    ##
    return(invisible(gamobj1))
  }


