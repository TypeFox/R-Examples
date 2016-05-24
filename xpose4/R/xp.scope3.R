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

"xp.scope3"<-
  function(object,
           covnam=xvardef("covariates", object),
           nmods = 3,
           smoother1 = 0,
           arg1 = NULL,
           smoother2 = 1,
           arg2 = NULL,
           smoother3 = "ns",
           arg3 = "df=2",
           smoother4 = "ns",
           arg4 = "df=3", 
	   excl1 = NULL,
           excl2 = NULL,
           excl3 = NULL,
           excl4 = NULL,
           extra = NULL,
           subset=xsubset(object),
           ...)
{
  ##Get data
  data <- Data(object,subset=subset)
  if(any(is.null(data))) return("The subset expression is invalid.")    
  
  ## parnam <- xvardef("parms", object)
  ## covnams <- xvardef("covariates", object)
  step.list <- as.list(covnam)
  names(step.list) <- covnam

  if(is.null(extra)) extra <- list()	## Set up the smoother lists

  for(cov in covnam) {
    ## Check that the covariate is not excluded completely
    if(!is.na(match(cov, excl1)) && !is.na(match(cov, excl2)) && !
       is.na(match(cov, excl3)) && !is.na(match(cov, excl4))) 
      {
	stop("A covariate cannot be excluded from all levels in the GAM scope.\n")
      }

    ## check that categorical covariates have more than one factor
    if(is.factor(data[, cov]) && nlevels(data[,cov])==1){
      step.list=step.list[names(step.list)!=cov]
    } else
    { # create the scope
      tmp.scope <- character(0)
      for(i in 1:nmods) {
        excl <- eval(parse(text = paste("excl", i, sep = "")))
        if(!is.na(match(cov, excl)))
          next
        if(is.factor(data[, cov]) && i == 1) {
          tmp.scope <- c(tmp.scope, "1")
          next
        }	else if(is.factor(data[, cov]) && i == 2) {
          tmp.scope <- c(tmp.scope, cov)
          next
        }  else if(is.factor(data[, cov]) && i > 2) {
          next
        }
        ## Check if we have any specific covariate settings for smoother
        if(!is.null(extra[[cov]][[paste("sm", i, sep = "")]])) 
          {
            sm <- extra[[cov]][[paste("sm", i, sep = "")]]
          } else {
            sm <- eval(parse(text = paste("smoother", i, 
                               sep = "")))
          }
        if(sm == 0) {
          tmp.scope <- c(tmp.scope, "1")
        }	else if(sm == 1) {
          tmp.scope <- c(tmp.scope, cov)
        }  else {
          ## Check if we have any specific settings for  arg
          if(!is.null(extra[[cov]][[paste("arg", i, sep
                                          = "")]])) {
            if(extra[[cov]][[paste("arg", i, sep = "")]] == 
               "") {
              arg <- NULL
            } else {
              arg <- extra[[cov]][[paste("arg", i, sep = 
                                         "")]]
            }
          } else {
            arg <- eval(parse(text = paste("arg", i, sep
                                = "")))
          }
          if(is.na(pmatch("I((",sm)))
            tmp.scope <- c(tmp.scope, paste(sm, "(", cov, 
                                            if(is.null(arg)) ")" else paste(",", arg, ")", 
                                                                            sep = "")))
          else
            tmp.scope <- c(tmp.scope,sm)
          
        }
      }
      tmp.scope <- eval(parse(text = paste("~", paste(tmp.scope, collapse = "+"
                                ))))
      step.list[[cov]] <- tmp.scope
    }
  }
  step.list
}
