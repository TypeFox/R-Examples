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

change.cat.cont <-
  function(object,
           listall=TRUE,
           classic=FALSE,
           to.cat.vec=NULL,
           to.cont.vec=NULL,
           change.type.vec=NULL,
           ...){

    if(!is.null(change.type.vec)){
      for(j in change.type.vec){
        
        if(is.na(pmatch(j, names(object@Data)))) {
          cat("Couldn't find",j,"in the current database. Skipping!\n")
          next
        }

        if(is.factor(object@Data[,j])) {
          to.cont.vec <- c(to.cont.vec,j)
        } else {
          to.cat.vec <- c(to.cat.vec,j)
        }
      }
    }
    
    if(is.null(to.cont.vec) & is.null(to.cat.vec) & is.null(change.type.vec)){
      if(listall) {
        db.names(object)
      }
    }
    
    if(is.null(to.cont.vec) & is.null(to.cat.vec) & is.null(change.type.vec)){
      cat("\nPlease type the names of any CATEGORICAL variables you\n")
      cat("want to change into continuous, one per line, and end with a\n")
      cat("blank line:\n")
      cats <- scan(what=character())
    } else {
      cats <- to.cont.vec
    }

    if(length(cats)){
      for(i in cats) {
        
        if(is.na(pmatch(i, names(object@Data)))) {
          cat("Couldn't find",i,"in the current database. Skipping!\n")
          next
        }
        
        if(is.factor(object@Data[,i])){
          cat("\n  Transforming",i,"from categorical to continuous\n",sep=" ")
          object@Data[,i] <- as.numeric(levels(object@Data[,i]))[object@Data[,i]]
        } else {
          cat("  ",i," is already continuous, no changes made", sep="")
        }
      }
    }


    if(is.null(to.cont.vec) & is.null(to.cat.vec) & is.null(change.type.vec)){
      cat("Please type the names of any CONTINUOUS variables you\n")
      cat("want to change into categorical, one per line, and end with a\n")
      cat("blank line:\n")
      conts <- scan(what=character())
    } else {
      conts <- to.cat.vec
    }
    
    if(length(conts)){
      for(i in conts) {
        
        if(is.na(pmatch(i, names(object@Data)))) {
          cat("Couldn't find",i,"in current database. Skipping!\n")
          next
        }
        
        if(!is.factor(object@Data[,i])){
          cat("\n  Transforming",i,"from continuous to categorical\n",sep=" ")
          object@Data[,i] <- as.factor(object@Data[,i])
        } else {
          cat("  ",i," is already categorical, no changes made\n", sep="")
        }
      }
    }
    
    if (classic==TRUE) {
      c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
      eval(c1)
      c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
      eval(c2)
    } 
    return(object)
  }

"change.cat.cont<-" <-
  function(object,
           listall=TRUE,
           classic=FALSE,
           to.cat.vec=NULL,
           to.cont.vec=NULL,
           ...,
           value){
    object <- change.cat.cont(object,listall=listall,calssic=classic,
                              to.cat.vec=to.cat.vec,
                              to.cont.vec=to.cont.vec,
                              change.type.vec=value,
                              ...)
    return(object)
  }
