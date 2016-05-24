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

"change.xvardef"<-
  function(object,
           var=".ask", # a one element vector
           def=".ask", # a vector (can have mutiple elements)
           listall=TRUE,
           classic=FALSE,
           check.var=FALSE,
           ...){


    if(is.null(def)|is.null(var)){
      if(listall) {
        db.names(object)
      }
    }

    if(var==".ask"){
      tmp.txt <- paste("\nPlease type the name of the xpose",
                       "variable you would like to change or add.",
                       "Current variables are (press return to exit without",
                       "changing):\n")
      cat(unlist(strsplit(tmp.txt," ")),fill=60,file="")
      ##xpose.string.print(tmp.txt)
      cat(names(object@Prefs@Xvardef),"\n",fill=60)
      var <- readline()
    }
    if(var==""){
      if (classic==TRUE) {
        return(cat(""))
      } else {
        return(object)
      }
    } else {
      if(classic | check.var){
        tmp <- is.na(pmatch(var, names(object@Prefs@Xvardef)))
        if(tmp) {
          cat("Adding", var, "to the current database!\n")
          ##return(Recall(object, listall = F, var=NULL,def=def,classic=classic,...))
        }
      }
    }

    if(!is.null(def)){
      if(length(def)== 1 && def==".ask"){
        tmp.txt <- paste("\nThe variable",
                         paste("\"", var, "\"",sep=""),
                         "is currently associated with the following observed items:\n")
        cat(unlist(strsplit(tmp.txt," ")),fill=60,file="")
        ##xpose.string.print(tmp.txt)
        cat(object@Prefs@Xvardef[[var]],fill=60)

        tmp.txt <-
          paste("\nPlease type the name(s) of the observed items in the database",
                "(object@Data) you would like to associate with the xpose",
                "variable",paste("\"", var, "\"",sep=""),
                "and finish with a blank line.",
                "Possible values are",
                "(press return to exit without changing, NULL to remove variable):\n")

        cat(unlist(strsplit(tmp.txt," ")),fill=60,file="")
        ##xpose.string.print(tmp.txt)
        cat(names(object@Data),"\n",fill=60)
        def <- scan(what=character())
      }
      if(length(def)==0){
        if (classic==TRUE) {
          return(cat(""))
        } else {
          return(object)
        }
      }
    }

    if(is.null(def) || def == "NULL" || def == "null") {
      object@Prefs@Xvardef[[var]] <- NULL
    } else {
      tmp <- is.na(match(def, names(object@Data)))
      if(any(tmp)) {
        cat("Couldn't find", def[tmp], "in the current database!\n")
        return(Recall(object,
                      listall = F,
                      var=var,
                      def=".ask",
                      classic=classic,
                      check.var=check.var,
                      ...))
      }
      object@Prefs@Xvardef[[var]] <- def
    }

    if (classic==TRUE) {
      c1<-call("assign",paste("xpdb", object@Runno, sep = ""), data, immediate=T, envir = .GlobalEnv)
      eval(c1)
      c2<-call("assign",pos = 1, ".cur.db", eval(as.name(paste("xpdb", object@Runno, sep = ""))))
      eval(c2)
      return(cat(""))
      
    } else {
      return(object)
    }
  }

"change.xvardef<-" <-
  function(object,
           var, # a one element vector
           listall=FALSE,
           classic=FALSE,
           check.var=FALSE,
           ...,
           value){

    object <- change.xvardef(object,
                             var=var,
                             def=value,
                             listall=listall,
                             classic=classic,
                             check.var=check.var,
                             ...)
    return(object)
  }
