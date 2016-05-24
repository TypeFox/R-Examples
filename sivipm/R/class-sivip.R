###################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
#--------------------------------------------
# Class sivip
# SLOTS:
# Class of the objects returned by the function 'sivipm'
# METHODS:
# print, show, summary
# CREATOR:
# the function 'sivipm'
# Exported
#--------------------------------------------
sivip <- setClass("sivip",
                  slots=c(
                      fo.isivip="vector", # ncol(X)
                      fo.isivip.percent="vector", # ncol(X)
                      tsivip="vector", # ncol(X)
                      tsivip.percent="vector", # ncol(X)
                      monosignif="vector", # nmono
                      correlalea="matrix",  # ncol(Y) X ncol(Y)
                      simca.signifcomponents="matrix", #  ncompo X ncol(Y) +1
                       lazraq.signifcomponents="vector", # ncompo
                       output="list"),
                  # default value
                 prototype = list(
                   fo.isivip=NULL,
                   fo.isivip.percent=NULL,
                   tsivip=NULL,
                   tsivip.percent=NULL,
                   monosignif=NULL,
                   correlalea=NULL,
                   simca.signifcomponents=NULL,
                   simca.signif2components=NULL,
                   lazraq.signifcomponents=NULL,
                   output=NULL
                  ))

#--------------------------------------------
# Method 'print':
#--------------------------------------------
print.sivip <- function (x, all=FALSE, ...) {
    # all=T : print of all the components.
  # Otherwise, the output slot is not written
  for (nom in slotNames(x)) {
    if (!is.null(slot(x,nom)) &&  (nom != "output")) {
        if ((nom =="fo.isivip") &&  (all==FALSE)) {
            next
        }
        if ((nom =="tsivip") &&  (all==FALSE)) {
            next
        }
      cat( nom, "\n")
      print(slot(x,nom), ...)
      cat("\n")
    }
  }
  
  if (!is.null(x@output)) {
    if (all==TRUE) {
      cat("Component output\n")
      print(x@output, ...)
      cat("\n")
    } else {
      warning("The component 'output' is only written when the option 'all' of the 'print' function is TRUE\n")
    }
  }
} # fin print.sivip

#--------------------------------------------
# Method 'show':
#--------------------------------------------
 setMethod("show",  signature(object="sivip"),
           function(object) {
           print(object, all=TRUE)
         })

                       
#--------------------------------------------
# Method 'summary':
#--------------------------------------------
summary.sivip <-  function(object, ...) {
           print(object, all=FALSE)
           return(invisible())
         }

setMethod("summary",  signature(object="sivip"),
           definition = summary.sivip)

#--------------------------------------------
# Method 'getNames'
# names of the non null components
#--------------------------------------------
setGeneric("getNames",
           function(x){
               value <- standardGeneric("getNames")
           })


                   
 setMethod("getNames", signature(x="sivip"),
           function(x) {
               cat("*** Names of the slots:\n")
               for (nom in slotNames(x)) {
                   if (length(slot(x, nom)) >0) {
                       cat(nom, " " )
                   }
               }
               cat("\n")


               output <- c()
               for (nom in names(x@output)) {
                   if (!is.null(x@output[[nom]]))
                       output <- c(output, nom)
               }
               if (length(output) > 0) {
                   cat("*** The slot 'output' is a list with components:\n")
                   print(output)
                   if (length(x@output$PLS) >0 ) {
                       cat("*** The component 'output$PLS' is a list  with components:\n")
                       print(names(x@output$PLS))
                   }
               }

            invisible()       
           }
           )
