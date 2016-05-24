####################################################################################
####				Generic functions				####
####################################################################################

# Generic functions for S3 methods
# Warning: the arguments in function must match exactly !
showall <- function(object) UseMethod("showall")
setGeneric(name="showall", def=function(object){standardGeneric("showall")})

gof.test <- function(object,which="npde",parametric=TRUE, ...) UseMethod("gof.test")
setGeneric(name="gof.test", def=function(object,which="npde",parametric=TRUE, ...){standardGeneric("gof.test")})

set.plotoptions <- function(object, ...) UseMethod("set.plotoptions")
setGeneric(name="set.plotoptions", def=function(object,...){standardGeneric("set.plotoptions")})

# Generic functions for S4 methods
#setGeneric(name="test",def=function(object,which="npde",parametric=TRUE, na.action="na.omit", ...){standardGeneric("test")})

# Functions specific to NpdeData objects
setGeneric(name="read.npdeData",def=function(object,header,sep,na.strings,detect,verbose) {standardGeneric("read.npdeData")}
)

setGeneric(name="read.npdeSimData",def=function(object,header,sep,na.strings, verbose){standardGeneric("read.npdeSimData")}
)

# Main algorithm functions
npde.main <- function(object) UseMethod("npde.main")
setGeneric(name="npde.main", def=function(object){standardGeneric("npde.main")})

setGeneric(name="npde.save",
  def=function(object, ...){standardGeneric("npde.save")}
)
setGeneric(name="npde.graphs",
  def=function(object,...){standardGeneric("npde.graphs")}
)

####################################################################################
####				Check if mclust installed			####
####################################################################################
# 
# .npde.mclust<-require("mclust",quietly=TRUE)
# if(.npde.mclust) {
# 	cat("Loading mclust\n")
# 	library(mclust)
# }

####################################################################################
