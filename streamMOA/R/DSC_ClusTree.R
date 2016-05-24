#######################################################################
# stream -  Infrastructure for Data Stream Mining
# Copyright (C) 2013 Michael Hahsler, Matthew Bolanos, John Forrest 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


# ClusTree (anytime clustering) options:
# IntOption("horizon", 'h', "Range of the window.", 1000)
# IntOption("maxHeight", 'H', "The maximal height of the tree", 8)

#Reclustering: suggests EM or k-means

DSC_ClusTree <- function(horizon=1000, maxHeight=8, lambda=NULL) {

  ### Java code does parameter checking
  paramList <- list(
    h=horizon,
    H=maxHeight
  )

  # converting the param list to a cli string to use in java
  cliParams <- convert_params(paramList)
  
  # initializing the clusterer
  clusterer <- .jcast(.jnew("moa/clusterers/clustree/ClusTree"),
    "moa/clusterers/AbstractClusterer")
  options <- .jcall(clusterer, "Lmoa/options/Options;", "getOptions")
  .jcall(options, "V", "setViaCLIString", cliParams)
  .jcall(clusterer, "V", "prepareForUse")

  # overwrite lambda
  if(!is.null(lambda))
  	.jfield(clusterer,"negLambda") <- -1*lambda

  # initializing the R object
  structure(
    list(
      description = "ClusTree",
      options = cliParams,
      javaObj = clusterer
    ),
    class = c("DSC_ClusTree","DSC_Micro","DSC_MOA","DSC")
  )
}
