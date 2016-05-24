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

DSC_CluStream <- function(
  m=100,
  horizon=1000, 
  t=2,
  k=NULL
) {
  
  ### Java code does parameter checking
  paramList <- list(
    h = as.integer(horizon),
    m = as.integer(m),
    t = t
    )
  
  # converting the param list to a cli string to use in java
  cliParams <- convert_params(paramList)
  
  # initializing the clusterer
#  clusterer <- .jnew("moa/clusterers/clustream/Clustream")
  clusterer <- .jcast(.jnew("moa/clusterers/clustream/WithKmeans"),
    "moa/clusterers/AbstractClusterer")
  options <- .jcall(clusterer, "Lmoa/options/Options;", "getOptions")
  .jcall(options, "V", "setViaCLIString", cliParams)
  .jcall(clusterer, "V", "prepareForUse")
  
  
  # initializing the R object
  clus <- structure(
    list(
      description = "CluStream",
      options = cliParams,
      javaObj = clusterer
    ),
    class = c("DSC_CluStream","DSC_Micro","DSC_MOA","DSC")
  )

  if(!is.null(k)) clus <- DSC_TwoStage(clus, 
    DSC_Kmeans(k=k, weighted=TRUE, nstart=5)) 
  
  clus
}
