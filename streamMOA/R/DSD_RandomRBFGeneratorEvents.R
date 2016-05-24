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


# RandomRBFGeneratorEvents options:
# -m modelRandomSeed
# -i instanceRandomSeed
# -K numCluster
# -k numClusterRange
# -R kernelRadius
# -r kernelRadiusRange
# -d densityRange
# -V speed
# -v speedRange
# -N noiseLevel
# -E eventFrequency
# -M eventMergeWeight
# -P eventSplitWeight
# -a numAtts (dimensionality)
# because there are so many parameters, let's only use a few key ones...
DSD_RandomRBFGeneratorEvents <- function(k=3, d=2, 
  numClusterRange=3L,
  kernelRadius=0.07,
  kernelRadiusRange=0,
  densityRange=0,
  speed=100L,
  speedRange=0L,
  noiseLevel=0.1,
  noiseInCluster=FALSE,
  eventFrequency=30000L,
  eventMergeSplitOption=FALSE,
  eventDeleteCreate=FALSE,
  modelSeed=NULL, 
  instanceSeed=NULL
) {
  #TODO: need error checking on the params
  
  if(is.null(modelSeed)) modelSeed <- as.integer(
    runif(1L, 0, .Machine$integer.max))
  if(is.null(instanceSeed)) instanceSeed <- as.integer(
    runif(1L, 0, .Machine$integer.max))
  
  paramList <- list(
    m=modelSeed,
    i=instanceSeed,
    K=k,
    k=as.integer(numClusterRange),
    R=kernelRadius,
    r=kernelRadiusRange,
    d=densityRange,
    V=speed,
    v=speedRange,
    N=noiseLevel,
    E=eventFrequency,
    n=noiseInCluster,
    M=eventMergeSplitOption,
    C=eventDeleteCreate,
    a=d
  )
  
  # converting the param list to a cli string to use in java
  cliParams <- convert_params(paramList)
  
  # initializing the clusterer
  strm <- .jnew("moa/streams/clustering/RandomRBFGeneratorEvents")
  options <- .jcall(strm, "Lmoa/options/Options;", "getOptions")
  .jcall(options, "V", "setViaCLIString", cliParams)
  .jcall(strm, "V", "prepareForUse")
  
  l <- list(
    description = "Random RBF Generator Events (MOA)",
    k = k,
    d = d,
    cliParams = cliParams,
    javaObj = strm
    )
  
  class(l) <- c("DSD_RandomRBFGeneratorEvents","DSD_MOA","DSD")
  l
}
