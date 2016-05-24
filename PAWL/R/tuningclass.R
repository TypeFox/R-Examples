###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################
setClass("tuningparameters",
         representation(nchains = "numeric", niterations = "numeric", 
                        computemean = "logical", computemeanburnin = "numeric",
                        saveeverynth = "numeric"))

setGeneric("tuningparameters", function(...)standardGeneric("tuningparameters"))
tuningparameters.constructor <- function(..., nchains, niterations, storeall, 
                                         computemean, computemeanburnin, saveeverynth){
    if (missing(nchains))
        nchains <- 1
    if (missing(niterations))
        stop(sQuote("niterations"), "has to be specified")
    if (missing(saveeverynth)){
        if (missing(storeall)){
            #cat("storeall unspecified: set to FALSE\n")
            storeall <- FALSE
            saveeverynth <- -1
        } else {
            if (storeall){
                saveeverynth <- 1
            } else {
                saveeverynth <- -1
            }
        }
    }
    if (missing(computemean)){
      computemean <- FALSE
      #cat("computemean unspecified: set to FALSE\n")
    }
    if (missing(computemeanburnin)){
        computemeanburnin <- 0
    }
    new("tuningparameters", 
        nchains = nchains, niterations = niterations, 
        computemean = computemean, computemeanburnin = computemeanburnin, 
        saveeverynth = saveeverynth)
}
setMethod("tuningparameters",
          definition = function(..., nchains, niterations, storeall, computemean, 
                                computemeanburnin, saveeverynth){
              tuningparameters.constructor(
                               nchains = nchains, niterations = niterations, 
                               storeall = storeall, computemean = computemean,
                               computemeanburnin = computemeanburnin,
                               saveeverynth = saveeverynth)
          })

setMethod(f = "show", signature = "tuningparameters", 
          def = function(object){
              cat("Object of class ", class(object), ".\n", sep = "")
              cat("*number of parallel chains:", object@nchains, "\n")
              cat("*number of iterations:", object@niterations, "\n")
              cat("*compute mean:", object@computemean, "\n")
              cat("*compute mean (burnin):", object@computemeanburnin, "\n")
              cat("*save every nth iteration:", object@saveeverynth, "\n")
          })


