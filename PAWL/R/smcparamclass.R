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
setClass("smcparameters",
         representation(nparticles = "numeric",  
                        temperatures = "numeric",
                        nmoves = "numeric",
                        ESSthreshold = "numeric",
                        movetype = "character",
                        movescale = "numeric",
                        resamplingscheme = "character"))

setGeneric("smcparameters", function(...)standardGeneric("smcparameters"))
smcparameters.constructor <- function(..., nparticles, temperatures, nmoves, ESSthreshold,
                                      movetype, movescale, resamplingscheme){
    if (missing(nparticles))
        nparticles <- 1
    if (missing(temperatures))
        temperatures <- seq(from = 0.01, to = 1, length.out = 100)
    if (missing(nmoves))
        nmoves <- 1
    if (missing(ESSthreshold))
        ESSthreshold <- 0.5
    if (missing(movetype)){
        movetype <- "independent"
    }
    if (movetype != "independent" && movetype != "randomwalk"){
        stop(sQuote("movetype"), " should be either 'independent' or 'randomwalk'!")
    }
    if (missing(movescale)){
        movescale <- 0.1
    }
    if (missing(resamplingscheme))
        resamplingscheme <- "systematic"
    if (resamplingscheme != "multinomial" && 
        resamplingscheme != "residual" && resamplingscheme != "systematic"){
        stop(sQuote("resamplingscheme"), 
             " should be either 'multinomial', 'residual' or 'systematic'!")
    }
    new("smcparameters", 
        nparticles = nparticles, temperatures = temperatures, nmoves = nmoves,
        ESSthreshold = ESSthreshold, movetype = movetype, movescale = movescale, 
        resamplingscheme = resamplingscheme)
}
setMethod("smcparameters",
          definition = function(..., nparticles, temperatures, nmoves, ESSthreshold, movetype, movescale,
                                resamplingscheme){
              smcparameters.constructor(
                                        nparticles = nparticles, temperatures = temperatures,
                                        nmoves = nmoves, ESSthreshold = ESSthreshold, movetype = movetype,
                                        movescale = movescale, resamplingscheme = resamplingscheme)
          })

setMethod(f = "show", signature = "smcparameters", 
          def = function(object){
            cat("Object of class ", class(object), ".\n", sep = "")
            cat("*number of particles:", object@nparticles, "\n")
            cat("*number of distributions:", length(object@temperatures), "\n")
            cat("*temperatures: ...", tail(object@temperatures), "\n")
            cat("*number of moves:", object@nmoves, "\n")
            cat("*ESS threshold:", object@ESSthreshold, "\n")
            cat("*move type:", object@movetype, "\n")
            if (object@movetype == "randomwalk") cat("*move scale:", object@movescale, "\n")
            cat("*resampling scheme:", object@resamplingscheme, "\n")
          })


