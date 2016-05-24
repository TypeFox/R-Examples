##############################################################################
#
# Copyright © 2005 Michel Grabisch and Ivan Kojadinovic    
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

## Class card.capacity extends card.game

##############################################################################

## Constructor from numeric
card.capacity <- function(object) {
    
    if (!is.numeric(object))
        stop("wrong argument")

    new("card.capacity", data = object, n = length(object) - 1)
}

## Constructor from card.set.func
setMethod("as.card.capacity", signature(object = "card.set.func"),
          function(object,...) {
              
              new("card.capacity", data = object@data, n = object@n)
          }
          )

## Constructor from capacity
setMethod("as.card.capacity", signature(object = "capacity"),
          function(object, ...) {
              
              csf <- as.card.set.func.internal(object)
              new("card.capacity", data = csf$data, n = csf$n)
          }
          )

## Constructor from card.capacity
## The conjugate of a capacity is a capacity
setMethod("conjugate", signature(object = "card.capacity"),
          function(object, ...) {
              
              card.capacity(object@data[object@n+1]
                            - object@data[(object@n+1):1])
          }
          )

## Constructors for simple capacities
lower.capacity <- function(n) {
    
    new("card.capacity", data = c(rep(0,n),1), n = n)
}

upper.capacity <- function(n) {
    
    new("card.capacity", data = c(0,rep(1,n)), n = n)
}

uniform.capacity <- function(n) {

    new("card.capacity", data = 0:n/n, n = n)
}

##############################################################################

## Tests whether the set function is normalized
setMethod("is.normalized", signature(object = "card.capacity"),
          function(object,...) {
              
              object@data[object@n+1] == 1
          }
          )

## Returns the normalized version of the capacity
setMethod("normalize", signature(object = "card.capacity"),
          function(object,...) {
              
              card.capacity(object@data / object@data[object@n+1])
          }
          )

##############################################################################

## Computes the orness of a capacity (normalized)
setMethod("orness", signature(object = "card.capacity"),
          function(object, ...) {
              return(sum(object@data[1:object@n])
                     / ((object@n - 1)
                        * object@data[object@n + 1]))
          }
          )

## Computes the veto indices, here equal to 1 - orness
setMethod("veto", signature(object = "card.capacity"),
          function(object, ...) {
              
              result <- 1 - orness(object)
              names(result) <- "All"
              result	
          }
          )

## Computes the favor indices, here equal to orness
setMethod("favor", signature(object = "card.capacity"),
          function(object, ...) {
              
              result <- orness(object)
              names(result) <- "All"
              result		
          }
          )

## Computes the normalized variance of a capacity
setMethod("variance", signature(object = "card.capacity"),
          function(object, ...) {
              
              (sum(((object@data[2:(object@n+1)] -
                     object@data[1:object@n])/object@data[object@n+1])^2) -
               1/object@n) * object@n / (object@n - 1)	
          }
          )

## Computes the normalized entropy of a capacity
setMethod("entropy", signature(object = "card.capacity"),
          function(object, ...) {	
          
              sum(sapply((object@data[2:(object@n+1)] - object@data[1:object@n])
                         / object@data[object@n+1],shannon.function)) /
                             log(object@n)		
          }
          )

##############################################################################
