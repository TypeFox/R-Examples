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

## Class capacity extends game

##############################################################################

## Constructor from numeric
capacity <- function(object) {
    
    sf <- set.func.internal(object)
    new("capacity", data = sf$data, subsets = sf$subsets, n = sf$n)
}


## Constructor from set.func
setMethod("as.capacity", signature(object = "set.func"),
          function(object,...) {
              
              new("capacity", data = object@data, subsets = object@subsets, 
                  n = object@n)
          }
          )

## Constructor
setMethod("as.capacity", signature(object = "card.capacity"),
          function(object,...) {
              
              sf <- as.set.func.internal(object)
              new("capacity", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from Mobius.game
setMethod("zeta", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              sf <- zeta.internal(object)
              new("capacity", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from capacity
## The conjugate of a capacity is a capacity
setMethod("conjugate", signature(object = "capacity"),
          function(object, ...) {
              
              sf <- conjugate.internal(object)
              new("capacity", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

##############################################################################

## Tests whether the set function is normalized
setMethod("is.normalized", signature(object = "capacity"),
          function(object,...) {
              
              object@data[2^object@n] == 1
          }
          )

## Returns the normalized version of the capacity
setMethod("normalize", signature(object = "capacity"),
          function(object,...) {
              
              mu <- object
              mu@data <- mu@data / mu@data[2^mu@n]
              mu
          }
          )

##############################################################################

## Computes the veto indices of a capacity
setMethod("veto", signature(object = "capacity"),
          function(object, ...) {
              
              result <- .C("veto_capacity", 
                           as.integer(object@n), 
                           as.double(object@data), 
                           result = double(object@n),
                           PACKAGE="kappalab")$result
              
              names(result) <- 1:object@n
              result	
          }
          )

## Computes the favor indices of a capacity
setMethod("favor", signature(object = "capacity"),
          function(object, ...) {

              result <- .C("favor_capacity", 
                           as.integer(object@n), 
                           as.double(object@data), 
                           result = double(object@n),
                           PACKAGE="kappalab")$result
              
              names(result) <- 1:object@n
              result	
          }
          )

## Computes the orness degree of a capacity
setMethod("orness", signature(object = "capacity"),
          function(object, ...) {
              
              .C("orness_capacity", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 res = double(1),
                 PACKAGE="kappalab")$res
          }
          )

## Computes the normalized variance of a capacity
setMethod("variance", signature(object = "capacity"),
          function(object, ...) {
              
              .C("variance_capacity", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 res = double(1),
                 PACKAGE="kappalab")$res
      }
          )

## Computes the normalized entropy of a capacity
setMethod("entropy", signature(object = "capacity"),
          function(object, ...) {
              .C("entropy_capacity", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 res = double(1),
                 PACKAGE="kappalab")$res
          }
          )

##############################################################################
