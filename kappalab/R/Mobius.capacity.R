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

## Class Mobius.capacity extends Mobius.game

##############################################################################

## Constructor from numeric
Mobius.capacity <- function(object, n, k) {
    
    msf <- Mobius.set.func.internal(object, n, k)
    new("Mobius.capacity", data = msf$data, subsets = msf$subsets, 
        n = msf$n, k=msf$k)
}

## Constructor from numeric
additive.capacity <- function(v) {
    
    if (!is.positive(v))
        stop("wrong argument")

    Mobius.capacity(c(0,v),length(v),1)	
}

## Constructor from Mobius.set.func
setMethod("as.Mobius.capacity", signature(object = "Mobius.set.func"),
          function(object, ...)	{
              
              new("Mobius.capacity", data = object@data, 
                  subsets = object@subsets, n = object@n, k = object@k)
          }
          )

## Constructor from capacity
setMethod("Mobius", signature(object = "capacity"),
          function(object, ...)	{
              
              msf <- k.truncate.Mobius.internal(object, object@n)
              new("Mobius.capacity", data = msf$data, subsets = msf$subsets, 
                  n = msf$n, k=msf$k)
          }
          )

##############################################################################

## Tests whether the set function is normalized
setMethod("is.normalized", signature(object = "Mobius.capacity"),
          function(object,...) {
              
              sum(object@data) == 1
          }
          )

## Returns the normalized version of the capacity
setMethod("normalize", signature(object = "Mobius.capacity"),
          function(object,...)  {
              
              Mobius.capacity(object@data / sum(object@data),
                              object@n, object@k)
          }
          )

##############################################################################

## Computes the veto indices of a capacity
setMethod("veto", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              result <- .C("veto_Mobius", 
                           as.integer(object@n),
                           as.integer(object@k), 
                           as.double(object@data),
                           as.integer(object@subsets), 
                           result = double(object@n),
                           PACKAGE="kappalab")$result
              names(result) <- 1:object@n
              result	
          }
          )

## Computes the favor indices of a capacity
setMethod("favor", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              result <- .C("favor_Mobius", 
                           as.integer(object@n),
                           as.integer(object@k), 
                           as.double(object@data),
                           as.integer(object@subsets), 
                           result = double(object@n),
                           PACKAGE="kappalab")$result

              names(result) <- 1:object@n
              result	
          }
          )

## Computes the orness degree of a capacity
setMethod("orness", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              .C("orness_Mobius", 
                 as.integer(object@n), 
                 as.integer(object@k),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 res = double(1),
                 PACKAGE="kappalab")$res
          }
          )

## Computes the normalized variance of a capacity
setMethod("variance", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              .C("variance_Mobius", 
                 as.integer(object@n), 
                 as.integer(object@k),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 res = double(1),
                 PACKAGE="kappalab")$res
          }
          )

## Computes the normalized entropy of a capacity
setMethod("entropy", signature(object = "Mobius.capacity"),
          function(object, ...) {
              
              .C("entropy_Mobius", 
                 as.integer(object@n), 
                 as.integer(object@k),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 res = double(1),
                 PACKAGE="kappalab")$res
          }
          )

##############################################################################
