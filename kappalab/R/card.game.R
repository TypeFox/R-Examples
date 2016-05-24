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
################################################################################

## Class card.game extends card.set.func

##############################################################################

## Constructor from numeric
card.game <- function(object) {
    
    if (!is.numeric(object))
        stop("wrong argument")

    new("card.game", data = object, n = length(object) - 1)
}

## Constructor from card.set.func
setMethod("as.card.game", signature(object = "card.set.func"),
          function(object,...) {
              
              new("card.game", data = object@data, n = object@n)
          }
          )

## Constructor from game
setMethod("as.card.game", signature(object = "game"),
          function(object, ...) {
              
              csf <- as.card.set.func.internal(object)
              new("card.game", data = csf$data, n = csf$n)
          }
          )

## Constructor from card.set.func
## The conjugate of a set function is a game
setMethod("conjugate", signature(object = "card.set.func"),
          function(object, ...) {
              
              card.game(object@data[object@n+1] - object@data[(object@n+1):1])
          }
          )

##############################################################################

## Computes the Choquet integral of f w.r.t mu (f >= 0)
setMethod("Choquet.integral",signature(object = "card.game", f = "numeric"),
          function(object,f,...) {

              if (!(object@n == length(f)))
                  stop("wrong arguments")
              
              sum(sort(f)
                  * (object@data[(object@n+1):2]-object@data[object@n:1]))
          }
          )

## Computes the Sugeno integral of f w.r.t mu (f >= 0)
setMethod("Sugeno.integral",signature(object = "card.game", f = "numeric"),
          function(object,f,...) {
              
              if (!(is.positive(f) && object@n == length(f)))
                  stop("wrong arguments")
              
              max(pmin(sort(f),object@data[(object@n+1):2]))
          }
          )

## Computes the Sipos integral of f w.r.t mu 
setMethod("Sipos.integral",signature(object = "card.game", f = "numeric"),
          function(object,f,...) {
              
              if (!(is.double(f) && object@n == length(f)))
                  stop("wrong arguments")

              return(Choquet.integral(object,pmax(f,0))
                     - Choquet.integral(object,pmax(-f,0)))
          }
          )

##############################################################################
