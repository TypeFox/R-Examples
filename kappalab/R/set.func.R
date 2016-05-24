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

## Class set.func extends superclass.set.func

##############################################################################

## Constructor from numeric
set.func <- function(object) {
    
    sf <- set.func.internal(object)
    new("set.func", data = sf$data, subsets = sf$subsets, n = sf$n)
}


## Constructor from Mobius.set.func
setMethod("zeta", signature(object = "Mobius.set.func"),
          function(object, ...) {
              
              sf <- zeta.internal(object)
              new("set.func", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from Mobius.card.set.func
setMethod("as.set.func", signature(object = "Mobius.set.func"),
          function(object, ...) {	
              
              set.func(c(object@data, rep(0,2^object@n
                                          - binom.sum(object@n,object@k))))
          }
          )

## Constructor from card.set.func
setMethod("as.set.func", signature(object = "card.set.func"),
          function(object, ...) {
              
              sf <- as.set.func.internal(object)
              new("set.func", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )

## Constructor from Mobius.card.set.func
setMethod("as.set.func", signature(object = "Mobius.card.set.func"),
          function(object, ...) {
              
              sf <- as.set.func.internal(object)
              new("set.func", data = sf$data, subsets = sf$subsets, n = sf$n)
          }
          )


##############################################################################

## Tests the monotonicity of a set function
setMethod("is.monotone", signature(object = "set.func"),
          function(object, verbose = FALSE, epsilon = 1e-9, ...) {	

              if (!is.logical(verbose)) 
                  stop("wrong arguments")
              
              .C("is_monotone_setfunction", 
                 as.integer(object@n), 
                 as.double(object@data), 
                 as.integer(verbose),
                 as.double(epsilon),
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0
          }
          )

## Tests whether the set function is cardinal
setMethod("is.cardinal", signature(object = "set.func"),
          function(object, ...) {
              
              .C("is_cardinal_setfunction", 
                 as.integer(object@n),
                 as.double(object@data),  
                 as.integer(object@subsets),
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0
          }
          )

## Tests whether the set function is k-additive
setMethod("is.kadditive", signature(object = "set.func",k = "numeric"),
          function(object, k, epsilon = 1e-9, ...) {
              
              if (!(k %in% 1:object@n))
                  stop("wrong arguments")

              .C("is_kadditive_setfunction", 
                 as.integer(object@n),
                 as.integer(k),
                 as.double(object@data),
                 as.integer(object@subsets),
                 as.double(epsilon),  
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0
              
          }
          )

##############################################################################

## Show method for object set.func
setMethod("show", signature(object = "set.func"),	
          function(object) {

              cat(paste("\t\t",is(object)[1],"\n",sep=""))
              .C("Rprint_setfunction", 
                 as.integer(object@n),
                 as.integer(object@n),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 as.integer(0),
                 PACKAGE="kappalab")
              cat("")
          }
          )

## Returns a data.frame representing the set function
setMethod("to.data.frame", signature(object = "set.func"),
          function(object, natural=TRUE, ...) {
              
              if (!is.logical(natural))
                  stop("wrong arguments")	
              
              if (natural) {
                  
                  mu.nat <- .C("binary2natural", 
                               as.integer(object@n), 
                               as.double(object@data),
                               as.integer(object@subsets),
                               mu = double(2^object@n),
                               PACKAGE="kappalab")$mu
                  
                  if (object@n < max.n.display) {
                      
                      subsets.char <- .C("k_power_set_char", 
                                         as.integer(object@n),
                                         as.integer(object@n), 
                                         as.integer(object@subsets),
                                         subsets = character(2^object@n),
                                         PACKAGE="kappalab")$subsets
                      
                      d <- data.frame(mu.nat, 
                                  row.names = subsets.char)
                  }
                  else
                      d <- data.frame(mu.nat,
                                      row.names = object@subsets)
              }
              else {
                  if (object@n < max.n.display) {
                      
                      subsets.char <- .C("power_set_binary_char", 
                                         as.integer(object@n), 
                                         subsets = character(2^object@n),
                                         PACKAGE="kappalab")$subsets
                      
                      d <- data.frame(object@data, 
                                      row.names = subsets.char)
                  }
              else
                  d <- data.frame(object@data)
                                }
              names(d) <- is(object)[1]
              d
          }
          )

##############################################################################

## Computes the Shapley value of a set function
setMethod("Shapley.value", signature(object = "set.func"),
          function(object, ...) {

              result <- .C("Shapley_value_setfunction", 
                           as.integer(object@n), 
                           as.double(object@data), 
                           phi = double(object@n),
                           PACKAGE="kappalab")$phi	
              
              names(result) <- 1:object@n
              result
          }
          )

## Computes the Shapley interaction indices of a set function
setMethod("interaction.indices", signature(object = "set.func"),
          function(object, ...) {

              phi <- .C("interaction_indices_setfunction", 
                        as.integer(object@n), 
                        as.double(object@data), 
                        result = double(object@n*object@n),
                        PACKAGE="kappalab")$result
              
              result <- matrix(phi,object@n,object@n)
              diag(result) <- rep(NA,object@n)
              dimnames(result) <- list(1:object@n,1:object@n)
              result
          }
          )

##############################################################################
