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
##################################################################################

## Class Mobius.set.func extends superclass.set.func

##############################################################################


## Constructor from numeric
Mobius.set.func <- function(object, n, k) {
    
    msf <- Mobius.set.func.internal(object, n, k)
    new("Mobius.set.func", data = msf$data, subsets = msf$subsets, 
        n = msf$n, k=msf$k)
}

## Constructor from set.func
setMethod("k.truncate.Mobius", signature(object = "set.func", k = "numeric"),
          function(object, k, ...) {	

              msf <- k.truncate.Mobius.internal(object, k)
              new("Mobius.set.func", data = msf$data, subsets = msf$subsets, 
                  n = msf$n, k=msf$k)
          }
          )

## Constructor from Mobius.set.func
setMethod("k.truncate.Mobius", signature(object = "Mobius.set.func",
                                         k = "numeric"),
          function(object, k, ...) {
              
              if (!(k %in% 1:object@n))
                  stop("wrong arguments")
              
              bs <- binom.sum(object@n, k)
              
              new("Mobius.set.func", data = object@data[1:bs], n = object@n, 
                  subsets = object@subsets[1:bs], k = k)
          }
          )

## Constructor from set.func
setMethod("Mobius", signature(object = "set.func"),
          function(object, ...)	 {
              
              msf <- k.truncate.Mobius.internal(object, object@n)
              new("Mobius.set.func", data = msf$data, subsets = msf$subsets, 
                  n = msf$n, k=msf$k)
          }
          )

## Constructor from set.func
setMethod("as.Mobius.set.func", signature(object = "set.func"),
          function(object, ...) {
              
              mu <- .C("binary2natural", 
                       as.integer(object@n),
                       as.double(object@data),
                       as.integer(object@subsets),
                       mu = double(2^object@n), 
                       PACKAGE="kappalab")$mu

              new("Mobius.set.func", data = mu, subsets = object@subsets, 
                  n = object@n, k = object@n)
          }
          )

## Constructor from card.set.func
setMethod("as.Mobius.set.func", signature(object = "card.set.func"),
          function(object, ...) {
              
              subsets <-  .C("k_power_set", 
                             as.integer(object@n),
                             as.integer(object@n),
                             subsets = integer(2^object@n),
                             PACKAGE="kappalab")$subsets
              
              mu <- .C("cardinal2setfunction", 
                       as.integer(object@n),
                       as.double(object@data),
                       mu = double(2^object@n), 
                       PACKAGE="kappalab")$mu 
              
              mu <- .C("binary2natural", 
                       as.integer(object@n),
                       as.double(mu),
                       as.integer(subsets),
                       mu = double(2^object@n), 
                       PACKAGE="kappalab")$mu

              new("Mobius.set.func", data = mu, subsets = subsets, 
                  n = object@n, k = object@n)
          }
          )

## Constructor from Mobius.card.set.func
setMethod("as.Mobius.set.func", signature(object = "Mobius.card.set.func"),
          function(object, ...) {
              
              subsets <-  .C("k_power_set", 
                             as.integer(object@n),
                             as.integer(object@n),
                             subsets = integer(2^object@n),
                             PACKAGE="kappalab")$subsets
              
              mu <- .C("cardinal2setfunction", 
                       as.integer(object@n),
                       as.double(object@data),
                       mu = double(2^object@n), 
                       PACKAGE="kappalab")$mu 
              
              mu <- .C("binary2natural", 
                       as.integer(object@n),
                       as.double(mu),
                       as.integer(subsets),
                       mu = double(2^object@n), 
                       PACKAGE="kappalab")$mu

              new("Mobius.set.func", data = mu, subsets = subsets, 
                  n = object@n, k = object@n)
          }
          )

##############################################################################

## Tests the monotonicity of a set function represented by its Mobius transform
setMethod("is.monotone", signature(object = "Mobius.set.func"),
          function(object, verbose = FALSE, epsilon = 1e-9, ...) {

              if (!is.logical(verbose))
                  stop("wrong arguments")
              
              .C("is_monotone_Mobius", 
                 as.integer(object@n),
                 as.integer(object@k), 
                 as.double(object@data),
                 as.integer(object@subsets),
                 as.integer(verbose),
                 as.double(epsilon), 
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0
          }
          )

## Tests whether the set function represented by its Mobius transform is cardinal
setMethod("is.cardinal", signature(object = "Mobius.set.func"),
          function(object, ...) {
              
              .C("is_kcardinal", 
                 as.integer(object@n),
                 as.integer(object@k), 
                 as.double(object@data), 
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0
          }
          )

## Tests whether the set function represented by
## its Mobius transform is k-additive
setMethod("is.kadditive", signature(object = "Mobius.set.func", k = "numeric"),
          function(object, k, epsilon = 1e-9, ...) {
              if (!(k %in% 1:object@n))
                  stop("wrong arguments")
              
              .C("is_kadditive_Mobius", 
                 as.integer(object@n),
                 as.integer(object@k),
                 as.integer(k),
                 as.double(object@data),
                 as.double(epsilon),
                 result = integer(1),
                 PACKAGE="kappalab")$result[1] == 0	
          }
          )

##############################################################################

## Show method for object Mobius.set.func
setMethod("show", signature(object = "Mobius.set.func"),	
          function(object) {

              cat(paste("\t\t",is(object)[1],"\n",sep=""))
              .C("Rprint_setfunction", 
                 as.integer(object@n),
                 as.integer(object@k),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 as.integer(1),
                 PACKAGE="kappalab")
              cat("")
          }
          )

## Displays the Mobius transform
setMethod("to.data.frame", signature(object = "Mobius.set.func"),
          function(object, ...) {

              if (object@n < max.n.display) {
                  
                  subsets <- .C("k_power_set_char", 
                                as.integer(object@n),
                                as.integer(object@k),
                                as.integer(object@subsets),  
                                subsets = character(binom.sum(object@n,object@k)),
                                PACKAGE="kappalab")$subsets
                  
                  d <- data.frame(c(object@n,object@k,object@data),
                                  row.names = c("  n", "  k", subsets))
              }
              else
                  d <- data.frame(c(object@n,object@k,object@data),
                                  row.names = c("n", "k", object@subsets))
              
              names(d)[1] <- is(object)[1]
              d
          }
          )

##############################################################################

## Computes the Shapley value of a set function represented
## by its Mobius transform
setMethod("Shapley.value", signature(object = "Mobius.set.func"),
          function(object, ...) {
              
              result <- .C("Shapley_value_Mobius", 
                           as.integer(object@n),	
                           as.integer(object@k),  
                           as.double(object@data),
                           as.integer(object@subsets), 
                           phi = double(object@n),
                           PACKAGE="kappalab")$phi
              
              names(result) <- 1:object@n
              result	
          }
          )

## Computes the Shapley interaction indices of a set function 
## represented by its Mobius transform
setMethod("interaction.indices", signature(object = "Mobius.set.func"),
          function(object, ...) {
              
              phi <- .C("interaction_indices_Mobius", 
                        as.integer(object@n),
                        as.integer(object@k), 
                        as.double(object@data),
                        as.integer(object@subsets), 
                        result = double(object@n*object@n),
                        PACKAGE="kappalab")$result
              
              result <- matrix(phi,object@n,object@n)
              diag(result) <- rep(NA,object@n)
              dimnames(result) <- list(1:object@n,1:object@n)
              result
          }
          )

##############################################################################
