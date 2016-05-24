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
###############################################################################

## Class card.set.func extends superclass.set.func

##############################################################################


## Constructor from numeric
card.set.func <- function(object) {

    if (!is.numeric(object))
        stop("wrong argument")

    new("card.set.func", data = object, n = length(object) - 1)
}

## Constructor from set.func
setMethod("as.card.set.func", signature(object = "set.func"),
          function(object, ...) {

              csf <- as.card.set.func.internal(object)
              new("card.set.func", data = csf$data, n = csf$n)
          }
          )

## Constructor from set.func
setMethod("as.card.set.func", signature(object = "Mobius.card.set.func"),
          function(object, ...) {

              new("card.set.func", data = object@data, n = object@n)
          }
          )

## Constructor from Mobius.set.func
setMethod("as.card.set.func", signature(object = "Mobius.set.func"),
          function(object, ...) {

              if (!is.cardinal(object))
                  stop("wrong argument")

              data <- .C("setfunction2cardinal",
                         as.integer(object@n),
                         as.integer(object@k),
                         as.double(object@data),
                         data = double(object@n+1),
                         PACKAGE="kappalab")$data

              new("card.set.func", data = data, n = object@n)
          }
          )


## Constructor from numeric
setMethod("zeta", signature(object = "Mobius.card.set.func"),
          function(object, ...) {

              mu <- numeric(object@n+1)
              for (i in 1:(object@n+1))
                  for (j in i:1)
                      mu[i] <- mu[i] + choose(i-1,j-1) * object@data[j]
              new("card.set.func", data = mu, n = object@n)
          }
          )

##############################################################################

## Tests the monotonicity of a set function
setMethod("is.monotone", signature(object = "card.set.func"),
          function(object, verbose = FALSE, epsilon = 1e-9, ...) {

              if (!is.logical(verbose))
                  stop("wrong arguments")

              monotone <- TRUE

              for (i in 1:object@n)
                  if (object@data[i] - object@data[i+1] > epsilon) {

                      if (verbose)
                          cat(paste("Violation of monotonicity constraints between subsets of cardinal",i-1,"and subsets of cardinal",i,".\n"))
                      monotone <- FALSE
                  }
              monotone
          }


          )

## Tests whether the set function is cardinal
setMethod("is.cardinal", signature(object = "card.set.func"),
          function(object,...)  {

              TRUE
          }
          )

## Tests whether the set function is k-additive
setMethod("is.kadditive", signature(object = "card.set.func",k = "numeric"),
          function(object, k, epsilon = 1e-9, ...) {

              if (!(k %in% 1:object@n))
                  stop("wrong arguments")

              a <- Mobius(object)

              if (abs(a@data[k+1]) < epsilon)
                  return(FALSE)

              for (i in (k+2):(object@n+1))
                  if (abs(a@data[i]) > epsilon)
                      return(FALSE)
              return(TRUE)
          }
          )

##############################################################################

## Show method for object card.set.func
#setMethod("show", signature(object = "card.set.func"),
#          function(object) {
#
#              show(to.data.frame(object))
#          }
#          )

## Displays the set function
setMethod("to.data.frame", signature(object = "card.set.func"),
          function(object, ...) {

              d <- as.data.frame(object@data,0:object@n)
              names(d)[1] <- is(object)[1]
              d
          }
          )

##############################################################################

## Computes the Shapley value of a set function
setMethod("Shapley.value", signature(object = "card.set.func"),
          function(object,...) {

              result <- mean(object@data[2:(object@n+1)]
                             - object@data[1:object@n])
              names(result) <- "All"
              result
          }
          )

## Computes the Shapley interaction indices of a set function
setMethod("interaction.indices", signature(object = "card.set.func"),
          function(object,...) {

              result <- matrix(mean(object@data[3:(object@n+1)]
                                    - 2 * object@data[2:object@n]
                                    + object@data[1:(object@n-1)]),1,1)
              dimnames(result) <- list("All","All")
              result
          }
          )

##############################################################################
