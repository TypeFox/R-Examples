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

## Internal functions

##############################################################################

## Internal util functions

## Max n display
max.n.display <- 7


## Computes the sum of binomial coefficients
binom.sum <- function(n,k) {

    bs <- 1
    for (i in 1:k)
        bs <- bs + choose(n,i)
    bs
}

## Simple test function
is.positive <- function(v) {
    
    is.double(v) && sum(v >= 0) == length(v)
}

##############################################################################

## set.func-class internals

## Creates a list from a numeric for creation of a set.func object
set.func.internal <- function(object) {
    
    if (!is.numeric(object))
        stop("wrong argument")

    n <- log2(length(object))
    if (!(as.integer(n) == n))		
        stop("wrong argument")

    subsets <-  .C("k_power_set", 
                   as.integer(n),
                   as.integer(n),
                   subsets = integer(2^n),
                   PACKAGE="kappalab")$subsets
    
    data <- .C("natural2binary", 
               as.integer(n),
               as.double(object),
               as.integer(subsets),
               data = double(2^n), 
               PACKAGE="kappalab")$data
    
    list(data = data, subsets = subsets, n = n)
}

## Creates a list from a Mobius.set.func for creation of a set.func object
zeta.internal <- function(object) {
    
    data <- .C("Mobius2setfunction", 
               as.integer(object@n),
               as.integer(object@k),
               as.double(object@data), 
               as.integer(object@subsets),
               data = double(2^object@n), 
               PACKAGE="kappalab")$data
    
    if (object@k < object@n)
        subsets <-  .C("k_power_set", 
                       as.integer(object@n),
                       as.integer(object@n),
                       subsets = integer(2^object@n),
                       PACKAGE="kappalab")$subsets
    else
        subsets <- object@subsets

    list(data = data, subsets = subsets, n = object@n)
}

## Creates a list from a card.set.func for creation of a set.func object
as.set.func.internal <- function(object) {
    
    subsets <-  .C("k_power_set", 
                   as.integer(object@n),
                   as.integer(object@n),
                   subsets = integer(2^object@n),
                   PACKAGE="kappalab")$subsets

    data <- .C("cardinal2setfunction", 
               as.integer(object@n),
               as.double(object@data), 
               data = double(2^object@n), 
               PACKAGE="kappalab")$data 

    list(data = data, subsets = subsets, n = object@n)
}

## Creates a list from a set.func for creation of a set.func object
conjugate.internal <- function(object) {
    
    data <- .C("setfunction2conjugate", 
               as.double(object@data),
               as.integer(object@n),
               data = double(2^object@n), 
               PACKAGE="kappalab")$data
    
    list(data = data, subsets = object@subsets, n = object@n)
}

################################################################################

## Mobius.set.func-class internals

## Creates a list from a numeric for creation of a Mobius.set.func object
Mobius.set.func.internal <- function(object, n, k) {
    
    if (!(is.numeric(object) && as.integer(n) == n && k %in% 1:n))
        stop("wrong arguments")
    
    subsets <-  .C("k_power_set", 
                   as.integer(n),
                   as.integer(k),
                   subsets = integer(binom.sum(n, k)),
                   PACKAGE="kappalab")$subsets
    
    list(data = object, n = n, subsets = subsets, k = k)
}

## Creates a list from a set.func for creation of a Mobius.set.func object
k.truncate.Mobius.internal <- function(object, k) {
    
    if (!(k %in% 1:object@n))
        stop("wrong arguments")

    bs <- binom.sum(object@n, k)

    data <- .C("k_truncation", 
               as.integer(object@n),
               as.integer(k),
               as.double(object@data),  								as.integer(object@subsets),
               data = double(bs),
               PACKAGE="kappalab")$data

    list(data = data, n = object@n, 
         subsets = object@subsets[1:bs], k = k)
}

################################################################################

## card.set.func-class internals

## Creates a list from a set.func for creation of a card.set.func object
as.card.set.func.internal <- function(object) {
    
    if (!is.cardinal(object))
        stop("wrong argument")

    mu <- .C("binary2natural", 
             as.integer(object@n),
             as.double(object@data),
             as.integer(object@subsets),
             mu = double(2^object@n), 
             PACKAGE="kappalab")$mu

    data <- .C("setfunction2cardinal", 
               as.integer(object@n),
               as.integer(object@n),
               as.double(mu),
               data = double(object@n+1),
               PACKAGE="kappalab")$data

    list(data = data, n = object@n)
}

##############################################################################
