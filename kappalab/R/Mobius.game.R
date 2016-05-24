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

## Class Mobius.game extends Mobius.set.func

##############################################################################

## Constructor from numeric
Mobius.game <- function(object, n, k)	{
    msf <- Mobius.set.func.internal(object, n, k)
    new("Mobius.game", data = msf$data, subsets = msf$subsets, 
        n = msf$n, k=msf$k)	
}

## Constructor from Mobius.set.func
setMethod("as.Mobius.game", signature(object = "Mobius.set.func"),
          function(object, ...)	{
              new("Mobius.game", data = object@data, subsets = object@subsets,
                  n = object@n, k = object@k)
          }
          )

## Constructor from game
setMethod("Mobius", signature(object = "game"),
          function(object, ...)	{
              msf <- k.truncate.Mobius.internal(object, object@n)
              new("Mobius.game", data = msf$data, subsets = msf$subsets, 
                  n = msf$n, k=msf$k)
          }
          )

##############################################################################

## Computes the Choquet integral of f w.r.t a
setMethod("Choquet.integral",signature(object = "Mobius.game", f = "numeric"),
          function(object,f,...) {
              
              if (!(object@n == length(f)))
                  stop("wrong arguments")
              
              .C("Choquet_integral_Mobius", 
                 as.integer(object@n),
                 as.integer(object@k), 
                 as.double(object@data), 
                 as.integer(object@subsets),
                 as.double(f),  
                 res = double(1),
                 PACKAGE="kappalab")$res[1]
          }
          )

## Computes the Sugeno integral of f w.r.t a (f >= 0)
setMethod("Sugeno.integral",signature(object = "Mobius.game", f = "numeric"),
          function(object,f,...) {
              
              if (!(is.positive(f) && object@n == length(f)))
                  stop("wrong arguments")

              if (is(object)[1] == "Mobius.capacity" && sum(f <=
                    sum(object@data)) != object@n)
                  stop("wrong arguments")
              
              .C("Sugeno_integral_Mobius",
                 as.integer(object@n),
                 as.integer(object@k), 
                 as.double(object@data), 
                 as.integer(object@subsets),
                 as.double(f),     
                 res = double(1),
                 PACKAGE="kappalab")$res[1]
          }
          )

## Computes the Sipos integral of f w.r.t a 
setMethod("Sipos.integral",signature(object = "Mobius.game", f = "numeric"),
          function(object,f,...) {
              if (!(is.double(f) && object@n == length(f)))
                  stop("wrong arguments")
              
              return(Choquet.integral(object,pmax(f,0))
                     - Choquet.integral(object,pmax(-f,0)))
          }
          )

##############################################################################

## Computes the expectation of the Choquet integral in the normal case
setMethod("expect.Choquet.norm",signature(object = "Mobius.game"),
          function(object,...) {
      
              .C("expectation_Choquet_norm_Mobius", 
                 as.integer(object@n), 
                 as.integer(object@k),
                 as.double(object@data),
                 as.integer(object@subsets),  
                 res = double(1),
                 PACKAGE="kappalab")$res
              
          }
          )
