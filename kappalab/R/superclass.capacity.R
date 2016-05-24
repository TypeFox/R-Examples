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
###################################################################################

## Virtual class superclass.capacity, superclass of all *.capacity classes

##############################################################################

## Summary method for object set.func
setMethod("summary", signature(object = "superclass.capacity"),
          function(object, ...) {
              
              new("summary.superclass.capacity", 
                  Shapley.value = Shapley.value(object),
                  interaction.indices = interaction.indices(object),
                  orness = orness(object), 
                  veto = veto(object), 
                  favor = favor(object), 
                  variance = variance(object), 
                  entropy = entropy(object)
                  )
          }
          )

## Show method for object summary.superclass.set.func
setMethod("show", signature(object = "summary.superclass.capacity"),
          function(object) {
              
              cat("Shapley value :\n")
              print(object@Shapley.value)
              cat("\nShapley interaction indices :\n")
              print(object@interaction.indices)
              cat("\nOrness :\n")
              cat(object@orness)
              cat("\n\nVeto indices :\n")
              print(object@veto)
              cat("\nFavor indices :\n")
              print(object@favor)
              cat("\nNormalized variance :\n")
              cat(object@variance)
              cat("\n\nNormalized entropy :\n")
              cat(object@entropy)
              cat("\n")
          }
          )

##############################################################################
