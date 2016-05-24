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

## virtual class, superclass of all *.set.func classes
setClass("superclass.set.func",
         representation(data = "numeric", n = "numeric"),
	 contains = "VIRTUAL"
         ) 

## virtual class, superclass of all *.capacity classes 
setClass("superclass.capacity", contains = "VIRTUAL")

###############################################################################

## class representing a set function
setClass("set.func",
         representation(subsets = "numeric"),
         contains = "superclass.set.func",
         validity = function(object) {
             
             2^object@n == length(object@data)
         }
         )

## class representing a game
setClass("game",
         contains = "set.func",
         validity = function(object) {
             
             object@data[1] == 0
         }
         )

## class representing a not necessarily normalized capacity
setClass("capacity",
         contains = c("game","superclass.capacity"),
         validity = function(object) {
             
             is.monotone(object)
         }
         )

###############################################################################

## class representing the Mobius transform of a set function
setClass("Mobius.set.func",
	 representation(subsets = "numeric", k = "numeric"),
         contains = "superclass.set.func",
         validity = function(object) {
             binom.sum(object@n, object@k) == length(object@data) &&
             length(object@data) == length(object@subsets)
         }
         )

## class representing the Mobius transform of a game
setClass("Mobius.game",
         contains = "Mobius.set.func",
         validity = function(object) {
             
             object@data[1] == 0
         }
         )

## class representing the Mobius transform of a not necessarily 
## normalized cardinal capacity
setClass("Mobius.capacity",
         contains = c("Mobius.game","superclass.capacity"),
         validity = function(object) {
             
             is.monotone(object)
         }
         )

###############################################################################

## class representing a cardinal set function
setClass("card.set.func",
         contains = "superclass.set.func",
         validity = function(object) {
             
             object@n + 1 == length(object@data)
         }
         )

## class representing a cardinal game
setClass("card.game",
         contains = "card.set.func",
         validity = function(object) {
             
             object@data[1] == 0
         }
         )

## class representing a not necessarily normalized cardinal capacity
setClass("card.capacity",
         contains = c("card.game","superclass.capacity"),
         validity = function(object) {
             
             is.monotone(object)
         }
         )

###############################################################################

## class representing the Mobius transform of a cardinal set function
setClass("Mobius.card.set.func",
         contains = "superclass.set.func",
         validity = function(object) {
             
             object@n + 1 == length(object@data)
         }
         )

###############################################################################

## summary object for superclass.set.func
setClass("summary.superclass.set.func",
         representation(Shapley.value = "numeric", interaction.indices = "matrix")
         )

## summary object for superclass.capacity
setClass("summary.superclass.capacity",
         representation(Shapley.value = "numeric", 
			interaction.indices = "matrix", 
			orness = "numeric", 
			veto = "numeric", 
			favor = "numeric", 
			variance = "numeric", 
			entropy = "numeric")
         )

##############################################################################
