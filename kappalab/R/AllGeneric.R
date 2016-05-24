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

## Define those generics that we need, if they don't exist

##############################################################################

if (!isGeneric("as.set.func")) {
    setGeneric("as.set.func", function(object, ...) standardGeneric("as.set.func"))
}

if (!isGeneric("as.game")) {
    setGeneric("as.game", function(object, ...) standardGeneric("as.game"))
}

if (!isGeneric("as.capacity")) {
    setGeneric("as.capacity", function(object, ...) standardGeneric("as.capacity"))
}

##############################################################################

if (!isGeneric("as.Mobius.set.func")) {
    setGeneric("as.Mobius.set.func", function(object, ...)
               standardGeneric("as.Mobius.set.func"))
}

if (!isGeneric("as.Mobius.game")) {
    setGeneric("as.Mobius.game", function(object, ...)
               standardGeneric("as.Mobius.game"))
}

if (!isGeneric("as.Mobius.capacity")) {
    setGeneric("as.Mobius.capacity", function(object, ...)
               standardGeneric("as.Mobius.capacity"))
}

##############################################################################

if (!isGeneric("as.card.set.func")) {
    setGeneric("as.card.set.func", function(object, ...)
               standardGeneric("as.card.set.func"))
}

if (!isGeneric("as.card.game")) {
    setGeneric("as.card.game", function(object, ...)
               standardGeneric("as.card.game"))
}

if (!isGeneric("as.card.capacity")) {
    setGeneric("as.card.capacity", function(object, ...)
               standardGeneric("as.card.capacity"))
}

##############################################################################

if (!isGeneric("as.Mobius.card.set.func")) {
    setGeneric("as.Mobius.card.set.func", function(object, ...)
               standardGeneric("as.Mobius.card.set.func"))
}

##############################################################################

if (!isGeneric("is.monotone")) {
    setGeneric("is.monotone", function(object, ...) standardGeneric("is.monotone"))
}

if (!isGeneric("is.cardinal")) {
    setGeneric("is.cardinal", function(object, ...)
               standardGeneric("is.cardinal"))
}

if (!isGeneric("is.kadditive")) {
    setGeneric("is.kadditive", function(object,k, ...)
               standardGeneric("is.kadditive"))
}

if (!isGeneric("is.normalized")) {
    setGeneric("is.normalized", function(object, ...)
               standardGeneric("is.normalized"))
}

##############################################################################

if (!isGeneric("conjugate")) {
    setGeneric("conjugate", function(object, ...) standardGeneric("conjugate"))
}

if (!isGeneric("Mobius")) {
    setGeneric("Mobius", function(object, ...) standardGeneric("Mobius"))
}

if (!isGeneric("zeta")) {
    setGeneric("zeta", function(object, ...) standardGeneric("zeta"))
}

if (!isGeneric("k.truncate.Mobius")) {
    setGeneric("k.truncate.Mobius", function(object,k, ...)
               standardGeneric("k.truncate.Mobius"))
}

if (!isGeneric("normalize")) {
    setGeneric("normalize", function(object, ...)
               standardGeneric("normalize"))
}

##############################################################################

if (!isGeneric("to.data.frame")) {
    setGeneric("to.data.frame", function(object, ...)
               standardGeneric("to.data.frame"))
}

##############################################################################

if (!isGeneric("rnd")) {
    setGeneric("rnd", function(x, digits = 0)
               standardGeneric("rnd"))
}

if (!isGeneric("summary")) {
    setGeneric("summary", function(object, ...)
               standardGeneric("summary"))
}

##############################################################################

if (!isGeneric("Shapley.value")) {
    setGeneric("Shapley.value", function(object, ...)
               standardGeneric("Shapley.value"))
}

if (!isGeneric("interaction.indices")) {
    setGeneric("interaction.indices", function(object, ...)
               standardGeneric("interaction.indices"))
}

##############################################################################

if (!isGeneric("Choquet.integral")) {
    setGeneric("Choquet.integral", function(object,f, ...)
               standardGeneric("Choquet.integral"))
}

if (!isGeneric("Sugeno.integral")) {
    setGeneric("Sugeno.integral", function(object,f, ...)
               standardGeneric("Sugeno.integral"))
}

if (!isGeneric("Sipos.integral")) {
    setGeneric("Sipos.integral", function(object,f, ...)
               standardGeneric("Sipos.integral"))
}

if (!isGeneric("cdf.Choquet.unif")) {
    setGeneric("cdf.Choquet.unif", function(object,y, ...)
               standardGeneric("cdf.Choquet.unif"))
}

if (!isGeneric("pdf.Choquet.unif")) {
    setGeneric("pdf.Choquet.unif", function(object,y, ...)
               standardGeneric("pdf.Choquet.unif"))
}

if (!isGeneric("pdf.Choquet.exp")) {
    setGeneric("pdf.Choquet.exp", function(object,y, ...)
               standardGeneric("pdf.Choquet.exp"))
}

if (!isGeneric("expect.Choquet.unif")) {
    setGeneric("expect.Choquet.unif", function(object, ...)
               standardGeneric("expect.Choquet.unif"))
}

if (!isGeneric("sd.Choquet.unif")) {
    setGeneric("sd.Choquet.unif", function(object, ...)
               standardGeneric("sd.Choquet.unif"))
}

if (!isGeneric("expect.Choquet.norm")) {
    setGeneric("expect.Choquet.norm", function(object, ...)
               standardGeneric("expect.Choquet.norm"))
}

if (!isGeneric("sd.Choquet.norm")) {
    setGeneric("sd.Choquet.norm", function(object, ...)
               standardGeneric("sd.Choquet.norm"))
}
##############################################################################

if (!isGeneric("veto")) {
    setGeneric("veto", function(object, ...)
               standardGeneric("veto"))
}

if (!isGeneric("favor")) {
    setGeneric("favor", function(object, ...)
               standardGeneric("favor"))
}

if (!isGeneric("orness")) {
    setGeneric("orness", function(object, ...)
               standardGeneric("orness"))
}

if (!isGeneric("variance")) {
    setGeneric("variance", function(object, ...)
               standardGeneric("variance"))
}

if (!isGeneric("entropy")) {
    setGeneric("entropy", function(object, ...)
               standardGeneric("entropy"))
}

##############################################################################
