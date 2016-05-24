#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  sybilSBML.R
#  Link to libSBML for sybil.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybilSBML.
#
#  SybilSBML is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SybilSBML is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SybilSBML.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#

versionLibSBML <- function() {

    version <- .Call("getLibSBMLversion", PACKAGE = "sybilSBML")
    return(version)

}


#------------------------------------------------------------------------------#

openSBMLfile <- function(fname, ptrtype = "sbml_doc") {

    if ( file.exists(fname) == FALSE ) {
        stop("file not found: ", sQuote(fname))
    }

    sbmlf <- .Call("readSBMLfile", PACKAGE = "sybilSBML",
                   as.character(normalizePath(fname)[1]),
                   as.character(ptrtype)
             )

    sbmlfP <- sbmlDocPointer(sbmlf)

    return(sbmlfP)
}


#------------------------------------------------------------------------------#

closeSBMLfile <- function(sbmlf) {

    invisible(
        .Call("delDocument", PACKAGE = "sybilSBML",
              sbmlPointer(sbmlf)
        )
    )

}


#------------------------------------------------------------------------------#

getSBMLmodel <- function(sbmlf, ptrtype = "sbml_mod") {

    sbmlm <- .Call("getSBMLmodel", PACKAGE = "sybilSBML",
                   sbmlPointer(sbmlf),
                   as.character(ptrtype)
             )

    sbmlmP <- sbmlModPointer(sbmlm, sbmlf)

    if (isTRUE(isNULLpointerSBML(sbmlmP))) {
        sbmlmP <- NULL
    }

    return(sbmlmP)
}


#------------------------------------------------------------------------------#

delSBMLmodel <- function(sbmlm) {

    invisible(
        .Call("delModel", PACKAGE = "sybilSBML",
              sbmlPointer(sbmlm)
        )
    )

}


#------------------------------------------------------------------------------#

getSBMLlevel <- function(sbmlf) {

    level <- .Call("getSBMLlevel", PACKAGE = "sybilSBML",
                   sbmlPointer(sbmlf)
             )

    return(level)
}



#------------------------------------------------------------------------------#

getSBMLversion <- function(sbmlf) {

    version <- .Call("getSBMLversion", PACKAGE = "sybilSBML",
                     sbmlPointer(sbmlf)
               )

    return(version)
}


#------------------------------------------------------------------------------#

validateSBMLdocument <- function(sbmlf) {

    if (is(sbmlf, "character")) {
        sbmlff <- openSBMLfile(fname = sbmlf)
    }
    else {
        sbmlff <- sbmlf
    }
    
    val <- .Call("validateDocument", PACKAGE = "sybilSBML",
                 sbmlPointer(sbmlff)
           )

    if (is(sbmlf, "character")) {
        val <- getSBMLerrors(sbmlff)
        closeSBMLfile(sbmlff)
    }

    return(val)
}


#------------------------------------------------------------------------------#

getSBMLerrors <- function(sbmlf) {

    err <- .Call("getSBMLerrors", PACKAGE = "sybilSBML",
                 sbmlPointer(sbmlf)
           )

    err <- sbmlError(err, sbmlf)

    return(err)
}


#------------------------------------------------------------------------------#

getSBMLmodId <- function(sbmlm) {

    modid <- .Call("getSBMLmodId", PACKAGE = "sybilSBML",
                   sbmlPointer(sbmlm)
             )

    return(modid)
}


#------------------------------------------------------------------------------#

getSBMLmodName <- function(sbmlm) {

    modn <- .Call("getSBMLmodName", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(modn)
}


#------------------------------------------------------------------------------#

getSBMLnumCompart <- function(sbmlm) {

    num <- .Call("getSBMLnumCompart", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(num)
}


#------------------------------------------------------------------------------#

getSBMLnumSpecies <- function(sbmlm) {

    num <- .Call("getSBMLnumSpecies", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(num)
}


#------------------------------------------------------------------------------#

getSBMLnumReactions <- function(sbmlm) {

    num <- .Call("getSBMLnumReactions", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(num)
}


#------------------------------------------------------------------------------#

getSBMLunitDefinitionsList <- function(sbmlm) {

    units <- .Call("getSBMLunitDefinitionsList", PACKAGE = "sybilSBML",
                   sbmlPointer(sbmlm)
             )

    return(units)
}


#------------------------------------------------------------------------------#

getSBMLCompartList <- function(sbmlm) {

    comp <- .Call("getSBMLCompartList", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(comp)
}


#------------------------------------------------------------------------------#

getSBMLSpeciesList <- function(sbmlm) {

    spec <- .Call("getSBMLSpeciesList", PACKAGE = "sybilSBML",
                  sbmlPointer(sbmlm)
            )

    return(spec)
}


#------------------------------------------------------------------------------#

getSBMLReactionsList <- function(sbmlm) {

    react <- .Call("getSBMLReactionsList", PACKAGE = "sybilSBML",
                   sbmlPointer(sbmlm)
             )

    return(react)
}









