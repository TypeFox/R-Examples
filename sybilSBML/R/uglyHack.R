#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  uglyHack.R
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


################################################
# Function: .uglyHack
#
#
#

.uglyHack <- function(filename, remapply = FALSE) {

    if (file.exists(filename) == FALSE ) {
        stop( c("cannot open file ", filename) )
    }

    #hackedFileName <- paste(filename, "_temp", sep = "")
    hackedFileName <- tempfile(pattern = basename(filename), fileext = "_temp")

    if (file.exists(hackedFileName) == TRUE ) {
        unlink(hackedFileName)
    }

    lof <- readLines(filename, warn = FALSE)

    #versionLine <- grep("<sbml xmlns=\"", lof, fixed = TRUE)
    #print(versionLine)

    # ------------------------------------------------------------------------ #
    # correct units

    # SBMLError
    # The units of the 'math' formula in a <kineticLaw> definition are expected
    # to be the equivalent of _substance per time_. Expected units are mole
    # (exponent = 1, multiplier = 1, scale = 0),
    # second (exponent = -1, multiplier = 1, scale = 0)
    # but the units returned by the <kineticLaw>'s <math> expression are
    # mole (exponent = 1, multiplier = 1, scale = -3),
    # gram (exponent = -1, multiplier = 1, scale = 0),
    # second (exponent = -1, multiplier = 0.00027777, scale = 0).
    gramLine <- grep("<unit kind=\"gram\"", lof, fixed = TRUE)
    if (length(gramLine) > 0) {
        lof <- lof[-gramLine]
    }
    
    ## SBMLWarning
    # As a principle of best modeling practice, the units of a <parameter>
    # should be declared rather than be left undefined. Doing so improves the
    # ability of software to check the consistency of units and helps make it
    # easier to detect potential errors in models.
#     lof <- sub("(parameter .+OBJECTIVE_COEFFICIENT[^/]+)>/",
#                "\\1 units=\"item\"", lof, fixed = FALSE)
#     lof <- sub("(parameter .+REDUCED_COST[^/]+)",
#                "\\1 units=\"dimensionless\"", lof, fixed = FALSE)

    ## SBMLWarning
    # As a principle of best modeling practice, the size of a <compartment>
    # should be set to a value rather than be left undefined. Doing so improves
    # the portability of models between different simulation and analysis
    # systems, and helps make it easier to detect potential errors in models.
#     lof <- sub("(compartment id=[^/]+)",
#                "\\1 spatialDimensions=\"0\"", lof, fixed = FALSE)

    #lof <- sub("<notes>",
    #           "<notes><body xmlns=\"http://www.w3.org/1999/xhtml\">",
    #           lof, fixed = TRUE)
    #lof <- sub("</notes>", "</body></notes>", lof, fixed = TRUE)


    # ------------------------------------------------------------------------ #
    # The following containers are all optional in a <reaction>,
    # but if any is present, it must not be empty:
    # <listOfReactants>,
    # <listOfProducts>,
    # <listOfModifiers>,
    # <kineticLaw>.

    # this is for Bs_iYO844_flux1.xml
    # this is for iNJ661_biomass9_0.052_middlebrook7H9_plus_glc_glyc_flux.xml
    # this is for Bs_iYO844_flux1.xml


    remLine <- rep(TRUE, length(lof))

    reactLineB  <- grep("<listOfReactants>",  lof, fixed = TRUE)
    reactLineE  <- grep("</listOfReactants>", lof, fixed = TRUE)
    prodLineB   <- grep("<listOfProducts>",   lof, fixed = TRUE)
    prodLineE   <- grep("</listOfProducts>",  lof, fixed = TRUE)
    modLineB    <- grep("<listOfModifiers>",  lof, fixed = TRUE)
    modLineE    <- grep("</listOfModifiers>", lof, fixed = TRUE)
    kinLawLineB <- grep("<kineticLaw>",       lof, fixed = TRUE)
    kinLawLineE <- grep("</kineticLaw>",      lof, fixed = TRUE)

    reactLine  <- reactLineE  - reactLineB
    prodLine   <- prodLineE   - prodLineB
    modLine    <- modLineE    - modLineB
    kinLawLine <- kinLawLineE - kinLawLineB

    noReact  <- which(reactLine  == 1)
    noProd   <- which(prodLine   == 1)
    noMod    <- which(modLine    == 1)
    noKinLaw <- which(kinLawLine == 1)

    remLine[c(reactLineE[noReact],   reactLineB[noReact])]   <-FALSE
    remLine[c(prodLineE[noProd],     prodLineB[noProd])]     <-FALSE
    remLine[c(modLineE[noMod],       modLineB[noMod])]       <-FALSE
    remLine[c(kinLawLineE[noKinLaw], kinLawLineB[noKinLaw])] <-FALSE

    lof <- lof[remLine]


    # ------------------------------------------------------------------------ #
    # The syntax of 'id' attribute values must conform to the syntax of the
    # SBML type 'SId'.
    
    # this is for iNJ661_biomass9_plus_vitamins_0.052_flux.xml (model id contains a .)
    
    lof <- gsub("(id=\"\\w+)\\.(\\w+\")", "\\1_\\2", lof, perl = TRUE)


    # ------------------------------------------------------------------------ #
    # A <reaction> definition must contain at least one <speciesReference>,
    # either in its <listOfReactants> or its <listOfProducts>.
    # A reaction without any reactant or product species is not permitted,
    # regardless of whether the reaction has any modifier species.

    # this is for iNJ661_biomass9_0.052_middlebrook7H9_plus_glc_glyc_flux.xml
    # this is for macModel.xml

    remReact   <- rep(TRUE, length(lof))
    reactBegin <- grep("<reaction.+id=", lof, perl = TRUE)
    reactEnd   <- grep("</reaction>", lof, fixed = TRUE)
    for (i in seq(along = reactBegin)) {
        tmp <- lof[reactBegin[i]:reactEnd[i]]
        noSpRef <- grep("<speciesReference", tmp, fixed = TRUE)
        if (length(noSpRef) == 0) {
            remReact[reactBegin[i]:reactEnd[i]] <- FALSE
        }
    }

    lof <- lof[remReact]


    # ------------------------------------------------------------------------ #
    # Outside of a <functionDefinition>, if a 'ci' element is the first element
    # within a MathML 'apply', then the 'ci''s value can only be chosen from the
    # set of identifiers of <functionDefinition>s defined in the SBML model. 
    # The formula 'LOWER_BOUND(UPPER_BOUND, OBJECTIVE_COEFFICIENT, FLUX_VALUE,
    # REDUCED_COST)' in the math element of the KineticLaw uses 'LOWER_BOUND'
    # which is not a function definition id.

    # this is for Ec_iJR904_GlcMM.xml

    if (isTRUE(remapply)) {
        apb <- grep("<apply>",  lof, fixed = TRUE)
        lof <- lof[-apb]
        ape <- grep("</apply>", lof, fixed = TRUE)
        lof <- lof[-ape]
    }


    # ------------------------------------------------------------------------ #
    # yeast model version < 4.05
    # An SBML XML document must conform to the XML Schema for the corresponding
    # SBML Level, Version and Release. The XML Schema for SBML defines the basic
    # SBML object structure, the data types used by those objects, and the order
    # in which the objects may appear in an SBML document.
    # A KineticLaw object must contain exactly one MathML <math> element.

    # this is for yeast_4.04.xml (and lower versions;
    # yeast_4.05 and higher are ok)

    # translate
    #     <kineticLaw>
    #         ...
    #     </kineticLaw>
    #
    # to
    #     <kineticLaw>
    #         <math xmlns="http://www.w3.org/1998/Math/MathML">
    #             <ci> FLUX_VALUE </ci>
    #         </math>
    #         ...
    #     </kineticLaw>


    # ------------------------------------------------------------------------ #
    # write the new, hopefully error free, file

    hackedFile <- file(hackedFileName, "w")

    #writeLines(lof, con = hackedFile)
    cat(lof, file = hackedFile, sep = "\n")
    close(hackedFile)

    return(hackedFileName)

}
