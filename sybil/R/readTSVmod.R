#  readTSVmod.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: readTSVmod
#
#
#

readTSVmod <- function(prefix, suffix,
                       reactList, metList = NA, modDesc = NA,
                       fielddelim = "\t", entrydelim = ", ", extMetFlag = "b",
                       excludeComments = TRUE,
                       oneSubSystem = TRUE,
                       mergeMet = TRUE,
                       balanceReact = TRUE,
                       remUnusedMetReact = TRUE,
                       singletonMet = FALSE,
                       deadEndMet = FALSE,
                       remMet = FALSE,
                       constrMet = FALSE,
                       tol = SYBIL_SETTINGS("TOLERANCE"),
                       fpath = SYBIL_SETTINGS("PATH_TO_MODEL"),
                       def_bnd = SYBIL_SETTINGS("MAXIMUM"),
                       arrowlength = NULL,
                       quoteChar = "",
                       commentChar = "",
                       ...) {


    #--------------------------------------------------------------------------#

    if (missing(suffix)) {
        fnEXT <- switch(fielddelim,
                        "\t" = { "tsv" },
                        ";"  = { "csv" },
                        ","  = { "csv" },
                        "|"  = { "dsv" },
                               { "dsv" }
        )
    }
    else {
        fnEXT <- suffix
    }

    # if argument prefix is empty, at least reactList must be not empty.
    if (missing(prefix)) {
        if (missing(reactList)) {
            stop("missing argument 'reactList' is required")
        }
        else {
            fnRL <- reactList
            fnML <- metList
            fnMD <- modDesc
        }
    }
    else {
        # filename
        fnRL <- paste(paste(prefix, "react", sep = "_"), fnEXT, sep = ".")
        fnML <- paste(paste(prefix, "met",   sep = "_"), fnEXT, sep = ".")
        fnMD <- paste(paste(prefix, "desc",  sep = "_"), fnEXT, sep = ".")
    }

    # file path
    fpRL <- file.path(fpath, fnRL)
    fpML <- file.path(fpath, fnML)
    fpMD <- file.path(fpath, fnMD)


    if (!isTRUE(file.exists(fpRL))) {
        stop("failed to open reaction list ", sQuote(fpRL))
    }

    if (!isTRUE(file.exists(fpML))) {
        fpML <- NULL
    }

    if (!isTRUE(file.exists(fpMD))) {
        fpMD <- NULL
    }


    #--------------------------------------------------------------------------#
    # some regular expressions
    #--------------------------------------------------------------------------#

    # delimiter for the compartment flag (round or square bracket)
    #compartDelimL <- "\\["
    #compartDelimR <- "\\]"
    compartDelimL  <- "(\\[|\\()"
    compartDelimR  <- "(\\]|\\))"
    compartCharSet <- "(\\w+)"

    # regular expression for the compartment flag
    compartRegEx    <- paste(compartDelimL, compartCharSet, compartDelimR, sep = "")
    compartRegExEnd <- paste(compartRegEx, "$", sep = "")

    # regular expression to identify external metabolites
    extMetRegEx <- paste(compartDelimL, extMetFlag, compartDelimR, sep = "")
    #extMetRegEx <- paste("\\[", extMetFlag, "\\]$", sep = "")

    
    # regular expression to identify the reaction arrow
    if (is.null(arrowlength)) {
        arrowregex <- "<?[-=]+>"
    }
    else {
        stopifnot(length(arrowlength) == 1)
        if (is.numeric(arrowlength)) {
            arrowregex <- paste("<?[-=]{", arrowlength, "}>", sep = "")
        }
        else if (is.character(arrowlength)) {
            arrowregex <- paste("<?[-=]", arrowlength, ">", sep = "")
        }
    }


    #--------------------------------------------------------------------------#
    # some functions
    #--------------------------------------------------------------------------#

    prepReact <- function(metabolites, educts, rowind) {

        components <- unlist(strsplit(metabolites, "+", fixed = TRUE))
        nc <- length(components)

        # data structures
        nNewMet  <- logical(nc)
        newMetId <- character(nc)
        newMetCo <- character(nc)

        newRowIA <- integer(nc)
        newRowRA <- integer(nc)

        # column indices (reactions, lines in the input file) in S
        newRowJA <- as.integer(rep(rowind, nc))

        # FALSE if a metabolite is not used, e.g. external metabolites,
        # or metabolites which are used more than once as educt or product.
        USE_MET  <- logical(nc)

        CURR_MET <- character(nc) # metabolites in components

        for (j in seq(along = components)) {

            components[j] <- paste(components[j], compfl, sep = "")

            # exclude external metabolites
            if (regexpr(extMetRegEx, components[j],
                        perl = TRUE, ignore.case = TRUE) == -1) {

                USE_MET[j] <- TRUE

                # stoichiometric coefficient must be in round bracket
                stoichpos <- regexpr("^\\([^)]+\\)", components[j], perl = TRUE)
                #stoichpos <- regexpr("^\\(.+\\)", components[j], perl = TRUE)
                if (stoichpos == 1) {
                    st_str <- substr(components[j], 2,
                                         (-1+attr(stoichpos, "match.length")))
                    stoich <- as.numeric(st_str)

                    # the stoichiometric coefficient must be a number
                    if (is.na(stoich)) {
                        msg <- paste("reaction no.", rowind,
                                      dQuote(reactABBR[rowind]),
                                      "stoichiometric coefficient",
                                      dQuote(st_str), "in file", dQuote(fpRL),
                                      "is not numeric, set to '1'.")
                        warning(msg, call. = FALSE)
                        stoich <- 1
                    }

                    currentMet <- substr(components[j],
                                         (1+attr(stoichpos, "match.length")),
                                         nchar(components[j]))
                }
                else {
                    stoich <- 1
                    currentMet <- components[j]
                }

                
                # put the compartment flag into a square bracket
                CURR_MET[j] <- sub(compartRegExEnd, "[\\2]", currentMet, perl = TRUE)
                #CURR_MET[j] <- currentMet

                stoich <- ifelse(isTRUE(educts), (stoich * -1), stoich)


                # cid contains the compartment abbreviation for the current
                # metabolite (excluding '[]')
                cidpos <- regexpr(compartRegExEnd, components[j], perl = TRUE)
                #cidpos <- regexpr("\\[\\w+\\]$", components[j], perl = TRUE)
                if (cidpos < 0) {
                    cid <- "unknown"
                }
                else {
                    cid <- substr(components[j], (cidpos+1),
                           (cidpos+attr(cidpos, "match.length")-2))
                }


                # check if the metabolite is used more than once in
                # the current reaction as educt or product
                if (isTRUE(mergeMet)) {
                    met_indexCURR <- match(CURR_MET[j], CURR_MET[-j])
                }
                else {
                    met_indexCURR <- NA
                }

                if (is.na(met_indexCURR)) {

                    # the non zero element in S for the current metabolite
                    newRowRA[j] <- stoich

                    # check if we have a new metabolite
                    met_index <- match(CURR_MET[j], Rmet)
                    if (is.na(met_index)) {

                        # new metabolite: new row in S
                        nNewMet[j]  <- TRUE
                        newMetId[j] <- CURR_MET[j]
                        newMetCo[j] <- paste("[", cid, "]", sep = "")
                        newRowIA[j] <- as.integer(sum(nNewMet) + NRmet)
                    }
                    else {
                        # existing metabolite
                        nNewMet[j]  <- FALSE
                        newRowIA[j] <- as.integer(met_index)
                    }

                }
                else {
                    # e.g.: if we have a + a, we merge it to (2) a
                    USE_MET[j] <- FALSE
                    newRowRA[met_indexCURR] <- newRowRA[met_indexCURR] + stoich
                    msg <- paste("reaction no.", rowind,
                                 dQuote(reactABBR[rowind]),
                                 "metabolite no.", j, dQuote(CURR_MET[j]),
                                 "was merged")
                    warning(msg, call. = FALSE)
                }


            }
            else {
            
            }
        }

        return(list(nNewMet  = sum(nNewMet),
                    newMetId = newMetId[nNewMet],
                    newMetCo = newMetCo[nNewMet],

                    newRowIA = newRowIA[USE_MET],
                    newRowJA = newRowJA[USE_MET],
                    newRowRA = newRowRA[USE_MET]))

    }


    #--------------------------------------------------------------------------#
    # model description
    #--------------------------------------------------------------------------#

    newModelName <- sub("(_react)?\\.[^.]+$", "", basename(fpRL), perl = TRUE)

    if (is.null(fpMD)) {
        
        modNM     <- newModelName   # name
        modID     <- newModelName   # id
        modDC     <- ""             # description
        modCO     <- NULL           # compartments
        modCA     <- NULL           # compartment abbreviations
        modNmet   <- NULL           # number of metabolites
        modNreact <- NULL           # number of reactions
        modNgenes <- NULL           # number of genes
        modNnnz   <- NULL           # number of non-zero elements in S
    }
    else {
        message("reading model description, ... ", appendLF = FALSE)

        md    <- read.table(fpMD,
                            header = TRUE,
                            sep = fielddelim,
                            quote = quoteChar,
                            comment.char = commentChar, ...)

        # The first line containing data (second line in file) is used.
        # Everything else will be ignored.
        modNM     <- try(as.character(md[ ,"name"]),         silent = TRUE)
        modID     <- try(as.character(md[ ,"id"]),           silent = TRUE)
        modDC     <- try(as.character(md[ ,"description"]),  silent = TRUE)
        modCO     <- try(as.character(md[ ,"compartment"]),  silent = TRUE)
        modCA     <- try(as.character(md[ ,"abbreviation"]), silent = TRUE)
        modNmet   <- try(as.integer(md[ ,"Nmetabolites"]),   silent = TRUE)
        modNreact <- try(as.integer(md[ ,"Nreactions"]),     silent = TRUE)
        modNgenes <- try(as.integer(md[ ,"Ngenes"]),         silent = TRUE)
        modNnnz   <- try(as.integer(md[ ,"Nnnz"]),           silent = TRUE)

        if (is(modNM, "try-error")) {
            warning("field 'name' is obligatory in model description")
            modNM <- newModelName
        }
        else {
            nal <- .checkEmptyField(modNM, "name")
            if (!is.null(nal)) {
                #stop(nal[["msg"]], call. = FALSE)
                modNM <- newModelName
            }
        }

        if (is(modID, "try-error")) {
            warning("field 'id' is obligatory in model description.")
            modID <- newModelName
        }
        else {
            nal <- .checkEmptyField(modID, "id")
            if (!is.null(nal)) {
                #stop(nal[["msg"]], call. = FALSE)
                modID <- newModelName
            }
        }

        if (is(modDC, "try-error")) {
            modDC <- ""
        }
        else {
            nal <- .checkEmptyField(modDC, "description")
            if (!is.null(nal)) {
                modDC <- ""
            }
        }

        if (is(modCO, "try-error")) {
            modCO <- NULL
        }
        else {
            nal <- .checkEmptyField(modCO, "compartment")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modCO <- unlist(strsplit(modCO[1], entrydelim, fixed = TRUE))
        }

        if (is(modCA, "try-error")) {
            modCA <- NULL
        }
        else {
            nal <- .checkEmptyField(modCA, "abbreviation")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modCA <- unlist(strsplit(modCA[1], entrydelim, fixed = TRUE))
        }

        if (is(modNmet, "try-error")) {
            modNmet <- NULL
        }
        else {
            nal <- .checkEmptyField(modNmet, "Nmetabolites")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modNmet <- modNmet[1]
        }

        if (is(modNreact, "try-error")) {
            modNreact <- NULL
        }
        else {
            nal <- .checkEmptyField(modNreact, "Nreactions")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modNreact <- modNreact[1]
        }

        if (is(modNgenes, "try-error")) {
            modNgenes <- NULL
        }
        else {
            nal <- .checkEmptyField(modNgenes, "Ngenes")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modNgenes <- modNgenes[1]
        }

        if (is(modNnnz, "try-error")) {
            modNnnz <- NULL
        }
        else {
            nal <- .checkEmptyField(modNnnz, "Nnnz")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            modNnnz <- modNnnz[1]
        }

        remove(md)

        message("OK")

    }


    #--------------------------------------------------------------------------#
    # new instance of class modelorg
    #--------------------------------------------------------------------------#

    model <- modelorg(modID[1], modNM[1])
    mod_desc(model) <- as.character(modDC[1])


    #--------------------------------------------------------------------------#
    # metabolites list
    #--------------------------------------------------------------------------#

    if (is.null(fpML)) {
        metABBR <- NULL
        metNAME <- NULL
        metCOMP <- NULL
    }
    else {
        message("reading metabolite list ... ", appendLF = FALSE)

        ml <- read.table(fpML,
                         header = TRUE,
                         sep = fielddelim,
                         quote = quoteChar,
                         comment.char = commentChar, ...)

        metABBR <- try(as.character(ml[ ,"abbreviation"]), silent = TRUE)
        metNAME <- try(as.character(ml[ ,"name"]),         silent = TRUE)
        metCOMP <- try(as.character(ml[ ,"compartment"]),  silent = TRUE)

        if (is(metABBR, "try-error")) {
            metABBR <- NULL
        }
        else {
            nal <- .checkEmptyField(metABBR, "abbreviation")
            if (!is.null(nal)) {
                stop(nal[["msg"]], call. = FALSE)
            }
            
            #  [...] containing leading whitespace(s), are treated as comments
            if (isTRUE(excludeComments)) {
                metABBR <- sub("\\s+\\[.+\\]$", "", metABBR, perl = TRUE)
            }
            else {
                metABBR <- gsub("\\s+", "", metABBR, perl = TRUE)
            }
            #metABBR <- gsub(" ", "", metABBR, fixed = TRUE)
        }

        if (is(metNAME, "try-error")) {
            warning("field 'name' does not exist in metabolite list")
            metNAME <- NULL
        }
        else {
            nal <- .checkEmptyField(metNAME, "name")
            if (!is.null(nal)) {
                metNAME[nal[["nalines"]]] <- ""
            }
        }

        if (is(metCOMP, "try-error")) {
            metCOMP <- NULL
        }
        else {
            nal <- .checkEmptyField(metCOMP, "compartment")
            if (!is.null(nal)) {
                metCOMP[nal[["nalines"]]] <- ""
            }
        }

        remove(ml)

        message("OK")
    }


    #--------------------------------------------------------------------------#
    # reactions list
    #--------------------------------------------------------------------------#

    message("parsing reaction list ... ", appendLF = FALSE)

    rl <- read.table(fpRL,
                     header = TRUE,
                     sep = fielddelim,
                     quote = quoteChar,
                     comment.char = commentChar, ...)

    # obligatory
    reactEQU  <- try(as.character(rl[ ,"equation"]), silent = TRUE)

    if (is(reactEQU, "try-error")) {
        stop("field 'equation' does not exist in the reactions list")
    }
    else {
        nal <- .checkEmptyField(reactEQU, "equation")
        if (!is.null(nal)) {
            stop(nal[["msg"]], call. = FALSE)
        }

        nreact <- length(reactEQU)          # number of reactions n

        if (!is.null(modNreact)) {
            if (nreact != modNreact) {
                msg <- paste(nreact," reactions detected,", modNreact,
                             "reactions expected, according to model",
                             "description file", dQuote(fpMD))
                warning(msg, call. = FALSE)
            }
        }
    }

    # optionally
    reactABBR <- try(as.character(rl[ ,"abbreviation"]), silent = TRUE)
    reactNAME <- try(as.character(rl[ ,"name"]),         silent = TRUE)
    reactREV  <- try(as.character(rl[ ,"reversible"]),   silent = TRUE)
    reactCOMP <- try(as.character(rl[ ,"compartment"]),  silent = TRUE)
    reactLOW  <- try(as.numeric(rl[ ,"lowbnd"]),         silent = TRUE)
    reactUPP  <- try(as.numeric(rl[ ,"uppbnd"]),         silent = TRUE)
    reactQBJ  <- try(as.numeric(rl[ ,"obj_coef"]),       silent = TRUE)
    reactRULE <- try(as.character(rl[ ,"rule"]),         silent = TRUE)
    reactSUBS <- try(as.character(rl[ ,"subsystem"]),    silent = TRUE)


    if (is(reactABBR, "try-error")) {
        msg <- paste("field 'abbreviation' does not",
                     "exist in reaction list")
        warning(msg, call. = FALSE)
        # If there is no field 'abbreviation', we use v1, v2, ..., vn as id's.
        reactABBR <- paste(rep("v", nreact), 1:nreact, sep = "")
    }
    else {
        nal <- .checkEmptyField(reactABBR, "abbreviation")
        if (!is.null(nal)) {
            stop(nal[["msg"]], call. = FALSE)
        }

        #  [...] containing leading whitespace(s), are treated as comments
        if (isTRUE(excludeComments)) {
            reactABBR <- sub("\\s+\\[.+\\]$", "", reactABBR, perl = TRUE)
        }
    }

    if (is(reactNAME, "try-error")) {
        reactNAME <- NULL
    }
    else {
        nal <- .checkEmptyField(reactNAME, "name")
        if (!is.null(nal)) {
            reactNAME[nal[["nalines"]]] <- ""
        }
    }

    if (is(reactREV, "try-error")) {      # if empty, use arrow symbol
        reactREV <- NULL                  # in reactEQU
    }
    else {
        nal <- .checkEmptyField(reactREV, "reversible")
        if (!is.null(nal)) {
            stop(nal[["msg"]], call. = FALSE)
        }
    }

    if (is(reactCOMP, "try-error")) {
        reactCOMP <- NULL
    }
    else {
        nal <- .checkEmptyField(reactCOMP, "compartment")
        if (!is.null(nal)) {
            #stop(nal[["msg"]], call. = FALSE)
            reactCOMP[nal[["nalines"]]] <- ""
        }
    }

    if (is(reactLOW, "try-error")) {      # If lowbnd and|or uppbnd is empty,
        reactLOW <- NULL                  # we use for irreversible reactions:
    }                                     # [0, SYBIL_SETTINGS("MAXIMUM")]
    else {
        nal <- .checkEmptyField(reactLOW, "lowbnd")
        if (!is.null(nal)) {
            #stop(nal[["msg"]], call. = FALSE)
            reactLOW[nal[["nalines"]]] <- -1 * def_bnd
        }
    }

    if (is(reactUPP, "try-error")) {      # and for reversible reactions
        reactUPP <- NULL                  # [-SYBIL_SETTINGS("MAXIMUM"),
    }                                     #   SYBIL_SETTINGS("MAXIMUM")]
    else {
        nal <- .checkEmptyField(reactUPP, "uppbnd")
        if (!is.null(nal)) {
            #stop(nal[["msg"]], call. = FALSE)
            reactUPP[nal[["nalines"]]] <- def_bnd
        }
    }

    if (is(reactQBJ, "try-error")) {
        reactQBJ <- rep(0, nreact)
    }
    else {
        nal <- .checkEmptyField(reactQBJ, "obj_coef")
        if (!is.null(nal)) {
            #stop(nal[["msg"]], call. = FALSE)
            reactQBJ[nal[["nalines"]]] <- 0
        }
    }

    if (is(reactRULE, "try-error")) {
        reactRULE <- rep("", nreact)
    }
    else {
        nal <- .checkEmptyField(reactRULE, "rule")
        if (!is.null(nal)) {
            reactRULE[nal[["nalines"]]] <- ""
        }
    }

    if (is(reactSUBS, "try-error")) {
        reactSUBS <- NULL
    }
    else {
        nal <- .checkEmptyField(reactSUBS, "subsystem")
        if (!is.null(nal)) {
            reactSUBS[nal[["nalines"]]] <- "none"
        }
    }

    remove(rl)


    #--------------------------------------------------------------------------#
    # parse reactions
    #--------------------------------------------------------------------------#

    # data structures
    Rrev   <- logical(nreact)
    #Rcomp  <- character(nreact)
    Rlow   <- numeric(nreact)
    Rupp   <- numeric(nreact)
    Rgenes <- vector(mode = "list", length = nreact)
    Rrules <- character(nreact)
    RruleL <- logical(nreact)

    allGenes <- character(0)

    NRmet <- 0                      # number of detected metabolites
    NRnnz <- 0                      # number of current non-zero elements

    if (is.null(modNmet)) {
        Rmet  <- character(0)       # metabolite id's
        RmetC <- integer(0)         # metabolite compartments
    }
    else {
        Rmet  <- character(modNmet)
        RmetC <- integer(modNmet)
    }

    if (is.null(modNreact)) {
        RmatIA <- integer(0)        # row indices
        RmatJA <- integer(0)        # column indices
        RmatRA <- numeric(0)        # non zero elements
    }
    else {
        RmatIA <- integer(modNnnz)
        RmatJA <- integer(modNnnz)
        RmatRA <- numeric(modNnnz)
    }

    SKIP_REACTION <- rep(TRUE, nreact)


    # exclude comments from metabolite id's in reaction strings
    if (isTRUE(excludeComments)) {
        requatC <- gsub("(\\w+)\\s+\\[[^]]+\\]", "\\1", reactEQU, perl = TRUE)
    }
    else {
        requatC <- reactEQU
    }

    # remove whitespaces from reaction string
    requatW <- gsub("\\s+", "", requatC, perl = TRUE)

    # check for possible transport reaction
    transpR <- regexpr(paste("^", compartRegEx, ":", sep = ""), requatW, perl = TRUE)
    #transpR <- regexpr("^\\[\\w+\\]:", requatW, perl = TRUE)


    # indices for the next metabolite, next nnz element
    nextMet <- 1
    nextNnz <- 1

    for (i in seq(along = reactEQU)) {

        # get the compartment flag
        if (transpR[i] > 0) {
            #equat   <- substr(requatW[i], 5, nchar(requatW[i]))
            #compfl  <- substr(requatW[i], 1, 3)

            equat  <- substr(requatW[i],
                             (1+attr(transpR, "match.length")[i]),
                             nchar(requatW[i]))
            compfl <- substr(requatW[i],
                             1, (attr(transpR, "match.length")[i]-1))
        }
        else {
            equat  <- requatW[i]
            compfl <- ""

        }

        # The reaction arrow splits the equation string into educts (on the
        # left) and products (on the right). Each reaction string must contain
        # exactly one reaction arrow. The last sign of the arrow must be a '>',
        # meaning the reaction is written in direction to products.
        # Reversible reactions can begin with a '<' sign.

        # get the position of the arrow symbol
        arrowpos  <- gregexpr(arrowregex, equat, perl = TRUE)[[1]]

        if (length(arrowpos) > 1) {
            msg <- paste("more than one arrow symbols found, skipping",
                         "reaction no.", i, dQuote(reactABBR[i]))
            warning(msg, call. = FALSE)
            SKIP_REACTION[i] <- FALSE
            next
        }

        if (arrowpos < 0) {
            msg <- paste("no reaction arrow found, skipping reaction no.",
                         i, dQuote(reactABBR[i]))
            warning(msg, call. = FALSE)
            SKIP_REACTION[i] <- FALSE
            next
        }

        arrowsymb <- substr(equat, arrowpos,
                            (arrowpos+attr(arrowpos, "match.length"))-1)


        # If reactREV is empty, the reaction arrow is used to determine
        # reversibility.

        if (is.null(reactREV)) {
            if (substr(arrowsymb, 1, 1) == "<") {
                Rrev[i] <- TRUE
            }
        }
        else {
            if ( (isTRUE(reactREV[i])) || (reactREV[i] == "Reversible") ) {
                Rrev[i] <- TRUE
            }
        }


        # The current reaction is treated as 'reversible', if the arrowsymbol
        # starts with a '<' (assuming the arrow is like <==>). If not and
        # 'reactREV' is not empty, the entry on reactREV[i] will be checked.

#        if (substr(arrowsymb, 1, 1) == "<") {
#            Rrev[i] <- TRUE
#        }
#        else {
#            if ( (!is.null(reactREV)) && (
#                 (isTRUE(reactREV[i])) || (reactREV[i] == "Reversible") ) ) {
#                Rrev[i] <- TRUE
#            }
#        }


        # Split reaction string in educts and products
        edprod <- unlist(strsplit(equat, arrowsymb, fixed = TRUE))

        if ( (length(edprod) > 2) || (length(edprod) < 1) ) {
            msg <- paste("something went wrong here, skipping reaction no.",
                         i, dQuote(reactABBR[i]))
            warning(msg)
            SKIP_REACTION[i] <- FALSE
            next
        }

        # Educts
        if (nchar(edprod[1]) > 0) {
            newR  <- prepReact(edprod[1], educts = TRUE, rowind = i)
            if (length(newR$newRowRA) > 0) {
                NRmet <- NRmet + newR$nNewMet
                NRnnz <- NRnnz + length(newR$newRowRA)
    
                if (is.null(modNmet)) {
                    Rmet  <- c(Rmet,  newR$newMetId)
                    RmetC <- c(RmetC, newR$newMetCo)
                }
                else {
                    if (newR$nNewMet > 0) {
                        Rmet[nextMet:NRmet]  <- newR$newMetId
                        RmetC[nextMet:NRmet] <- newR$newMetCo
    
                        nextMet <- NRmet + 1
                    }
                }
    
                if (is.null(modNnnz)) {
                    RmatIA <- c(RmatIA, newR$newRowIA)
                    RmatJA <- c(RmatJA, newR$newRowJA)
                    RmatRA <- c(RmatRA, newR$newRowRA)
                }
                else {
                    RmatIA[nextNnz:NRnnz] <- newR$newRowIA
                    RmatJA[nextNnz:NRnnz] <- newR$newRowJA
                    RmatRA[nextNnz:NRnnz] <- newR$newRowRA
    
                    nextNnz <- NRnnz + 1
                }
            }
        }

        # Products
        if ( (length(edprod) == 2) && (nchar(edprod[2]) > 0) ) {
            newR  <- prepReact(edprod[2], educts = FALSE, rowind = i)
            if (length(newR$newRowRA) > 0) {
                NRmet <- NRmet + newR$nNewMet
                NRnnz <- NRnnz + length(newR$newRowRA)
    
                if (is.null(modNmet)) {
                    Rmet  <- c(Rmet,  newR$newMetId)
                    RmetC <- c(RmetC, newR$newMetCo)
                }
                else {
                    if (newR$nNewMet > 0) {
                        Rmet[nextMet:NRmet]  <- newR$newMetId
                        RmetC[nextMet:NRmet] <- newR$newMetCo
    
                        nextMet <- NRmet + 1
                    }
                }
    
                if (is.null(modNnnz)) {
                    RmatIA <- c(RmatIA, newR$newRowIA)
                    RmatJA <- c(RmatJA, newR$newRowJA)
                    RmatRA <- c(RmatRA, newR$newRowRA)
                }
                else {
                    RmatIA[nextNnz:NRnnz] <- newR$newRowIA
                    RmatJA[nextNnz:NRnnz] <- newR$newRowJA
                    RmatRA[nextNnz:NRnnz] <- newR$newRowRA
    
                    nextNnz <- NRnnz + 1
                }
            }
        }


        # upper and lower bounds
        if (is.null(reactLOW)) {
            Rlow[i] <- ifelse(isTRUE(Rrev[i]), (def_bnd * -1), 0)
        }
        else {
            bound <- as.numeric(reactLOW[i])
            if (abs(bound) > def_bnd) {
                Rlow[i] <- ifelse(bound < 0, def_bnd * -1, def_bnd)
            }
            else {
                Rlow[i] <- bound
            }
            # Rlow[i] <- as.numeric(reactLOW[i])
        }

        if (is.null(reactUPP)) {
            Rupp[i] <- def_bnd
        }
        else {
            bound <- as.numeric(reactUPP[i])
            if (abs(bound) > def_bnd) {
                Rupp[i] <- ifelse(bound < 0, def_bnd * -1, def_bnd)
            }
            else {
                Rupp[i] <- bound
            }
            # Rupp[i] <- as.numeric(reactUPP[i])
        }


        # gpr association
        gene_rule <- .parseBoolean(reactRULE[i])
        Rgenes[[i]] <- gene_rule$gene                # list of involved genes
        Rrules[i]   <- gene_rule$rule                # the rule string
        if (gene_rule$rule != "") {
            allGenes  <- c(allGenes, gene_rule$gene)
            RruleL[i] <- TRUE
        }

    }

    message("OK")


    #--------------------------------------------------------------------------#
    # gene to reaction mapping
    #--------------------------------------------------------------------------#

    message("GPR mapping ... ", appendLF = FALSE)

    allGenes <- unique(allGenes)

    if (!is.null(modNgenes)) {
        if (modNgenes != length(allGenes)) {
            msg <- paste(length(allGenes), "genes detected,", modNgenes,
                         "genes expected, according to model description file",
                         dQuote(fpMD))
            warning(msg, call. = FALSE)
        }
    }

    rxnGeneMat <- Matrix::Matrix(FALSE,
                                 nrow = nreact,
                                 ncol = length(allGenes),
                                 sparse = TRUE)

    for (i in 1 : nreact) {

        if (isTRUE(RruleL[i])) {
            geneInd <- match(Rgenes[[i]], allGenes)
            rxnGeneMat[i, geneInd] <- TRUE

            for (j in 1 : length(geneInd)) {
                pat  <- paste("x(", j, ")", sep = "")
                repl <- paste("x[", geneInd[j], "]", sep = "")

                Rrules[i] <- gsub(pat, repl, Rrules[i], fixed = TRUE)
            }
        }

    }

    message("OK")


    #--------------------------------------------------------------------------#
    # subsystems
    #--------------------------------------------------------------------------#

    message("sub systems ... ", appendLF = FALSE)

    subSysdelim <- ifelse(isTRUE(oneSubSystem), NA, entrydelim)
    ssys <- .prepareSubSysMatrix(reactSUBS,
                                         nreact,
                                         entrydelim = subSysdelim)

    message("OK")


    #--------------------------------------------------------------------------#
    # prepare modelorg
    #--------------------------------------------------------------------------#

    message("prepare modelorg object ... ", appendLF = FALSE)

    if (!is.null(modNmet)) {
        if (modNmet != NRmet) {
            msg <- paste(NRmet, "metabolites detected,", modNmet,
                         "metabolites expected, according to model",
                         "description file", dQuote(fpMD))
            warning(msg, call. = FALSE)
        }
    }


    # Set compartment slot. If modCO is empty (NULL), use the compartment flags
    # from the reaction equation (the abbreviations).
    if (is.null(modCO)) {
        mod_compart(model) <- unique(RmetC)
    }
    else {
        mod_compart(model) <- as.character(modCO)
    }


    #--------------------------------------------------------------------------#
    # stoichiometric matrix

    ## this takes long, try to improve that!! ##
    RmatM <- Matrix::Matrix(0, nrow = NRmet, ncol = nreact, sparse = TRUE)

    for (k in seq(along = RmatRA)) {
        if (isTRUE(SKIP_REACTION[RmatJA[k]])) {

            if (isTRUE(balanceReact)) {
                if (RmatM[RmatIA[k], RmatJA[k]] != 0) {
                    msg <- paste("reaction no.", RmatJA[k],
                                 dQuote(reactABBR[RmatJA[k]]), "metabolite no.",
                                 RmatIA[k], dQuote(Rmet[RmatIA[k]]),
                                 "was balanced")
                    warning(msg, call. = FALSE)
                }

                # add up stoichiometric coefficients --> balancing
                RmatM[RmatIA[k], RmatJA[k]] <- RmatM[RmatIA[k], RmatJA[k]] + RmatRA[k]
            }
            else {
                RmatM[RmatIA[k], RmatJA[k]] <- RmatRA[k]
            }

        }
    }


    #--------------------------------------------------------------------------#
    # search for unused metabolites and unused reactions

    # binary matrix
    #RmatMb <- RmatM != 0
    RmatMb <- abs(RmatM) > tol

    SKIP_METABOLITE   <- Matrix::rowSums(RmatMb) != 0    # TRUE, if a metabolite is used
    UNUSED_REACTION   <- Matrix::colSums(RmatMb) != 0    # TRUE, if a reaction is used

    if (isTRUE(remUnusedMetReact)) {
        did <- "and therefore removed from S:"
    }
    else {
        did <- "in S:"
    }


    #--------------------------------------------------------------------------#
    # empty rows

    if (any(SKIP_METABOLITE == FALSE)) {
        met_list  <- paste(dQuote(Rmet[!SKIP_METABOLITE]), collapse = "\n\t")
        nmet_list <- sum(!SKIP_METABOLITE)
        msg_part  <- paste("not used in any reaction", did)
        msg <- sprintf(ngettext(nmet_list,
                                "%d metabolite is %s %s",
                                "%d metabolites are %s\n\t%s"),
                       nmet_list, msg_part, met_list)
        warning(msg, call. = FALSE)
    }


    #--------------------------------------------------------------------------#
    # empty columns

    if (sum(!UNUSED_REACTION) > 0) {
        ur_list  <- paste(dQuote(reactABBR[!UNUSED_REACTION]),
                          collapse = "\n\t")
        nur_list <- sum(!UNUSED_REACTION)
        msg_part <- paste("not used", did)
        msg <- sprintf(ngettext(nur_list,
                                "%d reaction is %s %s",
                                "%d reactions are %s\n\t%s"),
                       nur_list, msg_part, ur_list)
        warning(msg, call. = FALSE)
    }


    #--------------------------------------------------------------------------#
    # correct SKIP...

    if (!isTRUE(remUnusedMetReact)) {
        SKIP_METABOLITE[!SKIP_METABOLITE] <- TRUE
        #UNUSED_REACTION[!UNUSED_REACTION] <- TRUE
        SKIP_REACTION[!UNUSED_REACTION] <- TRUE
    }
    else {
        SKIP_REACTION[!UNUSED_REACTION]   <- FALSE
    }


    #--------------------------------------------------------------------------#
    # single metabolites

    sing_met   <- rep(NA, nrow(RmatM))
    sing_react <- rep(NA, ncol(RmatM))

    if (isTRUE(singletonMet)) {

        message("identifying reactions containing single metabolites ... ",
                appendLF = FALSE)

        singleton <- .singletonMetabolite(mat = RmatMb)

        sing_met[!singleton$smet]     <- FALSE
        sing_react[!singleton$sreact] <- FALSE
        sing_met[singleton$smet]      <- TRUE
        sing_react[singleton$sreact]  <- TRUE

        # singleton metabolites found?
        if (sum(singleton$smet) > 0) {

            if ( xor(isTRUE(constrMet), isTRUE(remMet)) ) {

                if (isTRUE(constrMet)) {
                    # set to zero
                    did_watm <- "identified"
                    did_watr <- "constrained"
                }
                else {
                    # remove
                    SKIP_METABOLITE[singleton$smet] <- FALSE
                    SKIP_REACTION[singleton$sreact] <- FALSE
                    did_watm <- "removed"
                    did_watr <- "removed"
                }

                met_list  <- paste(dQuote(Rmet[singleton$smet]),
                                   collapse = "\n\t")
                nmet_list <- sum(singleton$smet)
                react_list  <- paste(dQuote(reactABBR[singleton$sreact]),
                                   collapse = "\n\t")
                nreact_list <- sum(singleton$sreact)

                msgm <- sprintf(ngettext(nmet_list,
                                        "%s %d singleton metabolite: %s",
                                        "%s %d singleton metabolites:\n\t%s"),
                                did_watm, nmet_list, met_list)

                msgr <- sprintf(ngettext(nreact_list,
                         "%s %d reaction containing singleton metabolites: %s",
                         "%s %d reactions containing singleton metabolites:\n\t%s"),
                                did_watr, nreact_list, react_list)

                #warning(paste(msgm, msgr, sep = "\n\t "), call. = FALSE)
                warning(msgm, call. = FALSE)
                warning(msgr, call. = FALSE)

            }
            else {

                met_list  <- paste(dQuote(Rmet[singleton$smet]),
                                   collapse = "\n\t")
                nmet_list <- sum(singleton$smet)
                msg <- sprintf(ngettext(nmet_list,
                                   "%d metabolite is singleton in S: %s",
                                   "%d metabolites are singletons in S:\n\t%s"),
                               nmet_list, met_list)
                warning(msg, call. = FALSE)
            }
        }
        else {
            message("nothing found ... ", appendLF = FALSE)
            sing_met   <- logical(nrow(RmatM))
            sing_react <- logical(ncol(RmatM))
        }
    }


    #--------------------------------------------------------------------------#
    # dead end metabolites

    de_met   <- rep(NA, nrow(RmatM))
    de_react <- rep(NA, ncol(RmatM))

    if (isTRUE(deadEndMet)) {

        message("identifying reactions containing dead end metabolites ... ",
                appendLF = FALSE)

        demr <- .deadEndMetabolite(mat   = RmatM,
                                           lb    = Rlow,
                                           exclM = sing_met,
                                           exclR = sing_react,
                                           tol   = tol)

        de_met[!demr$dem]   <- FALSE
        de_react[!demr$der] <- FALSE
        de_met[demr$dem]    <- TRUE
        de_react[demr$der]  <- TRUE

        # dead end metabolites found?
        if (sum(demr$dem) > 0) {

            if ( xor(isTRUE(constrMet), isTRUE(remMet)) ) {

                if (isTRUE(constrMet)) {
                    # set to zero
                    did_watm <- "identified"
                    did_watr <- "constrained"
                }
                else {
                    # remove
                    SKIP_METABOLITE[demr$dem] <- FALSE
                    SKIP_REACTION[demr$der]   <- FALSE
                    did_watm <- "removed"
                    did_watr <- "removed"
                }

                met_list  <- paste(dQuote(Rmet[demr$dem]),
                                   collapse = "\n\t")
                nmet_list <- sum(demr$dem)
                react_list  <- paste(dQuote(reactABBR[demr$der]),
                                   collapse = "\n\t")
                nreact_list <- sum(demr$der)

                msgm <- sprintf(ngettext(nmet_list,
                                        "%s %d dead end metabolite: %s",
                                        "%s %d dead end metabolites:\n\t%s"),
                                did_watm, nmet_list, met_list)

                msgr <- sprintf(ngettext(nreact_list,
                                        "%s %d reaction containing dead end metabolites: %s",
                                        "%s %d reactions containing dead end metabolites:\n\t%s"),
                                did_watr, nreact_list, react_list)

                warning(msgm, call. = FALSE)
                warning(msgr, call. = FALSE)

            }
            else {

                met_list  <- paste(dQuote(Rmet[demr$dem]),
                                   collapse = "\n\t")
                nmet_list <- sum(demr$dem)
                msg <- sprintf(ngettext(nmet_list,
                                       "%d dead end metabolite in S: %s",
                                       "%d dead end metabolites in S:\n\t%s"),
                               nmet_list, met_list)
                warning(msg, call. = FALSE)
            }

        }
        else {
            message("nothing found ... ", appendLF = FALSE)
            de_met   <- logical(nrow(RmatM))
            de_react <- logical(ncol(RmatM))
        }
    }


    #--------------------------------------------------------------------------#
    # S

    RmatM <- RmatM[SKIP_METABOLITE, , drop = FALSE]
    RmatM <- RmatM[ , SKIP_REACTION, drop = FALSE]

    S(model) <- RmatM


    #--------------------------------------------------------------------------#
    # metabolites

    NRmet <- NRmet - sum(!SKIP_METABOLITE)
    Rmet  <- Rmet[SKIP_METABOLITE]

    met_id(model)     <- Rmet[1:NRmet]        # use only detected metabolites
    met_num(model)    <- as.integer(NRmet)

    met_single(model) <- sing_met[SKIP_METABOLITE]
    met_de(model)     <- de_met[SKIP_METABOLITE]


    # If there is a column with metabolite abbreviations, the corresponding
    # metabolite name should also exist. If not, the abbreviation will be used
    # as metabolite name.
    if (is.null(metABBR)) {
        met_name(model) <- Rmet[1:NRmet]
    }
    else {
        metAb <- sub(compartRegExEnd, "", Rmet[1:NRmet], perl = TRUE)
        #metAb <- sub("\\[\\w+\\]$", "", Rmet[1:NRmet], perl = TRUE)
        #metAb <- sub("^0+", "", metAb, perl = TRUE) # remove leading zeros
        metId <- match(metAb, metABBR)
        met_name(model) <- metNAME[metId]
        if (any(is.na(metId))) {
            msg <- paste("did not find some metabolite id's in the",
                         "list of metabolite names", dQuote(fpML),
                         "set to", sQuote(NA))
            warning(msg, call. = FALSE)
        }
    }

    # Set the metabolite <-> compartment accociations. If modCA is empty (NULL),
    # use the compartment flags from the reaction equation (the abbreviations).
    # is not empty.
    if (is.null(modCA)) {
        met_comp(model)    <- match(RmetC[1:NRmet], mod_compart(model))
    }
    else {
        met_comp(model)    <- match(RmetC[1:NRmet], modCA)
    }


    #--------------------------------------------------------------------------#
    # reactions

    react_id(model)   <- as.character(reactABBR[SKIP_REACTION])
    react_rev(model)  <- Rrev[SKIP_REACTION]
    react_num(model)  <- as.integer(sum(SKIP_REACTION))
    if (is.null(reactNAME)) {
        react_name(model) <- as.character(reactABBR[SKIP_REACTION])
    }
    else {
        react_name(model) <- as.character(reactNAME[SKIP_REACTION])
    }

    react_single(model) <- sing_react[SKIP_REACTION]
    react_de(model)     <- de_react[SKIP_REACTION]

    if (isTRUE(constrMet)) {
        Rupp[sing_react] <- 0
        Rlow[sing_react] <- 0
        Rupp[de_react]   <- 0
        Rlow[de_react]   <- 0
    }
    else {}

    uppbnd(model)     <- Rupp[SKIP_REACTION]
    lowbnd(model)     <- Rlow[SKIP_REACTION]
    obj_coef(model)   <- as.numeric(reactQBJ[SKIP_REACTION])


    #--------------------------------------------------------------------------#
    # genes

    allGenes(model)   <- allGenes
    genes(model)      <- Rgenes[SKIP_REACTION]
    gprRules(model)   <- Rrules[SKIP_REACTION]
    gpr(model)        <- reactRULE[SKIP_REACTION]
    rxnGeneMat(model) <- rxnGeneMat[SKIP_REACTION, , drop = FALSE]

    #subSys(model)     <- as.character(reactSUBS[SKIP_REACTION])
    subSys(model)     <- ssys[SKIP_REACTION, , drop = FALSE]


    message("OK")


    #--------------------------------------------------------------------------#
    # validate model
    #--------------------------------------------------------------------------#

    message("validating object ... ", appendLF = FALSE)

    check <- validObject(model, test = TRUE)

    if (check != TRUE) {
        msg <- paste("Validity check failed:", check, sep = "\n    ")
        warning(msg)
    }

    message("OK")


    #--------------------------------------------------------------------------#
    # return model
    #--------------------------------------------------------------------------#


    return(model)

}
