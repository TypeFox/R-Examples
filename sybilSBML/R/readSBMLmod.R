#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  readSBMLmod.R
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
#  along with sybilSBML.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: readSBMLmod
#
#
# The function readSBMLmod() is inspired by the function
# readCbModel() contained in the COBRA Toolbox.
# The algorithm is basically the same.


readSBMLmod <- function(filename, description,
                        def_bnd = SYBIL_SETTINGS("MAXIMUM"),
                        validateSBML = FALSE,
                        extMetFlag = "b",
                        bndCond = TRUE,
                        ignoreNoAn = FALSE,
                        mergeMet = TRUE,
                        balanceReact = TRUE,
                        remUnusedMetReact = TRUE,
                        singletonMet = FALSE,
                        deadEndMet = FALSE,
                        remMet = FALSE,
                        constrMet = FALSE,
                        tol = SYBIL_SETTINGS("TOLERANCE")) {

on.exit(expr = {
    if ( (exists("sbmldoc")) && (!isNULLpointerSBML(sbmldoc)) ) {
        closeSBMLfile(sbmldoc)
    }
} )

#------------------------------------------------------------------------------#
# open the model file

if ( file.exists(filename) == FALSE ) {
   stop("failed to open file ", sQuote(filename))
}

if (missing(description)) {
    mdesc <- filename
}
else {
    mdesc <- description
}


#------------------------------------------------------------------------------#
#                some functions we use to create the model                     #
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# X containes metabolite id's (an object of class "SpeciesReference"). The slot
# "species" contains metabolite id. The function entryforS() gets the array
# index of X from the vector "sybil::met_id(sbml)", which is the line number in
# the stoichiometric matrix S.
#------------------------------------------------------------------------------#

entryforS <- function(X) {

    n    <- length(X[["species"]])  # number of metabolites in X
    si   <- rep(i, n)    # the current column (reaction)
    sj   <- integer(n)   # the row (metabolite)
    s_ji <- numeric(n)   # the stoichiometric coefficient

    remMet <- logical(n) # metabolite removed from the initial metabolite list?

    CURR_MET <- character(n) # metabolites in X

    t <- 0
    for (i in seq(along = X[["species"]])) {
        t <- t + 1

        # This is possible, because the metabolite id's are unique.
        # Keep in mind!
        # The metabolite id's are removed from the metabolites list,
        # but not from the reactions list.

        CURR_MET[t] <- X[["species"]][i]
        if (isTRUE(mergeMet)) {
            met_indCURR <- match(CURR_MET[t], CURR_MET[-t])
        }
        else {
            met_indCURR <- NA
        }

        if (is.na(met_indCURR)) {
            sj[t]     <- match(X[["species"]][i], met_id_tmp)    # the row number
            s_ji[t]   <- X[["stoichiometry"]][i]
            remMet[t] <- ifelse(is.na(sj[t]), FALSE, TRUE)

        }
        else {
            remMet[t] <- FALSE
            s_ji[met_indCURR] <- s_ji[met_indCURR] + X[["stoichiometry"]][i]
            msg <- paste("reaction no.", i, dQuote(react_id_tmp[i]),
                         "metabolite no.", t, dQuote(CURR_MET[t]),
                         "was merged")
            warning(msg, call. = FALSE)
        }

#         if (is.na(sj[t])) {                       # if the current reaction is an
#             return(FALSE)                         # exchange reaction, sj[t] will
#         }                                         # be NA. So we'll leave s_ij[t]
#         else {                                    # at zero.
#             s_ji[t] <- el@stoichiometry           # the stoichiometric coefficient
#         }
    }

    #metUnq <- unique(sj[remMet])

    return(list(sj = sj[remMet], si = si[remMet], s_ji = s_ji[remMet]))

}


#------------------------------------------------------------------------------#
# Check the absolute value of a reaction bound (vmin, vmax). If it is larger
# than def_bnd, it will be replaced by def_bnd.
#------------------------------------------------------------------------------#

checkupplowbnd <- function(x) {

    if (abs(x) > def_bnd) {
        bound <- def_bnd
        if (x < 0) {
            bound <- bound * -1
        }
    }
    else {
        bound <- x
    }
    return(bound)

}


#------------------------------------------------------------------------------#
# A control function. For every part in the SBML file, the Id's are a mandatory
# argument. Here we check, if they are all brave and there.
#------------------------------------------------------------------------------#

missingId <- function(x) {

  mid <- which(x[["id"]] == "no_id")
  if (length(mid) > 0) {
      warning("id is missing in ", is(x)[1], ": ", paste(mid, collapse = ", "))
  }
  
  return(TRUE)

}


#------------------------------------------------------------------------------#
# beautify SBML id's
#------------------------------------------------------------------------------#

formatSBMLid <- function(idstr) {


    idstr <- gsub("-DASH-",        "-",   idstr, fixed = TRUE)
    idstr <- gsub("_DASH_",        "-",   idstr, fixed = TRUE)
    #idstr <- gsub("_FSLASH_",      "/",   idstr, fixed = TRUE)
    #idstr <- gsub("_BSLASH_",      "\\",  idstr, fixed = TRUE)
    idstr <- gsub("_LPAREN_",      "(",   idstr, fixed = TRUE)
    idstr <- gsub("_RPAREN_",      ")",   idstr, fixed = TRUE)
    idstr <- gsub("_LSQBKT_",      "[",   idstr, fixed = TRUE)
    idstr <- gsub("_RSQBKT_",      "]",   idstr, fixed = TRUE)
    idstr <- gsub("_COMMA_",       ",",   idstr, fixed = TRUE)
    idstr <- gsub("_PERIOD_",      ".",   idstr, fixed = TRUE)
    idstr <- gsub("_APOS_",        "'",   idstr, fixed = TRUE)
    idstr <-  sub( "_e_?$",        "(e)", idstr)   # nicer formatting of exchange reactions
    idstr <- gsub("-",             "_",   idstr, fixed = TRUE)
    #idstr <- gsub("&amp;",         "&",   idstr, fixed = TRUE)
    #idstr <- gsub("&lt;",          "<",   idstr, fixed = TRUE)
    #idstr <- gsub("&gt;",          ">",   idstr, fixed = TRUE)
    #idstr <- gsub("&quot;",        "\"",  idstr, fixed = TRUE)

    return(idstr)
}


#------------------------------------------------------------------------------#
# parse the notes field of the reactions
#------------------------------------------------------------------------------#

parseNotesReact <- function(notes) {

  if (regexpr("html:p", notes, fixed = TRUE) == -1) {
      tag <- "p"
  }
  else {
      tag <- "html:p"
  }

  split <- paste("<", tag, ">", sep = "")
  #split <- "\n"

  fields <- strsplit(notes, split, fixed = TRUE)
 # print(fields)

  start_tag  <- paste("<", tag, ">", sep = "")
  end_tag    <- paste("</", tag, ">", sep = "")
  regex      <- paste("^(?:[\\t]*\\Q", start_tag, "\\E)?", "(.*)", "\\Q", end_tag, "\\E", "(?s).*$", sep = "")
#  regex      <- paste("(.*)", end_tag, "(?s).*$", sep = "")
  #print(regex)

  fields_str <- sub(regex, "\\1", fields[[1]], perl = TRUE)
  #print(fields_str)

  subSyst   <- ""
  gpr       <- ""
  gene_rule <- NA

  for (j in 1:length(fields_str)) {
      if (grepl("GENE[_ ]?ASSOCIATION", fields_str[j])) {
      #if (charmatch("GENE", fields_str[j], nomatch = -1) != -1) {
          gpr <- sub("GENE[_ ]?ASSOCIATION: *", "", fields_str[j])
          gene_rule <- sybil:::.parseBoolean(gpr)
          #print(gene_rule)
      }
      if (charmatch("SUBSYSTEM", fields_str[j], nomatch = -1) != -1) {
          subSyst <- sub("SUBSYSTEM: *", "", fields_str[j])
          subSyst <- sub("^S_", "", subSyst, perl = TRUE)
          subSyst <- gsub("[_]+", " ", subSyst)
          if (nchar(subSyst) == 0) {
              subSyst <- "Exchange"
          }
          #print(subSyst)
      }
  }
  if (!is.list(gene_rule)) {
      gene_rule <- sybil:::.parseBoolean("")
  }

  return(list(sub_system = subSyst, genes = gene_rule$gene, rules = gene_rule$rule, gpr = gpr))

}



#------------------------------------------------------------------------------#
#                         main part of the script                              #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                            reading the model                                 #
#------------------------------------------------------------------------------#

message("reading SBML file ... ", appendLF = FALSE)

sbmldoc <- openSBMLfile(filename)

message("OK")


#------------------------------------------------------------------------------#
#                              check the model                                 #
#------------------------------------------------------------------------------#

if (isTRUE(validateSBML)) {
    message("validating SBML file ... ", appendLF = FALSE)
    check <- validateSBMLdocument(sbmldoc)                   # check for errors

    if (!isTRUE(check)) {

        err <- getSBMLerrors(sbmldoc)
        
        nerr <- getNumErrors(err)
        if (nerr["Errors"] > 0) {
        
            message("found errors, trying to fix ... ", appendLF = FALSE)

            closeSBMLfile(sbmldoc)
            hackedModel <- .uglyHack(filename)
            sbmldoc <- openSBMLfile(hackedModel)
            unlink(hackedModel)
            remove(hackedModel)
            check <- validateSBMLdocument(sbmldoc)
        
            if (!isTRUE(check)) {

                sbmlerr <- getSBMLerrors(sbmldoc)
                nerr <- getNumErrors(sbmlerr)
                if ((nerr["Errors"] > 0) ||
                    (nerr["Fatals"] > 0)) {
                    msg <- paste("FAILED: review file", dQuote(filename),
                                 "carefully, returning sbmlError object")
                    warning(msg)
                    return(sbmlerr)
                }
                else {
                    msg <- paste("found warnings and/or infos concerning SBML,",
                                 "you may want to check them with the command",
                                 sQuote(paste("validateSBMLdocument('", filename, "')", sep = "")), "... ")
                    message(msg, appendLF = FALSE)
                    #printSlot(sbmlerr, "Warnings")
                }
            }

        }
        else if (nerr["Fatals"] > 0) {
            msg <- paste("FAILED review file", dQuote(filename),
                         "carefully, returning sbmlError object")
            warning(msg)
            return(sbmlerr)
        }
        else {
            msg <- paste("found warnings and/or infos concerning SBML,",
                         "you may want to check them with the command",
                         sQuote(paste("validateSBMLdocument('", filename, "')", sep = "")), "... ")
            message(msg, appendLF = FALSE)
        }

    }
    message("OK")
}


#------------------------------------------------------------------------------#
#                         generate modelorg object                             #
#------------------------------------------------------------------------------#

message("getting the model ... ", appendLF = FALSE)

sbmlmod <- getSBMLmodel(sbmldoc)

if (is.null(sbmlmod)) {
    message("FAILED")
    stop("Could not get the SBML model. Run SBML validation by ",
         "validateSBMLdocument('", filename, "').")
}

mid   <- ifelse(length(getSBMLmodId(sbmlmod)) == 0,   filename, getSBMLmodId(sbmlmod))
mname <- ifelse(length(getSBMLmodName(sbmlmod)) == 0, filename, getSBMLmodName(sbmlmod))

mod <- sybil::modelorg(mid, mname)              # S4 object of class modelorg

message("OK")


#------------------------------------------------------------------------------#
#                             model description                                #
#------------------------------------------------------------------------------#

if (mdesc == filename) {
    mdesc  <- sub("\\.xml$", "", basename(filename))
}

sybil::mod_desc(mod) <- mdesc


#------------------------------------------------------------------------------#
#                                   units                                      #
#------------------------------------------------------------------------------#

# I have to do this, but later (2007-08-14)


#------------------------------------------------------------------------------#
#                               compartments                                   #
#------------------------------------------------------------------------------#

compartmentsList        <- getSBMLCompartList(sbmlmod)

if (is.null(compartmentsList)) {
    stop("file '", filename, "' has an empty listOfCompartments section")
}

missingId(compartmentsList)
sybil::mod_compart(mod) <- compartmentsList[["id"]]


#------------------------------------------------------------------------------#
#                           initial reactions list                             #
#------------------------------------------------------------------------------#

reactionsList <- getSBMLReactionsList(sbmlmod)

if (is.null(reactionsList)) {
    stop("file '", filename, "' has an empty listOfReactions section")
}

missingId(reactionsList)
react_id_tmp  <- reactionsList[["id"]]
numreact      <- getSBMLnumReactions(sbmlmod)


#------------------------------------------------------------------------------#
#                         initial metabolites list                             #
#------------------------------------------------------------------------------#

metabolitesList <- getSBMLSpeciesList(sbmlmod)

if (is.null(metabolitesList)) {
    stop("file '", filename, "' has an empty listOfSpecies section")
}

missingId(metabolitesList)
metSpIds        <- metabolitesList[["id"]]
#nummet          <- getSBMLnumSpecies(sbmlmod)

if (isTRUE(bndCond)) {
    metSpBnd <- metabolitesList[["boundaryCondition"]]
    met_id_pos <- !metSpBnd
}
else {
    # regular expression to identify external metabolites
    extMetRegEx <- paste("_", extMetFlag, "$", sep = "")
    met_id_pos  <- grep(extMetRegEx, metSpIds, invert = TRUE)
}

met_id_tmp <- metSpIds[met_id_pos]

# number of metabolites
nummet <- length(met_id_tmp)


#------------------------------------------------------------------------------#
#                            reversibilities                                   #
#------------------------------------------------------------------------------#

react_rev_tmp <- reactionsList[["reversible"]]


#------------------------------------------------------------------------------#
#                           data structures                                    #
#------------------------------------------------------------------------------#


#S <- matrix(0, nummet, numreact)
St <- Matrix::Matrix(0, nrow = nummet, ncol = numreact, sparse = TRUE)

lbnd <- numeric(numreact)        # v min
ubnd <- numeric(numreact)        # v max
ocof <- numeric(numreact)        # objective coefficients


#------------------------------------------------------------------------------#
#                       S matrix and constraints, gpr                          #
#------------------------------------------------------------------------------#

message("creating S and parsing constraints ... ", appendLF = FALSE)

# for the gpr stuff
subSys   <- character(numreact)
genes    <- list(numreact)
rules    <- character(numreact)
gpr      <- character(numreact)
#allGenes <- character(0)

# Only one entry, because if one reaction has a notes
# field, all others are supposed to have one. Otherwise
# the gpr stuff does not make sense.
hasNotes <- FALSE
hasAnnot <- FALSE

for (i in 1 : numreact) {

    # the notes/annotations field
    notes <- reactionsList[["notes"]][i]
    annot <- reactionsList[["annotation"]][i]

    if (nchar(notes) > 0) {

        hasNotes    <- TRUE
        notes_field <- parseNotesReact(notes)
        #print(notes_field)
        subSys[i]   <- notes_field$sub_system   # the reaction's sub system: glykolysis, TCA, ...
        genes[[i]]  <- notes_field$genes        # list of genes
        rules[i]    <- notes_field$rules        # rules
        gpr[i]      <- notes_field$gpr          # original gpr association
        #allGenes    <- c(allGenes, genes[[i]])

    }
    else {

        if (nchar(annot) > 0) {
            hasAnnot    <- TRUE
            pn <- regexpr("Pathway Name: [^<]+", annot, perl = TRUE)
            subSys[i] <- substr(annot, (pn+14), pn + ((attr(pn, "match.length"))-1))
        }

    }


    # Check here if reactants and products lists exist, same for the stoichiometry slot

    # Entries for S -- the reactants
    S_tmp <- entryforS(reactionsList[["reactants"]][[i]])
    #print(S_tmp)
    if (is.list(S_tmp) == TRUE) {
        St[S_tmp$sj, i] <- (S_tmp$s_ji * -1)
        #St[S_tmp$sj, S_tmp$si] <- (S_tmp$s_ji * -1)
    }

# Check here if S_tmp is FALSE. Should only be the case in
# the products slot due to the exclusion of external metabolites.
# In that case, the current reaction must be an exchange reaction.

#    else {
#        print(rsbml::reactions(rsbml::model(Mod))[[i]]@id)
#        print(S_tmp)
#        stop("something is wrong here")
#    }

    # Entries for S -- the products
    S_tmp <- entryforS(reactionsList[["products"]][[i]])
    if (length(S_tmp[["s_ji"]]) > 0) {
        #print(S_tmp)
        if (isTRUE(balanceReact)) {
            nnull <- St[S_tmp$sj, i] == 0
            St[S_tmp$sj, i] <- St[S_tmp$sj, i] + S_tmp$s_ji

            if ( any(nnull == FALSE) ) {
                msg <- paste("reaction no.", i,
                             dQuote(react_id_tmp[i]), sum(!nnull),
                             ngettext(sum(!nnull),
                                      "metabolite was balanced",
                                      "metabolites were balanced:\n\t"),
                             paste(dQuote(met_id_tmp[S_tmp$sj[!nnull]]),
                                   collapse = "\n\t "))
                warning(msg, call. = FALSE)
            }

        }
        else {
            St[S_tmp$sj, i] <- S_tmp$s_ji
        }
    }
#    else {
#        print(rsbml::reactions(rsbml::model(Mod))[[i]]@id)
#        print(S_tmp)
#        stop("something is wrong here")
#    }

    # the constraints
    parm <- reactionsList[["kinetic_law"]][[i]]
    if (is.null(parm)) {
        ubnd[i] <- def_bnd
        if (isTRUE(react_rev_tmp[i])) {
            lbnd[i] <- -1 * def_bnd
        }
        else {
            lbnd[i] <- 0
        }
        ocof[i] <- 0
    }
    else {
        for (j in seq(along = parm[["id"]])) {
            if (parm[["id"]][j] == "LOWER_BOUND") {
                lbnd[i] <- checkupplowbnd(parm[["value"]][j])
            }
            if (parm[["id"]][j] == "UPPER_BOUND") {
                ubnd[i] <- checkupplowbnd(parm[["value"]][j])
            }
            if (parm[["id"]][j] == "OBJECTIVE_COEFFICIENT") {
                ocof[i] <- parm[["value"]][j]
            }
            # flux value?   (sbml file)
            # reduced cost?  (sbml file)
        }
    }

}


# ---------------------------------------------------------------------------- #
# search for unused metabolites and unused reactions

# binary version of stoichiometric matrix
#Stb <- St != 0
Stb <- abs(St) > tol

SKIP_METABOLITE   <- rowSums(Stb) != 0
SKIP_REACTION     <- colSums(Stb) != 0


if (isTRUE(remUnusedMetReact)) {
    did <- "and therefore removed from S:"
}
else {
    did <- "in S:"
}


# ---------------------------------------------------------------------------- #
# empty rows

if ( any(SKIP_METABOLITE == FALSE) ) {
    met_list  <- paste(dQuote(met_id_tmp[!SKIP_METABOLITE]),
                       collapse = "\n\t")
    nmet_list <- sum(!SKIP_METABOLITE)
    msg_part  <- paste("not used in any reaction", did)
    msg <- sprintf(ngettext(nmet_list,
                            "%d metabolite is %s %s",
                            "%d metabolites are %s\n\t%s"),
                   nmet_list, msg_part, met_list)
    warning(msg, call. = FALSE)
}


# ---------------------------------------------------------------------------- #
# empty columns

if ( any(SKIP_REACTION == FALSE) ) {
    react_list  <- paste(dQuote(react_id_tmp[!SKIP_REACTION]),
                       collapse = "\n\t")
    nreact_list <- sum(!SKIP_REACTION)
    msg_part    <- paste("not used", did)
    msg <- sprintf(ngettext(nreact_list,
                            "%d reaction is %s %s",
                            "%d reactions are %s\n\t%s"),
                   nreact_list, msg_part, react_list)
    warning(msg, call. = FALSE)
}

if (!isTRUE(remUnusedMetReact)) {
    SKIP_METABOLITE[!SKIP_METABOLITE] <- TRUE
    SKIP_REACTION[!SKIP_REACTION]     <- TRUE
}


# ---------------------------------------------------------------------------- #
# single metabolites

sing_met   <- rep(NA, nrow(St))
sing_react <- rep(NA, ncol(St))

if (isTRUE(singletonMet)) {

    message("identifying reactions containing single metabolites ... ", appendLF = FALSE)

	singleton <- sybil:::.singletonMetabolite(mat = Stb)

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

            met_list  <- paste(dQuote(met_id_tmp[singleton$smet]),
                               collapse = "\n\t")
            nmet_list <- sum(singleton$smet)
            react_list  <- paste(dQuote(react_id_tmp[singleton$sreact]),
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

            #warning(paste(msgm, msgr, sep = "\n\t "))
            warning(msgm, call. = FALSE)
            warning(msgr, call. = FALSE)

        }
        else {

            met_list  <- paste(dQuote(met_id_tmp[singleton$smet]),
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
        sing_met   <- logical(nrow(St))
        sing_react <- logical(ncol(St))
    }
}


# ---------------------------------------------------------------------------- #
# dead end metabolites

de_met   <- rep(NA, nrow(St))
de_react <- rep(NA, ncol(St))

if (isTRUE(deadEndMet)) {

    message("identifying reactions containing dead end metabolites ... ", appendLF = FALSE)

	demr <- sybil:::.deadEndMetabolite(mat   = St,
					         		   lb    = lbnd,
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

            met_list  <- paste(dQuote(met_id_tmp[demr$dem]),
                               collapse = "\n\t")
            nmet_list <- sum(demr$dem)
            react_list  <- paste(dQuote(react_id_tmp[demr$der]),
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

			met_list  <- paste(dQuote(met_id_tmp[demr$dem]),
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
        de_met   <- logical(nrow(St))
        de_react <- logical(ncol(St))
    }
}


# ---------------------------------------------------------------------------- #
# S

St <- St[SKIP_METABOLITE, , drop = FALSE]
St <- St[ , SKIP_REACTION, drop = FALSE]

sybil::S(mod) <- St
#remove(St)

#sbml@S <- S
#remove(S)

# NNZ <- nonZeroElements(S(sbml))
#
# Sne(sbml) <- NNZ$ne
# Sia(sbml) <- NNZ$ia
# Sja(sbml) <- NNZ$ja
# Sar(sbml) <- NNZ$ar


numreact <- sum(SKIP_REACTION)

sybil::met_num(mod)   <- sum(SKIP_METABOLITE)
sybil::react_num(mod) <- numreact

sybil::met_single(mod)   <- sing_met[SKIP_METABOLITE]
sybil::react_single(mod) <- sing_react[SKIP_REACTION]

sybil::met_de(mod)   <- de_met[SKIP_METABOLITE]
sybil::react_de(mod) <- de_react[SKIP_REACTION]

if (isTRUE(constrMet)) {
    lbnd[sing_react] <- 0
    ubnd[sing_react] <- 0
    lbnd[de_react]   <- 0
    ubnd[de_react]   <- 0
}
else {}


sybil::lowbnd(mod)    <- lbnd[SKIP_REACTION]
sybil::uppbnd(mod)    <- ubnd[SKIP_REACTION]
sybil::obj_coef(mod)  <- ocof[SKIP_REACTION]


message("OK")


#------------------------------------------------------------------------------#
#                        gene to reaction mapping                              #
#------------------------------------------------------------------------------#

if (isTRUE(ignoreNoAn)) {
    sybil::gprRules(mod)   <- character(numreact)
    sybil::genes(mod)      <- vector(mode = "list", length = numreact)
    sybil::gpr(mod)        <- character(numreact)
    sybil::allGenes(mod)   <- character(numreact)
    sybil::rxnGeneMat(mod) <- Matrix::Matrix(FALSE, nrow = numreact, ncol = numreact, sparse = TRUE)
    sybil::subSys(mod)     <- Matrix::Matrix(FALSE, nrow = numreact, ncol = 1, sparse = TRUE)
}
else {

    subSys <- subSys[SKIP_REACTION]
    genes  <- genes[SKIP_REACTION]
    rules  <- rules[SKIP_REACTION]
    gpr    <- gpr[SKIP_REACTION]

    if (isTRUE(hasNotes)) {
        message("GPR mapping ... ", appendLF = FALSE)

        #allGenes <- unique(allGenes)
        #allGenesTMP <- unique(allGenes)
        allGenesTMP <- unique(unlist(genes))
        temp <- nchar(allGenesTMP)
        allGenes <- allGenesTMP[which(temp != 0)]


        rxnGeneMat <- Matrix::Matrix(FALSE,
                                     nrow = numreact,
                                     ncol = length(allGenes),
                                     sparse = TRUE)

        for (i in 1 : numreact) {

            if ( (length(genes[[i]] == 1)) && (genes[[i]] != "") ) {
                geneInd <- match(genes[[i]], allGenes)
                rxnGeneMat[i, geneInd] <- TRUE
    
                for (j in 1 : length(geneInd)) {
                    pat  <- paste("x(", j, ")", sep = "")
                    repl <- paste("x[", geneInd[j], "]", sep = "")
    
                    rules[i] <- gsub(pat, repl, rules[i], fixed = TRUE)
                }
            }
        }

        sybil::genes(mod)      <- genes
        sybil::gpr(mod)        <- gpr
        sybil::allGenes(mod)   <- allGenes
        sybil::gprRules(mod)   <- rules
        sybil::rxnGeneMat(mod) <- rxnGeneMat
        #sybil::subSys(mod)     <- subSys
        sybil::subSys(mod)     <- sybil:::.prepareSubSysMatrix(subSys, numreact)

        #sbml@gprRules <- rules
        #sbml@genes <- genes
        #sbml@gpr <- gpr
        #sbml@allGenes <- allGenes
        #sbml@subSys <- subSys

        message("OK")
    }
    else {
        sybil::rxnGeneMat(mod) <- Matrix::Matrix(NA, nrow = 0, ncol = 0)
        if (isTRUE(hasAnnot)) {
            #subSys(sbml)     <- subSys
            sybil::subSys(mod) <- sybil:::.prepareSubSysMatrix(subSys, numreact)
        }
        else {
            sybil::subSys(mod) <- Matrix::Matrix(FALSE,
                                                 nrow = numreact,
                                                 ncol = 1,
                                                 sparse = TRUE)
        }
    }

}


#------------------------------------------------------------------------------#
#                              reaction id's                                   #
#------------------------------------------------------------------------------#
message("cleaning up ... ", appendLF = FALSE)

react_id_tmp   <- sub( "^R[_]+",        "", react_id_tmp[SKIP_REACTION])   # remove the leading R_
#react_id_tmp   <- gsub("_LPAREN_",      "(",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_RPAREN_",      ")",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_LSQBKT_",      "[",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_RSQBKT_",      "]",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_COMMA_",       ",",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_APOS_",        "'",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_DASH_",        "-",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- sub( "_e_?$",         "(e)", react_id_tmp)   # nicer formatting of exchange reactions
#sybil::react_id(mod) <- gsub("-",      "_",   react_id_tmp, fixed = TRUE)
sybil::react_id(mod) <- formatSBMLid(react_id_tmp)


#------------------------------------------------------------------------------#
#                             reaction names                                   #
#------------------------------------------------------------------------------#

react_name_tmp   <- reactionsList[["name"]][SKIP_REACTION]
react_name_tmp   <- sub( "^R[_]+", "",  react_name_tmp)
react_name_tmp   <- gsub("[_]+",   " ", react_name_tmp)
react_name_tmp   <- sub( "\\s+$",  "",  react_name_tmp, perl = TRUE)
sybil::react_name(mod) <- react_name_tmp


#------------------------------------------------------------------------------#
#                             metabolite id's                                  #
#------------------------------------------------------------------------------#

met_id_tmp     <- sub( "^[MSms]_+", "", met_id_tmp[SKIP_METABOLITE])   # remove the leading M_ or S_
#met_id_tmp     <- gsub( "[_]+",           "_",    met_id_tmp)
met_id_tmp     <- sub( "_([A-Za-z0-9]+)$",     "[\\1]", met_id_tmp)   # put the compartment id into square brackets
sybil::met_id(mod) <- gsub("-",      "_",    met_id_tmp, fixed = TRUE)
#sybil::met_id(mod) <- gsub("[-_]+",  "_",     met_id_tmp)


#------------------------------------------------------------------------------#
#                        metabolite compartments                               #
#------------------------------------------------------------------------------#

met_comp_tmp <- metabolitesList[["compartment"]][met_id_pos][SKIP_METABOLITE]

sybil::met_comp(mod) <- match(met_comp_tmp, sybil::mod_compart(mod))


#------------------------------------------------------------------------------#
#                            metabolite names                                  #
#------------------------------------------------------------------------------#

met_name_tmp <- metabolitesList[["name"]][met_id_pos][SKIP_METABOLITE]
met_name_tmp   <- sub( "^[MS]?[_]+", "",  met_name_tmp)
met_name_tmp   <- gsub("[-_]+",   "-", met_name_tmp)
met_name_tmp   <- sub("-$",       "",  met_name_tmp)
met_name_tmp   <- sub( "\\s+$",   "",  met_name_tmp, perl = TRUE)
sybil::met_name(mod) <- met_name_tmp


#------------------------------------------------------------------------------#
#                             check reversibilities                            #
#------------------------------------------------------------------------------#

# check up with the matlab version
# check the reversibilities
react_rev_tmp <- react_rev_tmp[SKIP_REACTION]
isrev <- which(sybil::lowbnd(mod) < 0 & sybil::uppbnd(mod) > 0)
#print(isrev)
react_rev_tmp[isrev]  <- TRUE
sybil::react_rev(mod) <- react_rev_tmp

message("OK")


#------------------------------------------------------------------------------#
#                               validate the model                             #
#------------------------------------------------------------------------------#

message("validating object ... ", appendLF = FALSE)

check <- validObject(mod, test = TRUE)

if (check != TRUE) {
    msg <- paste("Validity check failed:", check, sep = "\n    ")
    warning(msg)
}

message("OK")


#------------------------------------------------------------------------------#
#                                return the model                              #
#------------------------------------------------------------------------------#

delSBMLmodel(sbmlmod)
closeSBMLfile(sbmldoc)

# Returns sbml, an object of the class modelorg
return(mod)

}
