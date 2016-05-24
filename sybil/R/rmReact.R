#  rmReact.R
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
# Function: rmReact
#
#
# The function rmReact() is inspired by the function
# removeRxns() contained in the COBRA Toolbox.
# The algorithm is the same.


rmReact <- function(model, react, rm_met = TRUE) {

  
#------------------------------------------------------------------------------#
#                           check model and react                              #
#------------------------------------------------------------------------------#

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }


  # check this, propably working wrong
  if (is.na(match(is(react)[1], c("reactId", "numeric", "integer", "character")))) {
      stop("argument react must be numeric, character, or of class reactId. Use checkReactId!")
  }

  # argument react comes from the function checkReactId()
  if (is(react, "reactId")) {
      rmReact <- react_pos(react)
  }
  else {
      checked_react <- checkReactId(model, react)
      #print(is(checked_react))
      if (!is(checked_react, "reactId")) {
          stop("Check your reaction Id's")
      }
      else {
          rmReact <- react_pos(checked_react)
      }
  }


#  if ((is(react, "numeric")) || (is(react, "integer"))) {
#      if (max(react) > react_num(model)) {
#          stop("check your reaction Id's!")
#      }
#      else {
#          rmReact <- react
#      }
#  }

#  if (is(react, "character")) {
#      checked_react <- checkReactId(model, react)
#      if ((is.logical(checked_react)) && (checked_react == FALSE)) {
#          stop("check your reaction Id's!")
#      }
#      else {
#          rmReact <- react_pos(checked_react)
#      }
#  }

  
#------------------------------------------------------------------------------#
#                              logical vector                                  #
#------------------------------------------------------------------------------#

  # logical vector: length = number of reactions, positions to be deleted
  # are set to FALSE
  #(numeric(react_num(model)) == 0)
  keepReact          <- rep(react_num(model), x=TRUE)
  keepReact[rmReact] <- FALSE


#------------------------------------------------------------------------------#
#                           construct the new model                            #
#------------------------------------------------------------------------------#
  
  if (is(model, "modelorg_irrev")) {
      mod_out <- modelorg_irrev(mod_id(model), mod_name(model))

      irrev(mod_out) <- TRUE

      for (i in 1:length(rmReact)) {
          remInd <- rmReact[i]
          if (matchrev(model)[remInd] > 0) {
              revInd <- matchrev(model)[remInd]
              react_rev(model)[revInd] <- FALSE

              #react_id(model)[revInd] <- substr(react_id(model)[revInd], 1, (nchar(react_id(model)[revInd]) - 2))
              react_id(model)[revInd] <- sub("_[bfr]$", "", react_id(model)[revInd])
           }
      }
  }
  else {
      mod_out <- modelorg(mod_id(model), mod_name(model))
  }
      
  # model description
  mod_desc(mod_out)     <- mod_desc(model)
  mod_compart(mod_out)  <- mod_compart(model)
  
  # reactions
  react_rev(mod_out)    <- react_rev(model)[keepReact]
  react_id(mod_out)     <- react_id(model)[keepReact]
  react_name(mod_out)   <- react_name(model)[keepReact]
  lowbnd(mod_out)       <- lowbnd(model)[keepReact]
  uppbnd(mod_out)       <- uppbnd(model)[keepReact]
  obj_coef(mod_out)     <- obj_coef(model)[keepReact]
  react_single(mod_out) <- react_single(model)[keepReact]
  react_de(mod_out)     <- react_de(model)[keepReact]

  react_num(mod_out)    <- length(react_id(mod_out))

  # stoichiometric matrix
  S(mod_out)            <- S(model)[ , keepReact, drop = FALSE]

  # GPR stuff
  if (!is.na(match("genes", slotNames(model)))) {

      gprRules(mod_out)     <- gprRules(model)[keepReact]
      genes(mod_out)        <- genes(model)[keepReact]
      gpr(mod_out)          <- gpr(model)[keepReact]
      rxnGeneMat(mod_out)   <- rxnGeneMat(model)[keepReact, , drop = FALSE]
      subSys(mod_out)       <- subSys(model)[keepReact, , drop = FALSE]

      ag                    <- unique(unlist(genes(mod_out)))
      # old code tried to assign NULL to allGenes, if no gene was left.
      if(length(ag)==0){
      	allGenes(mod_out)   <- character(0)
      }
      else {
      	ncag                <- nchar(ag)
     	ag					<- ag[which(ncag != 0)]
      	allGenes(mod_out)   <- ag
      }
      

      # reaction to gene mapping
      #SrGMbin     <- rxnGeneMat(mod_out) != 0

      #SrGMbindiag <- diag(crossprod(SrGMbin))

      #keepGenes   <- ifelse(SrGMbindiag == 0, FALSE, TRUE)
      keepGenes <- sapply(allGenes(model), function(x) match(x, allGenes(mod_out)))
      keepGenes <- ifelse(is.na(keepGenes), FALSE, TRUE)
      #print(keepGenes)

      rxnGeneMat(mod_out)   <- rxnGeneMat(mod_out)[, keepGenes, drop = FALSE]
      #print(dim(rxnGeneMat))
  }
  

  # check what to do with the reversible match stuff (model in irrev format)  
  if (is(model, "modelorg_irrev")) {

      matchrev(mod_out)     <- .reassignFwBwMatch(matchrev(model),
                                                          keepReact)

      # propably we do not need this:
      react_rev(mod_out)[matchrev(mod_out) == 0] <- FALSE

      #rev2irrev(mod_out)    <- integer(irrev_react_num)
      #irrev2rev(mod_out)    <- integer(irrev_react_num)
  }


#------------------------------------------------------------------------------#
#                         removed unused metabolites                           #
#------------------------------------------------------------------------------#

  if (rm_met == TRUE) {

      # binary matrix of S
      Sbin     <- S(mod_out) != 0

      Sbindiag <- diag(tcrossprod(Sbin))

      keepMet  <- ifelse(Sbindiag == 0, FALSE, TRUE)

      S(mod_out)          <- S(mod_out)[keepMet, ]
      met_num(mod_out)    <- sum(keepMet == TRUE)
      met_id(mod_out)     <- met_id(model)[keepMet]
      met_name(mod_out)   <- met_name(model)[keepMet]
      met_comp(mod_out)   <- met_comp(model)[keepMet]
      met_single(mod_out) <- met_single(model)[keepMet]
      met_de(mod_out)     <- met_de(model)[keepMet]
  }
  else {
      met_num(mod_out)  <- met_num(model)
      met_id(mod_out)   <- met_id(model)
      met_comp(mod_out) <- met_comp(model)
  }

  check <- validObject(mod_out, test = TRUE)

  if (check != TRUE) {
      warning(paste("Validity check failed:", check, sep = "\n    "), call. = FALSE)
  }

  return(mod_out)

}

