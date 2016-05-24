#  mod2irrev.R
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
# Function: mod2irrev
#
#
# The function mod2irrev() is inspired by the function
# convertToIrreversible() contained in the COBRA Toolbox.
# The algorithm is the same.


mod2irrev <- function(model, exex = FALSE) {

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (is(model, "modelorg_irrev")) {
      return(model)
  }

  # set the reversible flag of the exchange reactions to FALSE
  if (isTRUE(exex)) {
      exch <- findExchReact(model)
      ex <- react_pos(exch)
      react_rev(model)[ex] <- FALSE
  }


#------------------------------------------------------------------------------#
#                           data structures                                    #
#------------------------------------------------------------------------------#

  modelIr <- modelorg_irrev(mod_id(model), mod_name(model))

  irrev(modelIr)        <- TRUE

  # number of reactions in the irrev model
  irrev_react_num      <- react_num(model) + sum(react_rev(model) == TRUE)


  # data structures not different from the rev model
  mod_name(modelIr)     <- mod_name(model)
  mod_id(modelIr)       <- mod_id(model)
  mod_compart(modelIr)  <- mod_compart(model)
  met_id(modelIr)       <- met_id(model)
  met_num(modelIr)      <- met_num(model)
  met_name(modelIr)     <- met_name(model)
  met_comp(modelIr)     <- met_comp(model)
  
  met_single(modelIr)   <- met_single(model)
  met_de(modelIr)       <- met_de(model)
  
  mod_desc(modelIr)     <- paste(mod_desc(model), "irreversible")
  #if (!is.na(match("description", slotNames(model)))) {
  #    description(modelIr) <- paste(description(model), "irreversible")
  #}

  # data structures different from the rev model
  react_num(modelIr)    <- irrev_react_num
  react_rev(modelIr)    <- logical(irrev_react_num)
  react_id(modelIr)     <- character(irrev_react_num)
  react_name(modelIr)   <- character(irrev_react_num)
  lowbnd(modelIr)       <- numeric(irrev_react_num)
  uppbnd(modelIr)       <- numeric(irrev_react_num)
  obj_coef(modelIr)     <- numeric(irrev_react_num)

  matchrev(modelIr)     <- integer(irrev_react_num)
  rev2irrev(modelIr)    <- matrix(0, nrow = react_num(model), ncol = 2)
  irrev2rev(modelIr)    <- integer(irrev_react_num)
  #irrev2rev             <- integer(irrev_react_num)

  Sr <- S(model)
  Si <- Matrix::Matrix(0,
                       nrow = met_num(model),
                       ncol = irrev_react_num,
                       sparse = TRUE)
  #print(dim(Si))

#------------------------------------------------------------------------------#
#                          build the new model                                 #
#------------------------------------------------------------------------------#

  # counter for the reactions in the irrev model

  counti <- 0

  message("building the new model ...")
  progr <- .progressBar() # initialize the progressbar

  for (i in 1 : react_num(model)) {

      progr <- .progressBar(i, react_num(model), progr)

      counti <- counti + 1

      react_rev(modelIr)[counti] <- react_rev(model)[i]
      irrev2rev(modelIr)[counti] <- i

      # Reaction in negative direction in the rev model
      if ((uppbnd(model)[i] <= 0) && (lowbnd(model)[i] < 0)) {

          # turn the reaction in positive direction

          # the constraints
          #uppbnd(modelIr)[counti]     <- (uppbnd(model)[i] * -1)
          #lowbnd(modelIr)[counti]     <- (lowbnd(model)[i] * -1)
          uppbnd(modelIr)[counti]     <- (lowbnd(model)[i] * -1)
          lowbnd(modelIr)[counti]     <- (uppbnd(model)[i] * -1)

          # the stiochiometric coefficients
          Si[,counti]                <- (Sr[,i] * -1)
          #Si[,counti]               <- (S(model)[,i] * -1)

          # the objective coefficient
          obj_coef(modelIr)[counti]   <- (obj_coef(model)[i] * -1)

          # reaction name and id
          react_name(modelIr)[counti] <- react_name(model)[i]
          react_id(modelIr)[counti]   <- paste(react_id(model)[i], "_r", sep = "")
          react_rev(modelIr)[counti]  <- FALSE
          react_single(modelIr)[counti] <- react_single(model)[i]
          react_de(modelIr)[counti]     <- react_de(model)[i]

          # set reaction in rev model to irreversible
          react_rev(model)[i]       <- FALSE
      }

      # Reaction is in positive direction in rev model
      else {


          # the contraints
          uppbnd(modelIr)[counti]     <- uppbnd(model)[i]
          if (lowbnd(model)[i] < 0) {
              lowbnd(modelIr)[counti]     <- 0
          }
          else {
              lowbnd(modelIr)[counti]     <- lowbnd(model)[i]
          }

          # the stiochiometric coefficients
          Si[,counti]                <- Sr[,i]
          #Si[,counti]               <- S(model)[,i]

          # the objective coefficient
          obj_coef(modelIr)[counti]   <- obj_coef(model)[i]

          # reaction name and id
          react_name(modelIr)[counti] <- react_name(model)[i]
          react_id(modelIr)[counti]   <- react_id(model)[i]
          #react_rev(modelIr)[counti]  <- FALSE
          react_single(modelIr)[counti] <- react_single(model)[i]
          react_de(modelIr)[counti]     <- react_de(model)[i]
      }

      # current reaction is reversible
      if (react_rev(model)[i] == TRUE) {

          # ++ irreversible counter
          counti <- counti + 1

          # gedoens
          matchrev(modelIr)[counti]   <- as.integer(counti - 1)
          matchrev(modelIr)[counti-1] <- as.integer(counti)
          rev2irrev(modelIr)[i,]     <- c((counti - 1), counti)
          irrev2rev(modelIr)[counti]  <- i#as.integer(counti)
          #irrev2rev[counti]          <- i

          # reaction name
          react_name(modelIr)[(counti-1):counti] <- react_name(model)[i]

          # reaction id with forward, reverse indicator
          react_id(modelIr)[counti-1] <- paste(react_id(model)[i], "_f", sep = "")
          react_id(modelIr)[counti]   <- paste(react_id(model)[i], "_b", sep = "")

          # react_single/react_de
          react_single(modelIr)[counti]   <- react_single(model)[i]
          react_single(modelIr)[counti-1] <- react_single(model)[i]
          react_de(modelIr)[counti]       <- react_de(model)[i]
          react_de(modelIr)[counti-1]     <- react_de(model)[i]

          # stoichiometric coefficients
          Si[,counti]               <- (Sr[,i] * -1)
          #Si[,counti]              <- (S(model)[,i] * -1)

          # reverse inticator
          react_rev(modelIr)[counti]  <- TRUE
          #react_rev(modelIr)[i]       <- TRUE

          # lower bound
          lowbnd(modelIr)[counti]     <- 0

          # upper bound
          uppbnd(modelIr)[counti]     <- (lowbnd(model)[i] * -1)

          # objective coefficient
          obj_coef(modelIr)[counti]   <- 0
      }
      else {
          matchrev(modelIr)[counti]   <- as.integer(0)
          rev2irrev(modelIr)[i,]      <- as.integer(c(counti, counti))
          react_single(modelIr)[counti] <- react_single(model)[i]
      }

  }

  message("cleaning up ...")

  if (length(subSys(model)) > 0) {
      subSys(modelIr)   <- subSys(model)[irrev2rev(modelIr), , drop = FALSE]
      #subSys(modelIr)   <- subSys(model)[irrev2rev]
  }

  if (length(genes(model)) > 0) {
      genes(modelIr)      <- genes(model)[irrev2rev(modelIr)]
      gprRules(modelIr)   <- gprRules(model)[irrev2rev(modelIr)]
      gpr(modelIr)        <- gpr(model)[irrev2rev(modelIr)]
      #genes(modelIr)      <- genes(model)[irrev2rev]
      #gprRules(modelIr)   <- gprRules(model)[irrev2rev]
      #gpr(modelIr)        <- gpr(model)[irrev2rev]

      allGenes(modelIr)   <- allGenes(model)

      rxnG_temp          <- rxnGeneMat(model)
      rxnG_temp          <- rxnG_temp[irrev2rev(modelIr),]
      #rxnG_temp          <- rxnG_temp[irrev2rev,]
      rxnGeneMat(modelIr) <- rxnG_temp

      #modeli@rxnGeneMat <- model@rxnGeneMat[modeli@irrev2rev,]
  }


  if (irrev_react_num > counti) {
      react_num(modelIr)    <- as.integer(counti)
      react_rev(modelIr)    <- react_rev(modelIr)[1:counti]
      react_id(modelIr)     <- react_id(modelIr)[1:counti]
      react_name(modelIr)   <- react_name(modelIr)[1:counti]
      lowbnd(modelIr)       <- lowbnd(modelIr)[1:counti]
      uppbnd(modelIr)       <- uppbnd(modelIr)[1:counti]
      obj_coef(modelIr)     <- obj_coef(modelIr)[1:counti]

      matchrev(modelIr)     <- matchrev(modelIr)[1:counti]
      irrev2rev(modelIr)    <- irrev2rev(modelIr)[1:counti]
      Si                    <- Si[,1:counti]
  }

  S(modelIr) <- Si

  # reset the reversible flag of the exchange reactions to TRUE
  if (exex == TRUE) {
      ex_new <- which(irrev2rev(modelIr) %in% ex)
      react_rev(modelIr)[ex_new] <- TRUE
  }


  check <- validObject(modelIr, test = TRUE)

  if (check != TRUE) {
      warning(paste("Validity check failed:", check, sep = "\n    "), call. = FALSE)
  }

  return(modelIr)

}
