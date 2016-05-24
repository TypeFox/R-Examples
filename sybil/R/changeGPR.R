#  changeGPR.R
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
# Function: changeGPR
#
#
#
#
#     model must be an object of modelorg
#     react stands for the reaction, which will be changed
#     (Instance of checkReactId, positions or names)
#     gprRules stands for the new logical Expressions (GPR Rules)

changeGPR <- function(model, react, gprRules = "logicalExpression", verboseMode = 1) {

  # is the chosen model type of modelorg?
  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }
  
  # is the gprRules field empty?
  if (missing(gprRules)) {
      stop("please input an expression!")
  }

  if (any(!is.na(match(gprRules,"")))) {
      stop("\"\" is no valid logical expression!")
  }

  # check the reaction
  if (is(react, "reactId")) {

      reactNr <- react_pos(react)

  }
  else {

      if (!is(checkReactId(model, react), "reactId")) {

          stop("At least one reaction does not exist in the model!")

      }
      else {

          reactNr <- react_pos(checkReactId(model, react))

      }

  }

  if ( !identical(length(reactNr), length(gprRules)) ) {
      stop("not as many logical expressions as reactions!")
  }

  ###### is the logical Expression correct? ######

  checkV <- onlyCheckGPR(model, gprRules, reactNr, verboseMode=verboseMode)

  if (verboseMode > 1) { cat("Updating ... ") }

  ###### the expression is correct. Now change ######

  model <- onlyChangeGPR(model, gprRules, reactNr, verboseMode=verboseMode)

  #for (anz in seq(1, length(checkV), by = 1)) { 
  #    if (checkV[anz]) { model <- onlyChangeGPR(model, gprRules[anz],
  #                        reactNr[anz], verboseMode=verboseMode) }
  #}

  if (verboseMode > 1) { cat("OK.\n") }

  if (verboseMode > 0) { cat("Finished.\n") }


  return(model)

}
