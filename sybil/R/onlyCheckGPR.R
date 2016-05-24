#  onlyCheckGPR.R
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
# Function: onlyCheckGPR
#
#
#
#
#     model must be an object of modelorg
#     gprRules stands for the logical Expressions (GPR Rules)
#     reactNr must be numeric

onlyCheckGPR <- function(model,gprRules,reactNr,verboseMode = 1) {

  # is the chosen model type of modelorg?
  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (!missing(reactNr)) {
      if (!is(reactNr, "numeric")) {
          stop("argument reactNr must be numeric!")
      }
  }

  ###### is the logical Expression correct? ######

  # allowed logical symbols are stored in logSymb
  #logSymb <- c("and", "or", "not", "&", "|", "!")
  logSymb <- c("and", "or", "not", "AND", "OR", "NOT", "&", "|", "!")

  ### check every expression and change..
  for (anz in seq(1, length(gprRules), by = 1)) {

      if (verboseMode > 1 && !missing(reactNr)) { 
      cat(paste("Checking logical expression for reaction",
      react_id(model)[reactNr[anz]],"... "), sep="" )  }
      #check <- tolower(gprRules[anz])
      #check <- gsub("[()]","", check)
      check <- gsub("[()]","", gprRules[anz])
      check <- strsplit(check, "\\s+")
      check <- unlist(check)

      if (check[1] == "") {	
          check <- check[-1]
      }

      #print(check)

      # now start the checking
      if ( identical((length(check) %% 2), 0) ) {

          if (verboseMode > 0) {stop("=> Wrong logical expression!\n") }
          #errors <- TRUE
	  #changeVector <- c(changeVector, FALSE)
          #next

      }
      # only one gene?
      if ( identical(length(check),as.integer(1)) ) {

          if (!any(!is.na(match(logSymb, check[1])))) {

	      # do nothing here if the logical expression contains a Gene. Else stop
	      if ( any(is.na(match(check[1], allGenes(model)))) ) {

	          if (verboseMode > 0) {stop(paste("=> Gene ", check[1] ,
                  " does not exist in the model!\n", sep="")) }
	          #errors <- TRUE
		  #changeVector <- c(changeVector, FALSE)
	          #next

	      }

          } 
          else {

              if (verboseMode > 0) {stop("=> Wrong logical expression!\n") }
	      #errors <- TRUE
	      #changeVector <- c(changeVector, FALSE)
	      #next

          }

      }

      #print(check) 
  
      # go here if there is more than one splitted string. Check the logical expression
      else {

          x <- FALSE

          for (i in seq(2, length(check), by = 2)) {

	      if (any(!is.na(match(logSymb, check[i-1])))) {

	          if (verboseMode > 0) {stop("=> Wrong logical expression!\n") }
	          #errors <- TRUE
	          #x <- TRUE
	          #next

              }

              if ( any(is.na(match(check[i-1], allGenes(model)))) ) {

	          if (verboseMode > 0) { stop(paste("=> Gene ", check[i-1] ,
                  " does not exist in the model!\n", sep="")) }
	          #errors <- TRUE
	          #x <- TRUE
	          #next

	      }

	      if ((i+1) <= length(check)) {

	          if (any(!is.na(match(logSymb, check[i+1])))) {

	              if (verboseMode > 0) { stop("=> Wrong logical expression!\n") }
	 	      #errors <- TRUE
		      #x <- TRUE
		      #next

                  }

	          if ( any(is.na(match(check[i+1], allGenes(model)))) ) {

	              if (verboseMode > 0) { stop(paste("=> Gene ", check[i+1] ,
                      " does not exist in the model!\n", sep="")) }
		      #errors <- TRUE
		      #x <- TRUE
	              #next

	          }

	      }
	
	      if (!any(!is.na(match(logSymb, check[i])))) {

	          if (verboseMode > 0) { stop("=> Wrong logical expression!\n") }
	          #errors <- TRUE
	          #x <- TRUE
	          #next

              }      
	
          }

          # next, if an error occured
          #if (x) {
	  #    changeVector <- c(changeVector, FALSE)
          #    next
          #}	

      }

  #changeVector <- c(changeVector, TRUE)
  #if ( (verboseMode > 1) && !missing(reactNr) ) { cat("OK. \n") }

  }

  #if (errors) { cat("Check logical expression(s)! \n") }

  #if (!errors && (verboseMode > 0)) { cat("OK. \n") }

  #return(changeVector)

}
