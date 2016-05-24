#  onlyChangeGPR.R
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
# Function: onlyChangeGPR
#
#
#
#
#     model must be an object of modelorg
#     gprRules stands for the new logical Expressions (GPR Rules)
#     reactNr must be numeric

onlyChangeGPR <- function(model,gprRules,reactNr,verboseMode = 0) {

  # is the chosen model type of modelorg?
  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (!missing(reactNr)) {
      if (!is(reactNr, "numeric")) {
          stop("argument reactNr must be numeric!")
      }
  }

  # allowed logical symbols are stored in logSymb
  logSymb <- c("and", "or", "not", "AND", "OR", "NOT", "&", "|", "!")

  #if (verboseMode > 1) { cat("Updating ... ") }

  for (anz in seq(1, length(gprRules), by = 1)) {

      #check <- tolower(gprRules[anz])
      #check <- gsub("[()]","", check)
      check <- gsub("[()]","", gprRules[anz])
      check <- strsplit(check, "\\s+")
      check <- unlist(check)

      if (check[1] == "") {	
          check <- check[-1]
      }

      # update gpr()

      gpr(model)[reactNr[anz]] <- gprRules[anz]
      #print(gprRules)

      # update genes(), rxnGeneMat() and gprRules()

      genes(model)[[reactNr[anz]]] <- ""
      rxnGeneMat(model)[reactNr[anz],] <- 0
      gprRules(model)[reactNr[anz]] <- .parseBoolean(gprRules[anz])[[2]]
      gprRules_tmp <- as.numeric(unlist(strsplit(gprRules(model)[reactNr[anz]], "\\D+", perl = TRUE))[-1])
      #print(gprRules_tmp)
      #print(gprRules(model)[reactNr[anz]])
      #print(gprRules(model)[reactNr])
      #print(.parseBoolean(gprRules)[[2]])

      j <- 1

      for (i in seq(1, length(check), by = 2)) {

          genes(model)[[reactNr[anz]]][j] <- check[i]
          rxnGeneMat(model)[reactNr[anz],match(check[i], allGenes(model))] <- 1
          #print(gprRules(model)[reactNr[anz]])
#           gprRules(model)[reactNr[anz]] <- gsub( paste("\\(", deparse(j), "\\)", sep = ""),
#                                                  paste("[", deparse(as.double(match(check[i], allGenes(model)))),"]", sep = ""),
#                                                  gprRules(model)[reactNr[anz]], fixed = TRUE
#                                                )
          gprRules(model)[reactNr[anz]] <- gsub( paste("\\(", gprRules_tmp[j], "\\)", sep = ""),
                                                 paste("[", match(check[i], allGenes(model)),"]", sep = ""),
                                                 gprRules(model)[reactNr[anz]]
                                               )
          #print(gprRules(model)[reactNr[anz]])
          #print(paste("\\(", deparse(j), "\\)", sep = ""))
          #print(paste("[", match(check[i], allGenes(model)),"]", sep = ""))
          #print(check[i])
   	      j <- j + 1

      }
      
      genes(model)[[reactNr[anz]]] <- unique(genes(model)[[reactNr[anz]]])
            
      # edit the expression to make it look better ;)..
      # ..add blanks; add brackets around the expression

      if (length(check) > 1) {

          gpr(model)[reactNr[anz]] <- gsub(" ","  ",gpr(model)[reactNr[anz]]) 
          gprRules(model)[reactNr[anz]] <- gsub(" ","  ",gprRules(model)[reactNr[anz]]) 
          gpr(model)[reactNr[anz]] <- paste("(",gpr(model)[reactNr[anz]],")")
          gprRules(model)[reactNr[anz]] <- paste("(",gprRules(model)[reactNr[anz]],")")

      }

      #if (verboseMode > 1) { cat("OK. \n") }

  }

  return(model)

}
