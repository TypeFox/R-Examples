#  geneDel.R
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
# Function: geneDel
#
#
# The function geneDel() is inspired by the function
# deleteModelGenes() contained in the COBRA Toolbox.


geneDel <- function(model, genes, checkId = FALSE) {
#geneDel <- function(model, genes, lpmodel, solver = SYBIL_SETTINGS("SOLVER")) {

  #if (missing(lpmodel)) {
  #    solver = "none"
  #}

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }
  
  if (isTRUE(checkId)) {
      if (is(genes, "character")) {
          # Check if all genes are there
          geneExist <- which(is.na(match(genes, allGenes(model))))
        
          # if thats not the case ...
          if (length(geneExist) != 0) {
              stop(sprintf(ngettext(length(geneExist),
                                    "gene %s does not exist",
                                    "genes %s do not exist"
                                    ), paste(sQuote(genes[geneExist]), collapse = ", ")
                           )
                   )
          }
        
          geneInd <- match(genes, allGenes(model))
      }
      else {
          if (max(genes) > length(allGenes(model))) {
              stop("indices in argument genes do not exist")
          }
          else {
              geneInd <- genes
          }
      }
  }
  else {
      if (is(genes, "numeric")) {
          geneInd <- genes
      }
      else {
          stop("argument genes must be numeric, if checkId is FALSE")
      }
  }

  # Get the reaction id's of reactions in the gene list. Returns a list
  # Perhaps we find something better here.
#  reactInd <- as.list(apply(
#                    as.matrix(rxnGeneMat(model)[,geneInd]), 2, function(x)
#                                                        which(x != 0)
#
#                    ))
  reactInd <- apply(rxnGeneMat(model)[,geneInd, drop = FALSE], 2, function(x) which(x != 0) )

                    
 #   print(reactInd)  
  #print(unlist(reactInd))

  reactInd <- unlist(reactInd)
#return(reactInd)  

#print(reactInd)
  
  #x <- logical(length(allGenes(model)))
  x <- rep(TRUE, length(allGenes(model)))
  #print(x)
  x[geneInd] <- FALSE
  constReact <- logical(length(reactInd))
#print(constReact)


  # also do better here

  # Constrain a reaction if the corresponding gpr rule is FALSE.
  # If that's the case, the reaction needs gene bla.

  ru <- gprRules(model)[reactInd]
  for(i in 1:length(reactInd)) {
      #print(reactInd[i])
      #print(ru[i])
      #ev <- eval(parse(text = ru[i]))
      ev <- tryCatch(eval(parse(text = ru[i])), error = function(e) e)
      if (is(ev, "simpleError")) {
          stop("wrong gene association:",
               "\nreaction no. ", reactInd[i],
               "\nreaction id: ", react_id(model)[reactInd[i]],
               "\ngpr: ", gpr(model)[reactInd[i]])
      }
      if (any(is.na(ev))) {
          warning("reference to non existing gene id in gene association:",
                  "\nreaction no. ", reactInd[i],
                  "\nreaction id: ", react_id(model)[reactInd[i]],
                  "\nignoring gpr ", sQuote(gpr(model)[reactInd[i]]))
      }
      else {
          #if (eval(parse(text = gprRules(model)[reactInd[i]])) == FALSE) {
		  #if (eval(parse(text = ru[i])) == FALSE) {
		  if (ev == FALSE) {
			  #print(reactInd[i])
			  #print("cool")
			  constReact[i] <- TRUE
		  }
      }
  }

#  print(reactInd)
#  print(constReact)
#  print(reactInd[constReact])
  
  if (any(constReact)) {
     
      #if (solver == "none") {
      #    model <- changeBounds(model, reactInd[constReact])
      #
      #}

      #if (solver == "glpkAPI") {
      #     model <- changeBounds(lpmodel, reactInd[constReact], solver = solver)
      #}

      #return(model)

    return(unique(reactInd[constReact]))
  }

  #return(as.numeric(NA))
  return(NULL)

}
