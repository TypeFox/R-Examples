###################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################
# Build an object of class 'poly' from a character vector
# whose each element codes a monomial; in this vector, the inputs are either coded
# by their names (as they are in the colnames of the dataset, or provided - in
# case of transformed data), or by their number

vect2poly <- function(varnames, monomials) {
    nvar <- length(varnames)
    # calcul du degre max
    degmax <- function(X) length(strsplit(X, "*", fixed = T)[[1]])
    # transformation en matrice du vecteur des monomes pour pouvoir appliquer apply
    mmonomials = matrix(monomials, ncol = 1)
    d = max(apply(mmonomials, 1, degmax))
    P <- list()
    nmonomials <- length(monomials)
    indic <- matrix(0, nrow = nmonomials, ncol = nvar)
    indices <- rep(0, d)
    
    for (ilig in 1:nmonomials) {
        
        indices <- decodlemono(monomials[ilig], d, varnames)
        P[[ilig]] <- indices
        names(P[[ilig]]) <- paste("V", 1:d, sep = "")
        # dans la liste P, peu importe les labels
        
        degre <- 1
        while ((degre <= d) && (indices[degre] > 0)) {
            if (indices[degre] > nvar) {
                stop(paste("vect2poly: variable,", indices[degre], ", in the monomials is greater than the number of variables,", 
                  nvar))
            }
            
            
            indic[ilig, indices[degre]] <- 1
            degre <- degre + 1
        }  #fin while degre
        
    }  #fin ilig
    colnames(indic) <- varnames
    lePoly <- poly(P = P, indic = indic)
    return(lePoly)
}  # end vect2poly
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
decodlemono <- function(unmono, d, varnames) {
    
    # Decode an element of the  monomials vector
    # Internal function
    indices <- rep(0, d)
    
    # oter les blancs de debut et fin
    unmono <- trimSpace(unmono)  # trimSpace: library(seqinr)
    
    cde = ""
    zunmono <- strsplit(unmono, "*", fixed = T)[[1]]
    
    ielt <- 1
    for (iz in 1:length(zunmono)) {
        emono <- trimSpace(zunmono[iz])
        if ((emono != "") && (emono != "*")) 
            {
                # l'elt contient qqche
                cde <- paste(cde, " indices[", ielt, "] <- ")
                type <- grepRaw("[1-9]", unmono)
                if ((length(type) > 0) && (type == 1)) {
                  # le vecteur contient les nos de variables
                  cde <- paste(cde, as.integer(emono), "; ")
                } else {
                  # le vecteur contient les noms des variables
                  novar <- match(emono, varnames)
                  if (is.na(novar)) 
                    stop(paste("Variable", zunmono[iz], "not known"))
                  cde <- paste(cde, novar, "; ")
                }
                ielt <- ielt + 1
            }  # fin elt contient qqche
    }  # fin iz
    
    # evaluer l'expression pour affecter les indices des variables
    eval(parse(text = cde))
    if (length(indices) > d) 
        stop(paste("Monomial", unmono, " is of degree greater than ", d))
    return(indices)
}  # end decodlemono
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build an object of class 'polyX' from raw data and
# a character vector whose each element codes a monomial; in this vector, the
# inputs are either coded by their names (as they are in the colnames of the
# dataset, or provided - in case of transformed data), or by their number

vect2polyX <- function(dataX, monomials) {
    varnames <- colnames(dataX)
    Pindic <- vect2poly(varnames, monomials)
    lePolyX <- expand(Pindic, dataX)
    return(lePolyX)
}  # end vect2polyX
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build an object of class 'polyX' from transformed data and
# a character vector whose each element codes a monomial; in this vector, the
# inputs are either coded by their names (as they are in the colnames of the
# dataset, or provided - in case of transformed data), or by their number
vect2polyXT <- function(varnames, dataXT, monomials) {
    nvar <- length(varnames)
    Pindic <- vect2poly(varnames, monomials)
    lePolyX <- new("polyX", dataX.exp = dataXT, Pindic = Pindic)
    return(lePolyX)
}  # end vect2polyXT
