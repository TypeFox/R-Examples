#  doubleReact.R
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
# Function: doubleReact
#
# Returns a list with identical reactions:
# The every list element contains two entries:
#     1) the number of metabolites of the reactions,
#     2) an integer vector containing the reaction id's.
#
# So the final result is grouped in classes of reactions;
# a class is defined by the number of corresponding
# metabolites.


doubleReact <- function(model, checkRev = TRUE, linInd = FALSE) {

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  # test, if there are any duplicated columns
#  dup <- duplicated(as.matrix(S(model)), MARGIN = 2)
#  if (sum(dup) == 0) {
#      return(FALSE)
#  }
#  else {
#      # mat:   stoichiometric matrix
#      # nnzpc: number of non-zeros per column
#      # remc:  reaction types (in terms of number of metabolites)
#      mat   <- as(S(model), "CsparseMatrix")
#      nnzpc <- mat@p[-1] - mat@p[-length(mat@p)]
#      remc  <- sort(unique(nnzpc[dup]))
#  }
#  
#  print(sum(dup))
#  #print(nnzpc)
#  print(remc)
#  print(nnzpc[dup])
#
#  return(shrinkMatrix(model, j = (1:react_num(model))[dup]))
#
#  return(TRUE)
#------------------------------------------------------------------------------#
#                        reaction adjacency matrix                             #
#------------------------------------------------------------------------------#

  # binary matrix of S
  Sbin <- S(model) != 0

  # diagonal of reaction adjacency matrix, number of compounds in each reaction
  #Sbindiag <- diag(t(Sbin) %*% Sbin)
  Sbindiag <- diag(crossprod(Sbin))


#------------------------------------------------------------------------------#
# Next, we compute a vector containing all different numbers of compunds in the
# reactions, e.g.:
#     Sbindiag <- c(3, 6, 3, 4, 4, 3)
# so recatlength will be
#     reactlength <- c(3, 4, 6)
#------------------------------------------------------------------------------#
  
  # vector with all different reaction length
  # (length in the meaning of "number of metabolites")
  reactlength <- sort(unique(Sbindiag))


  # columns in S corresponding to reactlength
  # get the columns of S having a unique reactlength[i] <- j
  columns <- sapply(reactlength, function(x) which(x == Sbindiag))

  
#------------------------------------------------------------------------------#
#                         compare the columns of S                             #
#------------------------------------------------------------------------------#

  double_tmp  <- integer(0)  # a vector where all identical reactions are stored
  double      <- list()      # the main list
  clm         <- integer(1)  # counter for the main list (reactions with m metabolites)
  clr         <- integer(1)  # counter for the sub lists (identical reactions)
  checked     <- logical(1)
  
  clm <- 0

  for (i in 1:length(reactlength)) {

      # number of reactions containing reactlength[i] metabolites
      numr <- length(columns[[i]])
    
      # only if columns[[i]] contains more than one recation (column index)
      if (numr > 1) {
          #print("next")
          message(
                  sprintf(ngettext(reactlength[i],
                          "processing %s reactions with one metabolite ...",
                          "processing %s reactions with %s metabolites ..."
                                   ), numr, reactlength[i]
                          )
                  )

          # Put the columns of S into a new matrix, which contain the same
          # number of metabolites
          Stmp <- cBind(S(model)[, columns[[i]] ])

          # translate Stmp into a matrix with the row indices of the
          # non zero elements
          Stmp_rows <- matrix(
                              apply(
                                    Stmp, 2, function(x)
                                                 which(x != 0, arr.ind = TRUE)
                                    )
                              , nrow = reactlength[i]
                              )
          #print(Stmp_rows)
          
          # A (temporary) list for identical reactions with
          # reactlength[i] metabolites
          double_react <- list()
          #clr <- 1
          clr <- 0

          # walk through Stmp
          for (k in 1:(numr-1)) {

              .progressDots(10, k, (numr-1))
            
              # a vector containing the column indices of S of reactions
              # that are identical to reaction k in Stmp
              ident_react <- integer(0)

              l <- k + 1
              while(l <= numr) {
                  checked <- FALSE

                  # If both columns were already compared with some
                  # other column, they are the same and need not to
                  # be checked again. Uses double_tmp.
                  if (all( match(
                                 c(columns[[i]][k], columns[[i]][l]),
                                 double_tmp, nomatch = 0
                                 )
                          != 0)
                      ) {

                      l <- l + 1
                      next
                  }

                  # First, we only check whether the row indices of the
                  # non zero entries in columns k and l are the same (the
                  # columns are much smaller if we only take the row indices).
                  if (identical(Stmp_rows[, k],Stmp_rows[, l])) {

                      # If they are the same, the stiochiometric coeficients
                      # are checked. If they are also the same, the two columns
                      # k and l (reactions) are considered to be the same.

                      if (isTRUE(linInd)) {
                          # We only chack for linear independence, if all
                          # stoichiometric coefficients heve the same sign,
                          # otherwise (-1, 1) and (1, -1) are the same.
                          a <- Stmp[Stmp_rows[, k], k]
                          b <- Stmp[Stmp_rows[, l], l]
                          if (identical(sign(a), sign(b))) {
                              decomp     <- qr(cbind(Stmp[Stmp_rows[, k], k], Stmp[Stmp_rows[, l], l]))
                              checkIdent <- decomp$rank != 2 # rank == 2 means 2 vectors are linear independent
                          }
                          else {
                              checkIdent <- FALSE
                          }
                      }
                      else {
                          checkIdent <- identical(Stmp[Stmp_rows[, k], k], Stmp[Stmp_rows[, l], l])
                      }

                      #if (identical(Stmp[Stmp_rows[, k], k], Stmp[Stmp_rows[, l], l])) {
                      if (isTRUE(checkIdent)) {

                          # check the reversibilities
                          if (checkRev == TRUE) {
                              if (identical(react_rev(model)[columns[[i]][k]], react_rev(model)[columns[[i]][l]])) {
                                  checked <- TRUE
                              }
                          }
                          else {
                              checked <- TRUE
                          }

                          if (checked == TRUE) {
                              double_tmp  <- c(double_tmp, columns[[i]][k], columns[[i]][l])
                              ident_react <- c(ident_react, columns[[i]][k], columns[[i]][l])
                          }
                      }
                  }

                  l <- l + 1
              }

              # if ident_react contains something, we put it into the
              # list double_react (only different entries - unique)
              if (length(ident_react != 0)) {
                  clr                 <- clr + 1
                  double_react[[clr]] <- unique(ident_react)
              }
          }

          # if double_react contains something, we put it into the
          # main list double
          if (length(double_react) != 0) {
              clm           <- clm + 1
              #double_react[[1]] <- reactlength[i]
              #double[[clm]] <- double_react
              double <- append(double, double_react)
          }
      }
  }

  # Return FALSE if there are no double recations,
  # otherwise return a list with double reactions.
  if (length(double) != 0) {
      return(double)
  }
  else {
      return(FALSE)
  }

}
