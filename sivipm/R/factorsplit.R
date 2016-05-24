
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

factorsplit <- function(data) {
  # Split each factor column of a dataframe into as many
  # 0-1 columns as occuring levels in the factor when it has more than 2 levels
  # or transform the column in -1/+1 column  when it has 2 levels. 
    label <- ""
  datares <- matrix(0, nrow=nrow(data), ncol=1)

   for (icol in 1:ncol(data)) {
     if (is.factor(data[, icol])) {
        nmoda <- nlevels(data[, icol])
       if (nmoda ==1) {
         stop(paste("factorsplit. Invalid variable ", icol,
                    ": Factors should have several levels"))
       }
        if (nmoda > 2) {
          # more than 2 levels => nmoda variables
       z <- matrix(0, nrow=nrow(data), ncol=nmoda)
       for (il in 1:nmoda) {
         lignes <- (as.integer(data[, icol])==il)
         z[lignes, il] <- 1
         label <- c(label,
           paste(colnames(data)[icol], levels(data[, icol])[il], sep="_"))
       } # fin il
     } else {
       # 2 levels => a  -1/+1 variable
       leslevels <- levels(data[, icol])
       z <- rep(1, nrow(data))
       z[data[, icol]==leslevels[1]] <- -1
       label <- c(label, colnames(data)[icol])
      }
        
       datares <- cbind(datares, z)
     } else {
         datares <- cbind(datares, data[, icol])
         label <- c(label,colnames(data)[icol])
       }

} # fin icol

  col <- apply(datares, 2, function(X) !all(X==0))
  datares <- datares[,col]
    dimnames(datares)[[2]] <- label[col]
    nomvar <- colnames(datares)
    cat("The names of the variables are:\n")
    names(nomvar) <- 1:length(nomvar)
    print(nomvar)
    cat( "\n")
    datares <- as.data.frame(datares)
  return(datares)
}
