# Copyright (C) 2015 
# Kim-Anh Le Cao, University of Queensland, Brisbane, Australia
# Benoit Gautier, University of Queensland, Brisbane, Australia
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


## note: this should probably be set up as an S3 function.
tune.splsdalevel2 <- function (X, 
                               design,
                               ncomp = 1, test.keepX = c(5, 10, 15), already.tested.X = NULL) {
  
  cat("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set", "\n")
  
  if (!is.null(already.tested.X)) 
    cat("Number of variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X)
  
  if (length(already.tested.X) != (ncomp - 1)) 
    stop("The number of already tested parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  # *BG* Start: check condition design matrix instead of cond
#   if ((is.matrix(cond)) && ncol(cond) > 2) 
#     stop("'cond' must be a matrix with 2 columns for the 2 factor analysis.")
#   
#   if (ncol(cond) == 2) {
#     if (!is.factor(cond[, 1])) {
#       warning("First cond response was set as a factor", call. = FALSE)
#     }
#     if (!is.factor(cond[, 2])) {
#       warning("Second cond response was set as a factor", call. = FALSE)
#     }
#     cond1 = as.factor(cond[, 1])
#     cond2 = as.factor(cond[, 2])
#   }

  if (is.data.frame(design[, -1]) && ncol(design[, -1]) > 2) 
    stop("'cond' must be a matrix with 2 columns for the 2 factor analysis.")

  if (ncol(design[, -1]) == 2){
    if (!is.factor(design[, 2])) {
      design[, 2] = as.factor(design[, 2])
      warning("First factor response was set as a factor", call. = FALSE)
    }
    if (!is.factor(design[, 3])) {
      design[, 3] = as.factor(design[, 3])
      warning("Second factor response was set as a factor", call. = FALSE)
    }
  }    
  
  cond.fact = as.factor(paste(design[, 2], design[, 3], sep = ".")) 

  Xw <- suppressMessages(withinVariation(X = X, design = design)) 
  cor.value = vector(length = length(test.keepX))
  names(cor.value) = paste("var", test.keepX, sep = "")
  
  for (i in 1:length(test.keepX)) {
    if (ncomp == 1) {
      spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = test.keepX[i])
    } else {
      spls.train = splsda(Xw, cond.fact, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]))
    }
    
    # Note: this is performed on the full data set
    # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
    cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
    # 
  }
  return(list(cor.value = cor.value))
}