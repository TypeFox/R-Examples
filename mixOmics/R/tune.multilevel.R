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
tune.multilevel <- function (X, Y = NULL, 
                             design, 
                             ncomp = 1, test.keepX = c(5, 10, 15), test.keepY = NULL, 
                             already.tested.X = NULL, already.tested.Y = NULL, 
                             method, mode, dist, validation = c("Mfold", "loo"), folds = 10) {
  
  X = as.matrix(X)
  if (length(dim(X)) != 2 || !is.numeric(X)) 
    stop("'X' must be a numeric matrix.")

  if(missing(design)) 
    stop('You forgot to input the design matrix.', call. = FALSE)
  
  if(is.null(design)) 
    stop("the 'design' matrix is missing.", call. = FALSE)

  if ((nrow(design) != nrow(X))) 
    stop("unequal number of rows in 'X' and 'design'.", call. = FALSE)
  
  # added condition for the spls case (no need to have the 2n and 3rd column in design)
  if ((ncol(design) < 2) & (method == 'splsda'))
    stop("'design' must be a matrix or data frame with at least 2 columns for method = splsda.",
         call. = FALSE)
  
  design = as.data.frame(design)

  if (length(design[, 1]) != nrow(X)) 
    stop("X and the vector sample should have the same number of subjects")
  
  if (is.factor(design[, 1])) {
    design[, 1] = as.numeric(design[, 1])
    warning("the vector sample was converted into a numeric vector", call. = FALSE)
  }

  if (length(summary(as.factor(design[, 1]))) == nrow(X)) 
    stop("Check that the vector sample reflects the repeated measurements")
  
  if (!any(names(summary(as.factor(design[, 1]))) == "1")) {
    cat("The vector sample includes the values: ", as.vector(names(summary(as.factor(sample)))), "\n")
    stop("sample vector", call. = FALSE)
  }

  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) 
    stop("invalid number of variates, 'ncomp'.")
  
  if (is.null(method)) 
    stop("Input method missing, should be set to splsda or spls", call. = FALSE)
  
  if ((method == "splsda") && (!is.null(already.tested.Y))) 
    warning("Unecessary already.tested.Y parameter")
  
  if (method == "splsda"){
    if (dim(design)[2] == 2) {
      result = tune.splsdalevel1(X = X, 
                                 design = design, 
                                 ncomp = ncomp, test.keepX = test.keepX, dist = dist, 
                                 already.tested.X = already.tested.X, validation, folds)
    } else {
      result = tune.splsdalevel2(X = X,  
                                 design = design, 
                                 ncomp = ncomp, test.keepX = test.keepX, already.tested.X = already.tested.X)
    }
  } else {
    result = tune.splslevel(X = X, Y = Y, 
                            design = design, 
                            mode = mode,
                            ncomp = ncomp, test.keepX = test.keepX, test.keepY = test.keepY, 
                            already.tested.X = already.tested.X, already.tested.Y = already.tested.Y)
  }
  return(result)
}
