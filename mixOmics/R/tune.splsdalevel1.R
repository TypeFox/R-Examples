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
tune.splsdalevel1 <- function (X, ncomp = 1, test.keepX, dist = NULL, 
                               design = NULL, 
                              already.tested.X = NULL, validation, folds) {
  
  if (validation == "Mfold") {
    cat(paste("For a one-factor analysis, the tuning criterion is based on", folds, "cross validation"), "\n")
  } else {
    cat(paste("For a one-factor analysis, the tuning criterion is based on leave-one-out cross validation"), "\n")
  }

  if (!is.factor(design[, 2])) {
    design[, 2] = as.factor(design[, 2])
    warning("design[, 2] was set as a factor", call. = FALSE)
  }
  
  if (is.null(dist)) 
    stop("Input dist is missing")
  
  if (length(already.tested.X) != (ncomp - 1)) 
    stop("The number of already tested parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
  
  if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X))) 
    stop("Expecting a numerical value in already.tested.X", call. = FALSE)
  
  if (!is.null(already.tested.X)) 
    cat("Number of variables selected on the first ", ncomp - 1, "component(s) was ", already.tested.X)
  
  sample = design[, 1]# *BG* add sample feature
  n = length(unique(sample))
  
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n) 
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n) 
        stop("Invalid folds.")
      M = length(folds)
    } else {
      if (is.null(folds) || !is.numeric(folds) || folds < 
            2 || folds > n) 
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n))
      }
    }
  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
  }
  
  vect.error = vector(length = length(test.keepX))
  names(vect.error) = paste("var", test.keepX, sep = "")
  error.sw = matrix(nrow = M, ncol = length(test.keepX), data = 0)
  # for the last keepX (i) tested, prediction combined for all M folds so as to extract the error rate per class
  prediction.all = vector(length = length(sample))
  
  # calculate first the within variation matrix, which is ok because Xw is not dependent on each subject (so there is not bias here)
  Xw = suppressMessages(withinVariation(X = X, design = design)) # Estimation of the design matrix. Outside the loop because Xw is "subject specific"

  for (j in 1:M) {
    omit = which(sample %in% folds[[j]] == TRUE)
    trainxw = Xw[-c(omit), ]
    xwtest = Xw[omit, , drop = FALSE]

    
    cond.train = design[-c(omit), 2] 
    cond.test = design[c(omit), 2] 

    remove.zero = nearZeroVar(trainxw)$Position
    if (length(remove.zero) != 0) {
      trainxw = trainxw[, -c(remove.zero)]
      xwtest = xwtest[, -c(remove.zero)]
    }
    
    for (i in 1:length(test.keepX)) {
      if (ncomp == 1) {
        result.sw <- splsda(trainxw, cond.train, keepX = test.keepX[i], ncomp = ncomp)
      }
      else {
        result.sw <- splsda(trainxw, cond.train, keepX = c(already.tested.X, test.keepX[i]), ncomp = ncomp)
      }
      test.predict.sw <- predict(result.sw, xwtest, dist = dist)
      Prediction.sw <- levels(design[, 2])[test.predict.sw$class[[dist]][, ncomp]] 
      error.sw[j, i] <- sum(as.character(cond.test) != Prediction.sw)
    } # end i
    
    # for the last keepX (i) tested, prediction combined for all M folds:
    prediction.all[omit] = Prediction.sw
    
  } # end J (M folds)
  result <- apply(error.sw, 2, sum)/length(design[, 2]) 

  names(result) = paste("var", test.keepX, sep = "")
return(list(
  error = result,
  # error per class for last keepX tested
  prediction.all = prediction.all))
}