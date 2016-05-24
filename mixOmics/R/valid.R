# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh La Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#
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


valid <-function(object, ...){    #UseMethod("valid")
 cat('The valid function has been supersided by the perf function to evaluate the performance of the PLS, sPLS, PLS-DA and sPLS-DA methods. 
     See ?perf for more details', '\n')
}

# #------------------------------------------------------#
# #-- announces deprectation of function -------------- --#
# #------------------------------------------------------#
# 
# # ---------------------------------------------------
# # valid for pls object
# # ---------------------------------------------------
# valid.pls <-
#   function(object, 
#            validation = c("Mfold", "loo"),
#            folds = 10,
#            max.iter = 500, 
#            tol = 1e-06, ...)
#   {
#     
#     #-- validation des arguments --#
#     X = object$X
#     Y = object$Y
#     
#     if (length(dim(X)) != 2) 
#       stop("'X' must be a numeric matrix for validation.")
#     
#     # this will be used only for loocv computation (see below)
#     means.Y = attr(scale(Y), "scaled:center")
#     sigma.Y = attr(scale(Y), "scaled:scale")
#     
#     validation = match.arg(validation)
#     
#     mode = object$mode
#     ncomp = object$ncomp
#     n = nrow(X)
#     p = ncol(X)
#     q = ncol(Y)
#     res = list()
#     
#     if (ncomp == 1)
#       warning("'ncomp' must be > 2 to compute the Q2 criterion.")
#     
#     
#     if (any(is.na(X)) || any(is.na(Y))) 
#       stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
#     
#     
#     #-- define  M fold or loo cross validation --------------------#
#     #- define the folds
#     if (validation == "Mfold") {
#       if (is.list(folds)) {
#         if (length(folds) < 2 | length(folds) > n)
#           stop("Invalid number of folds.")
#         if (length(unique(unlist(folds))) != n)
#           stop("Invalid folds.")
#         
#         M = length(folds)
#       }
#       else {
#         if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
#           stop("Invalid number of folds.")
#         else {
#           M = round(folds)
#           folds = split(sample(1:n), rep(1:M, length = n)) 
#         }
#       }
#     } else { 
#       folds = split(1:n, rep(1:n, length = n)) 
#       M = n
#     }
#     
#     #-- compute the criteria --#
#     RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
#     RSS.indiv = array(0, c(n, q, ncomp+1))
#     PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
#     press.mat = Ypred = array(0, c(n, q, ncomp))
#     MSEP = R2 = matrix(0, nrow = q, ncol = ncomp)
#     
#     # in case the test set only includes one sample, it is better to advise the user to perform loocv
#     stop.user = FALSE
#     for (i in 1:M) {
#       omit = folds[[i]]
#       
#       # see below, we stop the user if there is only one sample drawn on the test set using MFold
#       if(length(omit) == 1) stop.user = TRUE
#       
#       # the training set is scaled
#       X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
#       Y.train = scale(Y[-omit, ], center = TRUE, scale = TRUE)
#       
#       # the test set is scaled either in the predict function directly (for X.test)
#       # or below for Y.test
#       X.test = matrix(X[omit, ], nrow = length(omit))
#       Y.test = matrix(Y[omit, ], nrow = length(omit))
#       
#       
#       #-- pls --#
#       result = pls(X = X.train, Y = Y.train, ncomp = ncomp, 
#                    mode = mode, max.iter = max.iter, tol = tol)
#       
#       if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
#       # in the predict function, X.test is already normalised w.r.t to training set X.train, so no need to do it here
#       Y.hat = predict(result, X.test)$predict
#       
#       for (h in 1:ncomp) {
#         Ypred[omit, , h] = Y.hat[, , h]
#         
#         # compute the press and the RSS
#         # par definition de tenenhaus, RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
#         if(validation == 'Mfold'){
#           # for Mfold, Y.test is simply scaled (seemed to be ok when there are enough samples per fold)
#           press.mat[omit, , h] = (scale(Y.test) - Y.hat[, , h])^2
#           RSS.indiv[omit, ,h+1] = (scale(Y.test) - Y.hat[, , h])^2
#         } else{ 
#           # in the case of loo we need to scale w.r.t the parameters in Y.train
#           Y.test = sweep(Y.test, 2, means.Y, FUN = "+")
#           Y.test = sweep(Y.test, 2, sigma.Y, FUN = "*")
#           press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
#           RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
#         }
#       } # end h
#     } #end i (cross validation)
#     # warn the user that at least test set had a length of 1
#     if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
#     
#     
#     for (h in 1:ncomp) { 
#       
#       
#       
#       # compute pRESS and Q2
#       if(q > 1){
#         MSEP[, h] = apply(press.mat[, , h], 2, mean, na.rm = TRUE)
#         R2[, h] = (diag(cor(scale(Y), Ypred[, , h], use = "pairwise")))^2
#         RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum))
#         PRESS.inside[h, ] = colSums(press.mat[, , h], na.rm = TRUE)
#       } else {  # if q == 1
#         MSEP[q, h] = mean(press.mat[, q, h], na.rm = TRUE)
#         R2[1, h] = (diag(cor(scale(Y), Ypred[, , h], use = "pairwise")))^2
#         
#         RSS[h+1,] = sum(RSS.indiv[,q,h+1], na.rm = TRUE)
#         PRESS.inside[h, ] = colSums(as.matrix(press.mat[, , h]), na.rm = TRUE)
#         ##Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
#       } # end if q
#       
#       Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
#     } # end h
#     
#     
#     colnames(MSEP) = colnames(R2) = paste('ncomp', c(1:ncomp), sep = " ")
#     rownames(MSEP) = rownames(R2) = colnames(Y)
#     
#     if (q == 1){
#       rownames(MSEP) = rownames(R2) = ""      
#       Q2.total = 1 - rowSums(as.matrix(PRESS.inside), na.rm = TRUE)/rowSums(as.matrix(RSS[-(ncomp+1), ]), na.rm = TRUE)
#     }
#     
#     if (q > 1){
#       Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
#     }
#     
#     Y.names = dimnames(Y)[[2]]
#     if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
#     
#     if (q > 1) {
#       colnames(Q2.inside) = Y.names
#       rownames(Q2.inside) = paste('comp', 1:ncomp, sep = " ")
#       names(Q2.total) = paste('comp', 1:ncomp, sep = " ")    
#     }
#     else {
#       #colnames(Q2.inside) = ""
#       names(Q2.inside) = names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
#     }
#     
#     
#     method = "pls.mthd"
#     class(res) = c("valid", method)
#     return(      
#       list(MSEP  = MSEP,
#            R2=  R2,
#            Q2 = t(Q2.inside),
#            Q2.total = Q2.total))
#   }
# 
# 
# # ===========================================================================================
# # ---------------------------------------------------
# # valid for spls object
# # ---------------------------------------------------
# valid.spls <-
#   function(object, 
#            validation = c("Mfold", "loo"),
#            folds = 10,
#            max.iter = 500, 
#            tol = 1e-06, ...)
#   {
#     
#     #-- validation des arguments --#
#     X = object$X
#     Y = object$Y
#     # tells which variables are selected in X and in Y:
#     keepX = (object$loadings$X != 0) 
#     keepY = (object$loadings$Y != 0)
#     
#     if (length(dim(X)) != 2) 
#       stop("'X' must be a numeric matrix for validation.")
#     
#     # these parameters are only used for the specific case of normalising the test sets for LOO CV
#     means.Y = attr(scale(Y), "scaled:center")
#     sigma.Y = attr(scale(Y), "scaled:scale")
#     
#     validation = match.arg(validation)
#     
#     mode = object$mode
#     ncomp = object$ncomp
#     n = nrow(X)
#     p = ncol(X)
#     q = ncol(Y)
#     res = list()
#     
#     if (ncomp == 1)
#       warning("'ncomp' must be > 2 to compute the Q2 criterion.")
#     
#     if (any(is.na(X)) || any(is.na(Y))) 
#       stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
#     
#     #-- define the folds M fold or loo cross validation -------#
#     if (validation == "Mfold") {
#       if (is.list(folds)) {
#         if (length(folds) < 2 | length(folds) > n)
#           stop("Invalid number of folds.")
#         if (length(unique(unlist(folds))) != n)
#           stop("Invalid folds.")
#         
#         M = length(folds)
#       }
#       else {
#         if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
#           stop("Invalid number of folds.")
#         else {
#           M = round(folds)
#           folds = split(sample(1:n), rep(1:M, length = n)) 
#         }
#       }
#     } 
#     else { 
#       folds = split(1:n, rep(1:n, length = n)) 
#       M = n
#     }
#     
#     #-- compute the different criteria --#
#     RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
#     RSS.indiv = array(NA, c(n, q, ncomp+1))
#     PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
#     press.mat = Ypred = array(NA, c(n, q, ncomp))
#     MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
#     
#     # in case the test set only includes one sample, it is better to advise the user to
#     # perform loocv
#     stop.user = FALSE
#     
#     for (i in 1:M) {
#       omit = folds[[i]]
#       
#       # see below, we stop the user if there is only one sample drawn on the test set using MFold
#       if(length(omit) == 1) stop.user = TRUE
#       
#       # the training set is scaled
#       X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
#       Y.train = scale(Y[-omit, ], center = TRUE, scale = TRUE)
#       # the test set is scaled either in the predict function directly (for X.test)
#       # or below for Y.test
#       X.test = matrix(X[omit, ], nrow = length(omit))
#       Y.test = matrix(Y[omit, ], nrow = length(omit))
#       
#       #-- spls --#
#       result = spls.model(X.train, Y.train, ncomp, mode, 
#                           max.iter, tol, keepX, keepY)
#       
#       if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
#       Y.hat = predict(result, X.test)$predict
#       
#       for (h in 1:ncomp) {
#         Ypred[omit, , h] = Y.hat[, , h]
#         
#         # compute the press and the RSS
#         # Tenenhaus definition: RSS[h+1,] = (y_i - y.hat_(h-1)_i)^2
#         if(validation == 'Mfold'){
#           # for Mfold, Y.test is simply scaled (seemed to be ok when there are enough samples per fold)
#           press.mat[omit, , h] = (scale(Y.test) - Y.hat[, , h])^2
#           RSS.indiv[omit, ,h+1] = (scale(Y.test) - Y.hat[, , h])^2
#         }else{ 
#           # in the case of loo we need to scale w.r.t the parameters in Y.train
#           Y.test = sweep(Y.test, 2, means.Y, FUN = "-")
#           Y.test = sweep(Y.test, 2, sigma.Y, FUN = "/")
#           press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
#           RSS.indiv[omit, ,h+1] = (Y.test - Y.hat[, , h])^2
#         }
#       } # end h
#     } #end i (cross validation)
#     
#     # warn the user that at least test set had a length of 1
#     if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
#     
#     for (h in 1:ncomp) { 
#       if(q >1){
#         MSEP[keepY[, h], h] = apply(as.matrix(press.mat[, keepY[, h], h]), 2, mean, na.rm = TRUE)
#         RSS[h+1,] = t(apply(RSS.indiv[,,h+1], 2, sum, na.rm = TRUE))
#       }else{
#         MSEP[keepY[, h], h] = mean(press.mat[, keepY[, h], h], na.rm = TRUE)
#         RSS[h+1,] = sum(RSS.indiv[,,h+1], na.rm = TRUE)
#       }
#       
#       if(sum(keepY[,h]==TRUE) >1){
#         R2[keepY[, h], h] = (diag(cor(Y[, keepY[, h]], Ypred[, keepY[, h], h], use = "pairwise")))^2
#       } else{
#         R2[keepY[, h], h] = cor(Y[, keepY[, h]], Ypred[, keepY[, h], h], use = "pairwise")^2
#       }
#       PRESS.inside[h, ] = colSums(as.matrix(press.mat[, , h]), na.rm = TRUE)
#       Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
#       
#     }
#     
#     colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste('ncomp', c(1:ncomp), sep = " ")
#     rownames(MSEP) = rownames(R2) = colnames(Q2.inside)  = colnames(Y)
#     
#     if (q == 1){
#       rownames(MSEP) = rownames(R2) = ""
#       Q2.total = 1 - rowSums(as.matrix(PRESS.inside), na.rm = TRUE)/rowSums(as.matrix(RSS[-(ncomp+1), ]), na.rm = TRUE)
#     }
#     
#     if (q > 1) {
#       Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp+1), ], na.rm = TRUE)
#     }
#     
#     Y.names = dimnames(Y)[[2]]      
#     if (is.null(Y.names)) Y.names = paste("Y", 1:q, sep = "")
#     
#     if (q > 1) {
#       colnames(Q2.inside) = Y.names
#       rownames(Q2.inside) = paste('comp', 1:ncomp, sep = " ")
#       names(Q2.total) = paste('comp', 1:ncomp, sep = " ")    
#     }
#     else {
#       #colnames(Q2.inside) = ""
#       names(Q2.inside) = names(Q2.total) = paste('comp', 1:ncomp, sep = " ")
#     }
#     
#     method = "pls.mthd"
#     class(res) = c("valid", method)
#     return(      
#       list(MSEP  = MSEP,
#            R2=  R2,
#            Q2 = t(Q2.inside),
#            Q2.total = Q2.total))
#   }
# 
# 
# # ---------------------------------------------------
# # valid for plsda object
# # ---------------------------------------------------
# valid.plsda <-
#   function(object,
#            method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
#            validation = c("Mfold", "loo"),
#            folds = 10,
#            max.iter = 500, 
#            tol = 1e-06, ...)
#   {
#     
#     #-- validation des arguments --#
#     X = object$X
#     lev = object$names$Y
#     Y = object$ind.mat
#     Y = map(Y)
#     Y = as.factor(lev[Y])
#     ncomp = object$ncomp
#     n = nrow(X)
#     
#     method = match.arg(method, several.ok = TRUE)
#     if (any(method == "all")) nmthd = 3 
#     else nmthd = length(method)
#     
#     error.fun = function(x, y) {
#       error.vec = sweep(x, 1, y, FUN = "-")
#       error.vec = (error.vec != 0)
#       error.vec = apply(error.vec, 2, sum) / length(y)
#       return(error.vec)
#     }
#     
#     #-- define the folds --#
#     if (validation == "Mfold") {
#       if (is.list(folds)) {
#         if (length(folds) < 2 | length(folds) > n)
#           stop("Invalid number of folds.")
#         if (length(unique(unlist(folds))) != n)
#           stop("Invalid folds.")
#         
#         M = length(folds)
#       }
#       else {
#         if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
#           stop("Invalid number of folds.")
#         else {
#           M = round(folds)
#           folds = split(sample(1:n), rep(1:M, length = n)) 
#         }
#       }
#     } 
#     else { 
#       folds = split(1:n, rep(1:n, length = n)) 
#       M = n
#     }
#     
#     error.mat = array(0, dim = c(ncomp, nmthd, M))
#     
#     # in case the test set only includes one sample, it is better to advise the user to
#     # perform loocv
#     stop.user = FALSE
#     
#     
#     for (i in 1:M) {
#       omit = folds[[i]]
#       
#       # see below, we stop the user if there is only one sample drawn on the test set using MFold
#       if(length(omit) == 1) stop.user = TRUE
#       
#       # the training set is scaled
#       X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
#       Y.train = Y[-omit]
#       
#       X.test = matrix(X[omit, ], nrow = length(omit))
#       
#       # run pls-da
#       result = plsda(X = X.train, Y = Y.train, ncomp = ncomp, 
#                      max.iter = max.iter, tol = tol)
#       
#       if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
#       Y.hat = predict(result, X.test, method = method)$class
#       error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit]))
#     }
#     
#     # warn the user that at least test set had a length of 1
#     if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
#     
#     
#     #-- compute the error --#
#     res = apply(error.mat, 1:2, mean)
#     rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
#     colnames(res) = names(Y.hat)
#     
#     method = "plsda.mthd"
#     class(res) = c("valid", method)
#     return(invisible(res))
#   }
# 
# # ---------------------------------------------------
# # valid for splsda object
# # ---------------------------------------------------
# valid.splsda <-
#   function(object,
#            method = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
#            validation = c("Mfold", "loo"),
#            folds = 10,
#            max.iter = 500, 
#            tol = 1e-06, ...)
#   {
#     
#     #-- validation des arguments --#
#     X = object$X
#     lev = object$names$Y
#     Y = object$ind.mat
#     Y = map(Y)
#     Y = as.factor(lev[Y])
#     ncomp = object$ncomp
#     n = nrow(X)
#     keepX = (object$loadings$X != 0)
#     keepY = (object$loadings$Y != 0)
#     
#     method = match.arg(method, several.ok = TRUE)
#     if (any(method == "all")) nmthd = 3 
#     else nmthd = length(method)  
#     
#     error.fun = function(x, y) {
#       error.vec = sweep(x, 1, y, FUN = "-")
#       error.vec = (error.vec != 0)
#       error.vec = apply(error.vec, 2, sum) / length(y)
#       return(error.vec)
#     }
#     
#     #-- define the folds --#
#     if (validation == "Mfold") {
#       if (is.list(folds)) {
#         if (length(folds) < 2 | length(folds) > n)
#           stop("Invalid number of folds.")
#         if (length(unique(unlist(folds))) != n)
#           stop("Invalid folds.")
#         
#         M = length(folds)
#       }
#       else {
#         if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
#           stop("Invalid number of folds.")
#         else {
#           M = round(folds)
#           folds = split(sample(1:n), rep(1:M, length = n)) 
#         }
#       }
#     } 
#     else { 
#       folds = split(1:n, rep(1:n, length = n)) 
#       M = n
#     }
#     
#     error.mat = array(0, dim = c(ncomp, nmthd, M))
#     
#     # in case the test set only includes one sample, it is better to advise the user to
#     # perform loocv
#     stop.user = FALSE
#     
#     for (i in 1:M) {
#       omit = folds[[i]]
#       
#       # see below, we stop the user if there is only one sample drawn on the test set using MFold
#       if(length(omit) == 1) stop.user = TRUE
#       X.train = scale(X[-omit, ], center = TRUE, scale = TRUE)
#       Y.train = Y[-omit]
#       X.test = matrix(X[omit, ], nrow = length(omit))
#       
#       # run sPLS-DA  
#       result = splsda.model(X.train, Y.train, ncomp, 
#                             max.iter, tol, keepX, keepY)
#       
#       if (!is.null(result$nzv$Position)) X.test = X.test[, -result$nzv$Position]
#       Y.hat = predict(result, X.test, method = method)$class
#       error.mat[, , i] = sapply(Y.hat, error.fun, y = as.numeric(Y[omit]))
#     }
#     
#     # warn the user that at least test set had a length of 1
#     if(stop.user == TRUE & validation == 'Mfold') stop('The folds value was set too high to perform cross validation. Choose validation = "loo" or set folds to a lower value')
#     
#     
#     #-- compute the error --#
#     res = apply(error.mat, 1:2, mean)
#     rownames(res) = paste('ncomp', 1:ncomp, sep = " ")
#     colnames(res) = names(Y.hat)
#     
#     method = "plsda.mthd"
#     class(res) = c("valid", method)
#     return(invisible(res))
#   }
# 
# 
