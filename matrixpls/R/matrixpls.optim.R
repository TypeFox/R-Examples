# =========== Optimization criterion functions ===========

#'@title Optimization criteria functions
#'
#'@description Optimization criterion functions calculate various optimization criterion values
#'from \code{matrixpls} objects. 
#'
#'@param matrixpls.res An object of class \code{matrixpls} from which the
#'criterion function is calculated
#'
#'@return Value of the optimization criterion.
#'
#'@seealso \code{\link{weight.optim}}
#'@name optimCriterion
NULL


#'@describeIn optimCriterion maximizes the sum of R2 statistics of the \code{inner} matrix
#'@export

optim.maximizeInnerR2 <- function(matrixpls.res){
  -sum(r2(matrixpls.res))
}


#'@describeIn optimCriterion maximizes the sum of R2 statistics of the \code{reflective} matrix.
#'@export

optim.maximizeIndicatorR2 <- function(matrixpls.res){
  lambda <- loadings(matrixpls.res)
  IC <- attr(matrixpls.res,"IC")
  -sum(diag(lambda %*% IC))
}

#'@describeIn optimCriterion maximizes the sum of R2 statistics of the \code{inner} and \code{reflective} matrices.
#'@export


optim.maximizeFullR2 <- function(matrixpls.res){
  optim.maximizeIndicatorR2(matrixpls.res) + optim.maximizeInnerR2(matrixpls.res)
}

#'@describeIn optimCriterion minimies the generalized structured component analysis criterion. See \link{GSCA}
#'@export

optim.gsca <- function(matrixpls.res){
  
  C <- attr(matrixpls.res,"C")
  IC <- attr(matrixpls.res,"IC")
  nativeModel <- attr(matrixpls.res,"model")
  
  reflective <- nativeModel$reflective
  reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
  
  #  formative <- nativeModel$formative
  #  formative[which(formative==1)] <- matrixpls.res[grep("<~",names(matrixpls.res))]
  
  #  f <- apply(nativeModel$formative != 0,1,any)
  r <- apply(nativeModel$reflective != 0,1,any)
  endo <- apply(nativeModel$inner != 0,1,any)
  
  inner_resid <- (1 - r2(matrixpls.res)[endo])
  #  form_resid <- (1 - rowSums(IC[f,] * formative[f,]))
  refl_resid <- (1 - rowSums(t(IC[,r]) * reflective[r,]))
  
  #  sum(inner_resid, form_resid, refl_resid)
  sum(inner_resid, refl_resid)
}


# #'@title Optimization criterion for maximal prediction 
# #'
# #'@details Calculates the predicted variances of reflective indicators. The
# #'prediction criterion is negative of the sum of predicted variances.
# #'
# #'@param matrixpls.res An object of class \code{matrixpls} from which the
# #'criterion function is calculated
# #'
# #'@return Mean squared prediction error.
# #'
# #'@family Weight optimization criteria
# #'
# #'@export
# #'
# 
# optim.maximizePrediction <- function(matrixpls.res){
#   
#   # TODO: convert to using the predict function
#   
#   C <- attr(matrixpls.res,"C")
#   nativeModel <- attr(matrixpls.res,"model")
#   exog <- rowSums(nativeModel$inner)==0
#   W <- attr(matrixpls.res,"W")
#   inner <- attr(matrixpls.res,"inner")
#   
#   reflective <- nativeModel$reflective
#   reflective[which(reflective==1)] <- matrixpls.res[grep("=~",names(matrixpls.res))]
#   
#   # Predict endog LVs using reduced from equations
#   Beta <- inner[! exog, ! exog]
#   Gamma <- inner[! exog, exog]
#   
#   invBeta <- solve(diag(nrow(Beta)) - Beta)
#   endogC <- invBeta %*% Gamma %*% C[exog,exog] %*% t(Gamma) %*% t(invBeta)
#   C[!exog, !exog] <- endogC
#   
#   -sum(diag(reflective %*% C %*% t(reflective)))
# }
