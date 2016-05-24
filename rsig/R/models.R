fit_lasso = function(surv, X, cv.folds = 10L, sparse.output=TRUE) {
  prof = cv.glmnet(X, surv, family = "cox", nfolds = cv.folds, alpha = 1)
  lambda = prof$lambda.min
  beta = coef(prof, s = "lambda.min") # this is a sparse object, so return it as it is.
  
  intercept = 0
  if(rownames(beta)[1] == "(Intercept)") {
    intercept = beta[1]
    beta = beta[2:nrow(beta),,drop=FALSE]
  }
  
  if(!sparse.output) {
    beta = as.matrix(beta)
    beta = beta[which(abs(beta) > sqrt(.Machine$double.eps)),,drop=FALSE]
  }

  return(list(beta=beta, intercept=intercept))
}

fit_prlasso = function(surv, X, cv.folds = 10L, sparse.output=FALSE, ncomp = 1L) {

  data = list(x = t(X), y = surv[, 1L], censoring.status = surv[, 2L],
              featuenames = colnames(X))

  ### Part 1: preconditioning
#   options(warn=-1) # to suppress warning from svd(x, LINPACK=TRUE) in superpc.

  
  suppressAll({
    fit = superpc.train(data, type = "survival")
    cv = superpc.cv(fit, data, n.components = ncomp)
  })
  
  maxidx = which.max(cv$scor)
  maxcol = floor((maxidx-1)/ncomp) + 1L  # which threshold?
  maxrow = ifelse(maxidx %% ncomp == 0, ncomp, maxidx%%ncomp) # up to which order of pc's?
  thres = cv$threshold[ maxcol ]
  ncomp.cv = maxrow
#   thres.new = cv$threshold[ max.col(cv$scor) ][ncomp]
#   cat('thres = ', thres, ', thres.new = ', thres.new, '\n')
  pred.tr = superpc.predict(fit, data, data, threshold = thres, n.components=ncomp.cv, prediction.type="continuous")$v.pred.1df

#   options(warn=0)
  
  ### Part 2: linear regression with Lasso penality
  prof = cv.glmnet(t(data$x), pred.tr, family="gaussian", nfolds=cv.folds, alpha=1, standardize=TRUE)
 
  lambda = prof$lambda.min
  beta = coef(prof, s = "lambda.min") # this is a sparse object, so return it as it is.
  
  intercept = 0
  if(rownames(beta)[1] == "(Intercept)") {
    intercept = beta[1]
    beta = beta[2:nrow(beta),,drop=FALSE]
  }
  
  if(!sparse.output) {
    beta = as.matrix(beta)
    beta = beta[which(abs(beta) > sqrt(.Machine$double.eps)),,drop=FALSE]
  }

  return(list(beta=beta, intercept=intercept))
}
