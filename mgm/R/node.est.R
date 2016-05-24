
# fits a regularized glm using CV/EBIC for lambda selection
# as an extra function, because used twice, depending on method='linear' or 'glm'

node.est <- function(lambda.sel, y, fam, folds, emp_lev, v, n, nadj, gam, X, method, weights, type) {
  
  # lambda selection with EBIC
  
  if(lambda.sel == "EBIC") {
    
    fit <- glmnet(X, y, family = fam, alpha = 1, weights=weights)
    EBIC_obj <- calc.EBIC(fit, y, type, emp_lev, v, n, nadj, gam, X, method, weights) # calculates EBIC for a set of lambdas
    lambda_select <- EBIC_obj$lambda
    EBIC <- EBIC_obj$EBIC
    coefs <- coef(fit, s=lambda_select) #lambda value with lowest EBIC
    
    
    # lambda selection with CV
    
  } else {
    
    fit <- cv.glmnet(X, y, family=fam, alpha=1, nfolds=folds, type.measure = "deviance") # , weights=weights)
    lambda_select <-  fit$lambda.min
    EBIC <- NA
    coefs <- coef(fit, s=lambda_select)
  } 
  
  mod_list <- list('coefs' = coefs, 'EBIC' = EBIC, 'lambda_select' = lambda_select)
  
  return(mod_list)
}