effectInfomod10 <-
function(object)
 {
  aux <- summary(object$model)$coefficients
  beta <- t(as.matrix(aux[2, ]))
  rownames(beta) <- rownames(aux)[2]
  Xincrease <- "additive of c units"
  effecttype <- "percent change in the mean of Y"
  effectsize <- "100 * [exp(c * beta) - 1]%"
  furtherinfo <- "\nFurther details can be obtained using effect(), providing the increase in X, 'c', and the\nlevel for the confidence interval, 'level'."  
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res
 }
