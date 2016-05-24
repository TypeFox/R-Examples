effectInfomod5 <-
function(object)
 {
  aux <- summary(object$model)$coefficients
  beta <- t(as.matrix(aux[2, ]))
  rownames(beta) <- rownames(aux)[2] 	
  Xincrease <- "multiplicative of factor q (equivalently, adding an r = 100 * (q - 1)% to X)"
  effecttype <- "additive change in the mean of Y"
  effectsize <- "beta * log(q) units of Y"
  furtherinfo <- "\nFurther details can be obtained using effect(), providing either the multiplicative ('q') or\nthe percent ('r') change in X, and the level for the confidence interval, 'level'."
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res 	
 }
