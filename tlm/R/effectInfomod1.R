effectInfomod1 <-
function(object)
 {
  aux <- summary(object$model)$coefficients
  beta <- t(as.matrix(aux[2, ]))
  rownames(beta) <- rownames(aux)[2]
  Xincrease <- "additive of c units"
  effecttype <- "additive change in the mean of Y"
  effectsize <- "c * beta units of Y"
  furtherinfo <- "\nFurther details can be obtained using effect(), providing the increase in X\n, 'c',\nand the level for the confidence interval, 'level'."	
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res 	
 }
