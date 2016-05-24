effectInfomod12 <-
function(object)
 {
  aux <- summary(object$model)$coefficients
  beta <- t(as.matrix(aux[2, ]))
  rownames(beta) <- rownames(aux)[2]
  Xincrease <- "multiplicative of factor q (equivalently, adding an r = 100 * (q - 1)% to X)"
  effecttype <- "percent change in the mean of Y"
  effectsize <- "100 * (q^beta - 1)%"
  res <- list(beta = beta, Xincrease = Xincrease, effecttype = effecttype, effectsize = effectsize)
  res
 }
