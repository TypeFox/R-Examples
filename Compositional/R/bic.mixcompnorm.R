################################
#### BIC for normal mixture models for compositional data
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
#### References: Ryan P. Browne, Aisha ElSherbiny and
#### Paul D. McNicholas (2015)
#### R package mixture: Mixture Models for Clustering and Classification
################################

bic.mixcompnorm <- function(x, A, type = "alr") {
  ## x is the compositional data
  ## A is the maximum number of components to be considered
  ## type is either 'alr' or 'ilr'
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure x is compositional
  p <- ncol(x)  ## dimensionality of the data

  if (type == "alr") {
    y <- log(x[, -p]/x[, p])
  } else {
    y0 <- log(x)
    y1 <- y0 - rowMeans( y0 )
    y <- tcrossprod( y1, helm(p) )

  }

  mod <- mixture::gpcm(y, G = 1:A, start = 0, atol = 0.01)
  bic <- mod$BIC[, , 3]  ## BIC for all models
  ## Next, we plot the BIC for all models
  plot(1:A, bic[, 1], type = "b", pch = 9, xlab = "Number of components",
  ylab = "BIC values", ylim = c(min(bic, na.rm = T), max(bic, na.rm = T)))
  for (i in 2:nrow(bic)) lines(1:A, bic[, i], type = "b", pch = 9, col = i)
  list(mod = mod, BIC = bic)
}
