setMethod("worth", "btmix", function(object, component = NULL){
  cf0 <- parameters(object, component = component)
  rownames(cf0) <- gsub("coef.", "", rownames(cf0), fixed = TRUE)
  cf <- matrix(0, nrow = nrow(cf0) + 1, ncol = ncol(cf0))
  colnames(cf) <- colnames(cf0)
  rownames(cf) <- labels(object)
  cf[rownames(cf0), ] <- cf0
  apply(cf, 2, function(x) exp(x)/sum(exp(x)))
})

labels.btmix <- function(object, ...) object@labels
mscale.btmix <- function(object, ...) object@mscale
