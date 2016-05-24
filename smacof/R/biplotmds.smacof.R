biplotmds.smacof <- function(object, extvar, scale = TRUE) {
  
  if (any(class(object) == "smacofID")) X <- object$gspace else X <- object$conf
  p <- ncol(X)
  
  if (is.data.frame(extvar)) extvar <- as.matrix(extvar)
  ext <- scale(extvar, scale = TRUE)
  
  regfit <- lm(ext ~ -1 + X)
  rownames(regfit$coefficients) <- colnames(X)
  regsum <- summary(regfit)
  R2vec <- sapply(regsum, `[[`, "r.squared")
  names(R2vec) <- gsub("Response ", "", names(R2vec))
  regfit$R2vec <- R2vec
  class(regfit) <- c("mdsbi", "mlm", "lm")
  return(regfit)
}

biplotmds.smacofID <- biplotmds.smacof 