perc.cis <- function(object, ncomp = object$ncomp, conf = .95, 
                     type = c("coefficients",
                              "loadings", "weights")) {
  if ((object$val.method == "none" | object$val.method == "loo")) {
    stop("No bootstrapping was done for this model")
  }
  conf <- conf
  Bootstraps <- do.call("rbind", object$validation[names(object$validation) == type][[1]])
  Bootstraps <- as.matrix(Bootstraps[, ncomp]) 
  Upper <- 1 - (((1 - conf)/2))
  Lower <- 1 - Upper
  Order <- as.factor(row.names(data.frame((object[names(object) == type])[[1]][, 1])))
  if(length(ncomp) > 1) {
    Bootstraps.boot.cis <- llply(ncomp, function(y) {
      do.call("rbind", as.list(
        by(Bootstraps[, y], list(row.names(Bootstraps)), function(x){
          c(ncomp = y, boot.mean = mean(x), quantile(x, c(Lower, Upper), na.rm = T))
        }
        )))
    })
  } else {
    Bootstraps.boot.cis <- llply(ncomp, function(y) {
      do.call("rbind", as.list(
        by(Bootstraps[, 1], list(row.names(Bootstraps)), function(x){
          c(ncomp = y, boot.mean = mean(x), quantile(x, c(Lower, Upper), na.rm = T))
        }
        )))
    })
  }
  llply(1:length(Bootstraps.boot.cis), function(x) {
    Bootstraps.boot.cis2 <- as.data.frame(Bootstraps.boot.cis[[x]])
    Bootstraps.boot.cis2$variables <- row.names(Bootstraps.boot.cis[[x]])
    row.names(Bootstraps.boot.cis2) <- NULL
    Out <- Bootstraps.boot.cis2[Order, ]
    row.names(Out) <- NULL
    Out[, c(1, 5, 2:4)]
  })
}