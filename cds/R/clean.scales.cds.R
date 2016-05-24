#' @method clean.scales cds
#' @rdname clean.scales
#' @export
clean.scales.cds <- function(object, data, K, col.subset = NULL, ...) {
  if(!inherits(data, 'cdsdata')) data <- createcdsdata(data)
  alphamat <- object$alphamat
  grp <- object$grp
  q <- object$q
  pts <-  (1:q - 0.5)/q
  Mmat <- ispline(pts, tvec = c(0, 0.5, 1), intercept = TRUE)
  opt.scales <- t(tcrossprod(Mmat, alphamat))
  if(!is.null(col.subset)) {
    data <- data$postrs[, col.subset, drop = FALSE]
  } else {
    data <- data$postrs
  }
  tmp <- opt.scales[grp, ]
  clean.ratings <- t(apply(cbind(1:nrow(tmp), tmp), 1, function(x, y) 
                      x[-1][as.numeric(y[x[1], ])], y = data))
  
  dimnames(clean.ratings) <- dimnames(data)

  return(list(clean.data = clean.ratings, optimal.scores = opt.scales))
  
}
