##' safe.as.character
##'
##' @keywords internal
##' @param x x
safe.as.character <- function(x) {
  # preserves dim of x
  xx <-  sapply(as.data.frame(x, stringAsFactor = FALSE, check.names = FALSE), as.character)
  if (is.null(dim(x))) {
    dim(xx) <- NULL
  } else if (nrow(x) == 1) {
    xx <- t(xx)
  } else if (ncol(x) == 1) {
    x <- t(t(xx))
  }
  xx
}

##' trim
##'
##' @keywords internal
##' @param x x
trim <- function (x) {
  x <- sub("^ +", "", x)
  x <- sub(" +$", "", x)
  x
}

##' formatCustom
##'
##' @keywords internal
##' @param x x
##' @param format format
##' @param digits digits
##' @param ... ...
formatCustom <- function(x, format = "nice", digits = 2, ...) {
  sapply(x, function(x) {
    if (!(format %in% c("nice", "plim"))) {
      formatC(x, format = format, digits = digits, ...)
    } else if (format == "nice") {
      formatC(x, format = "f", digits = digits, drop0trailing = TRUE, ...)
    } else {
      plim(x, digits = digits)
    }
  })
}

##' tocharac
##'
##' @export
##' @keywords internal
##' @param x x
##' @param include.rownames include.rownames
##' @param include.colnames include.colnames
##' @param rownames rownames
##' @param colnames colnames
##' @param format format
##' @param digits digits
##' @param decimal.mark decimal.mark
##' @param na.print na.print
tocharac <- function(x, include.rownames = FALSE, include.colnames = FALSE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "") {
  if (is.factor(x))
    x <- as.character(x)
  if (is.vector(x))
    x <- t(x)

  numerics <- sapply(as.data.frame(x), is.numeric)
  
  format <- expand(format, nrow(x), ncol(x), drop = FALSE)
  digits <- expand(digits, nrow(x), ncol(x), drop = FALSE)
  
  rn <- rep(rownames, length = nrow(x))
  cn <- rep(colnames, length = ncol(x))
  if (include.rownames & is.null(rn)) {
    rn <- rownames(x, do.NULL = F, prefix = "")
  }
  if (include.colnames & is.null(cn)) {
    cn <- colnames(x, do.NULL = F, prefix = "")
  }

  xx <- safe.as.character(x)
  
  for (i in 1:ncol(xx)) {
    if (numerics[i]) {
      xx[,i] <- apply(as.matrix(as.numeric(xx[,i])), 2, Vectorize(formatCustom), digits = digits[,i], format = format[,i], decimal.mark = decimal.mark)
    }
    xx[,i] <- trim(xx[, i]) 
  }
  xx[is.na(xx)] <- na.print
  xx[xx == "NA"] <- na.print
  xx[xx == "NaN"] <- na.print
  
  if (include.colnames) {
    xx <- rbind(cn, xx)
  }
  if (include.rownames) {
    if (include.colnames)
      xx <- cbind(c("", rn), xx)
    if (!include.colnames)
      xx <- cbind(rn, xx)
  }

  names(xx) <- NULL
  rownames(xx) <- NULL
  colnames(xx) <- NULL

  xx
}
