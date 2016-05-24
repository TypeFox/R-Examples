## Calculate the residuals from regression of X.j on all the other X variables.
X.residuals <- function(x, ...)
  UseMethod("X.residuals")

X.residuals.default <- function(x, y.name, na.action=na.exclude, ...) {
  nnames <- names(x)
  nn.x <- seq(along=nnames)
  if (missing(y.name))
    y.number <- 0
  else {
    y.number <- match(y.name, nnames, 0)
    nn.x <-  nn.x[-y.number]
  }
  if (length(nn.x) < 2) stop("X.residuals requires two or more X-variables.")
  X.residuals.j <- rep(list(), length=length(nn.x))
  names(X.residuals.j) <- nnames[nn.x]
  for (i in nn.x) {
    tmp.lm <- lm(x[,i] ~
                 data.matrix(x[,-c(y.number, i)]),
                 na.action=na.action)
  X.residuals.j[[nnames[i]]] <- resid(tmp.lm)
  }
  as.data.frame(X.residuals.j)
}

X.residuals.formula <- function(x, data, na.action=na.exclude, ...) {
  X.residuals(lm(x, data, na.action=na.action, x=TRUE))
}

X.residuals.lm <- function(x, na.action=na.exclude, ...) {
  if(length(x[["x"]])==0) { ## '[["x"]]' is needed to avoid partial matching in R
    x <- try(update(x, x = TRUE), silent=TRUE)
    if (class(x) == "Error" || class(x)=="try-error") ## S-Plus || R
      stop("Please recompute the 'lm' object with 'x=TRUE'.")
  }
  x <- as.data.frame(unclass(x[["x"]]))[-1]
  X.residuals(x, na.action=na.action)
}
