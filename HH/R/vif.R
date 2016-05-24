"vif" <-
function(xx, ...)
  UseMethod("vif")

"vif.default" <-
function(xx, y.name, na.action=na.exclude, ...) {
  nnames <- names(xx)
  nn.x <- seq(along=nnames)
  if (missing(y.name))
    y.number <- 0
  else {
    y.number <- match(y.name, nnames, 0)
    nn.x <-  nn.x[-y.number]
  }
  r2 <- nn.x
  names(r2) <- nnames[nn.x]
  if (length(r2) < 2) stop("vif requires two or more X-variables.")
  for (i in nn.x) {
    tmp.lm <- lm(xx[,i] ~
                 data.matrix(xx[,-c(y.number, i)]),
                 na.action=na.action)
  r2[nnames[i]] <- summary(tmp.lm)$r.squared
  }
  1/(1-r2)
}

"vif.formula" <-
function(xx, data, na.action=na.exclude, ...) {
  vif(lm(xx, data, na.action=na.action, x=TRUE))
}

"vif.lm" <-
function(xx, na.action=na.exclude, ...) {
  xxx <- xx  ## deal with scope problem
  if(length(xxx$x)==0 ||
     !(class(xxx$x) == "model.matrix" || class(xxx$x) == "matrix")) {
    xxx <- try(update(xxx, x = TRUE), silent=TRUE)
    if (class(xxx) == "Error" || class(xx)=="try-error") ## S-Plus || R
      stop("Please recompute the 'lm' object with 'x=TRUE'.")
  }
  xx <- as.data.frame(unclass(xxx$x))[-1]
  vif(xx, na.action=na.action)
}
