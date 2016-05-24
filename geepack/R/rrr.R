## This is a wrapper function for relative risk regression
## for binary data with log link using the copy method
relRisk<- function(formula, id, waves = NULL,
                   data = parent.frame(), subset = NULL,
                   contrasts = NULL, na.action = na.omit,
		   corstr = "indep",
                   ncopy = 1000, control = geese.control(),
                   b = NULL, alpha = NULL) {
  
  family <- binomial("log") ## fixed
  
  scall <- match.call()
  mnames <- c("", "formula", "data", "offset", "subset", "na.action", "id", "waves")
  cnames <- names(scall)
  cnames <- cnames[match(mnames,cnames,0)]
  mcall <- scall[cnames]
  if (is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())

  y <- model.extract(m, "response")
  if (is.null(dim(y))) N <- length(y) else N <- dim(y)[1]
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts)
  offset <- model.extract(m, "offset")
  if (is.null(offset)) offset <- rep(0, N)

  w <- rep(1 - 1 / ncopy, N)
  w.copy <- rep(1 / ncopy, N)

  y.copy <- 1 - y
  id <- model.extract(m, id)
  waves <- model.extract(m, waves)

  ## augmented data
  Y <- c(y, y.copy)
  W <- c(w, w.copy)
  X <- rbind(x, x)
  ID <- c(id, id + max(id))
  Waves <- c(waves, waves)
  Offset <- c(offset, offset)
  Freq <- rep(c(2, 1), each = N)
  
  ## get initial values
  fit0 <- glm.fit(X, Y, offset = Offset, weights = Freq, family = family)
  fit1 <- glm.fit(X, Y, offset = Offset, family = family, weights = W, start = fit0$coefficients)

  ## feed geese
  ans <- geese.fit(X, Y, ID, Offset, weights = W, waves = Waves,
            family = family, control = control, corstr = corstr,
            b = fit1$coefficients, scale.fix = TRUE)
  ans <- c(ans, list(call=scall, formula=formula)) 
  class(ans) <- "geese"
  ans

}
