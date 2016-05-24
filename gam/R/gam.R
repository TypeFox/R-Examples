"gam" <-
  function(formula, family = gaussian, data, 
           weights, subset, na.action, start = NULL, etastart, mustart, control = gam.control(...),
           model = TRUE, method="glm.fit", x = FALSE, y = TRUE, ...)
{
  call <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mt <- if(missing(data)) terms(formula, gam.slist) else terms(formula,gam.slist,data = data)
  mf$formula<-mt                                                          
  mf <- eval(mf, parent.frame())
   switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1, 
         stop("invalid `method': ", method))


  Y <- model.response(mf, "any")
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0)
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if (!is.null(weights) && any(weights < 0)) 
    stop("Negative wts not allowed")
  if (!is.null(offset) && length(offset) != NROW(Y)) 
    stop("Number of offsets is ", length(offset), ", should equal ", 
         NROW(Y), " (number of observations)")
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
fit<-gam.fit(x=X,y=Y,smooth.frame=mf,weights=weights,start=start,
             etastart=etastart,mustart=mustart,
             offset=offset,family=family,control=control)
  
### If both an offset and intercept are present, iterations are needed to
### compute the Null deviance; these are done here
###
  if(length(offset) && attr(mt, "intercept")>0) {
    fit$null.dev <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
               y = Y, weights = weights, offset = offset, family = family, 
               control = control[c("epsilon","maxit","trace")], intercept = TRUE)$deviance
  }
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if(x) fit$x <- X
    if(!y) fit$y <- NULL
   fit <- c(fit, list(call = call, formula = formula,
		       terms = mt, data = data,
		       offset = offset, control = control, method = method,
		       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("gam","glm", "lm")
  if(!is.null(fit$df.residual) && !(fit$df.residual > 0))
    warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available."
            )
  fit
}

