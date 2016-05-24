"preplot.gam" <-
  function(object, newdata, terms = labels.gam(object),...)
{
  ## this labels.gam above is because there does not seem to be a label method for glms
  Terms <- object$terms
  a <- attributes(Terms)
  Call <- object$call
  all.terms <- labels(Terms)
  xvars <- parse(text=all.terms)
  names(xvars) <- all.terms
  terms <- sapply(terms,match.arg, all.terms)
  Interactions <- a$order > 1
 if(any(Interactions)) {
    all.terms <- all.terms[!Interactions]
    TM <- match(terms, all.terms, 0)
    if(!all(TM)) {
      terms <- terms[TM > 0]
      warning("No terms saved for \"a:b\" style interaction terms"
              )
    }
  }
   xvars <- xvars[terms]
  xnames <- as.list(terms)
  names(xnames) <- terms
  modes <- sapply(xvars, mode)
   for(term in terms[modes != "name"]) {
    evars <- all.names(xvars[term], functions = FALSE, unique = TRUE)
    if(!length(evars))
      next
    xnames[[term]] <- evars
    evars <- parse(text = evars)
    if(length(evars) == 1)
      evars <- evars[[1]]
    else {
      evars <- c(as.name("list"), evars)
      mode(evars) <- "call"
    }
    xvars[[term]] <- evars
  }
  xvars <- c(as.name("list"), xvars)
  mode(xvars) <- "call"

  if(!missing(newdata))
    xvars <- eval(xvars, newdata)
  else {
    if(!is.null(Call$subset) | !is.null(Call$na.action) | !is.null(
                                                                   options("na.action")[[1]])) {
      Rownames <- names(object$fitted)

      if(!(Rl <- length(Rownames)))
        stop("need to have names for fitted.values when call has a subset or na.action argument"
             )
      form<-paste("~",unlist(xnames),collapse="+")
      Mcall <- c(as.name("model.frame"), list(formula = 
                                              terms(as.formula(form)),
                                              subset = Rownames, na.action = function(x)
                                              x))
      mode(Mcall) <- "call"
      Mcall$data <- Call$data
        env <- environment(Terms)##added 7/28/13
        if (is.null(env)) ##
            env <- parent.frame()##
      xvars <- eval(xvars, eval(Mcall,env))
    }
    else {
      ecall <- substitute(eval(expression(xvars)))
      ecall$local <- Call$data
      xvars <- eval(ecall)
    }
  }
  if(missing(newdata))
    pred <- predict(object, type = "terms", terms = terms,
			se.fit = TRUE)
  else pred <- predict(object, newdata, type = "terms", terms = terms,
                           se.fit = TRUE)
  if(is.list(pred)){# oneday predict might return se.fit with newdata
    fits <- pred$fit
    se.fits <- pred$se.fit
  }
  else{
    fits <- pred
    se.fits <- NULL
  }

  gamplot <- xnames
  for(term in terms) {
    x <- xvars[[term]]
    ## oldClass(x) <- unique(c(oldClass(x), data.class(unclass(x))))
    xlab <- xnames[[term]]
    ## Fix ylab for linear terms:
    ylab <- if(length(xlab) == 1 && term == xlab) paste(
                                      "partial for", term) else term
    TT <- list(x = x, y = fits[, term], se.y = if(is.null(se.fits)
                                          ) NULL else se.fits[, term], xlab = xlab, ylab = ylab)
    oldClass(TT) <- "preplot.gam"
    gamplot[[term]] <- TT
  }
  oldClass(gamplot) <- "preplot.gam"
  gamplot
}
