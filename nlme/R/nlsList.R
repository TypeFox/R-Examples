###                  Create a list of nls objects
###
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
### Copyright 2006-2016 The R Core team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

nlsList <-
  ## A list of nls objects
  function(model, data, start, control, level, subset, na.action = na.fail,
           pool = TRUE, warn.nls = NA) # Deprecation: will be 'TRUE'
      UseMethod("nlsList")

nlsList.selfStart <-
  function (model, data, start, control, level, subset, na.action = na.fail,
            pool = TRUE, warn.nls = NA) # Deprecation: will be 'TRUE'
{
  mCall <- as.list(match.call())[-1]
  if (!inherits(data, "groupedData")) {
    stop("second argument must be a groupedData object")
  }
  marg <- substitute(model)
  if (mode(marg) != "name") {
    stop("cannot use an anonymous function for the model")
  }
					# Build up a call to the model function
  m <- call(as.character(marg))
  args <- lapply(names(formals(eval(marg))), as.name)
  args[[1]] <- getCovariateFormula(data)[[2]]
  m[1 + seq_along(args)] <- args
  form <- formula(data)
  form[[3]][[2]] <- m
  mCall$model <- form
  do.call("nlsList.formula", mCall)
}

nlsList.formula <-
  function(model, data, start = NULL, control, level, subset,
           na.action = na.fail, pool = TRUE,
           warn.nls = NA) # Deprecation: will be 'TRUE'
{
  if (!missing(level) && length(level) > 1)
    stop("multiple levels not allowed")
  ## Deprecation: options(show.error.messages = FALSE) should continue to work for now
  if(is.na(warn.nls <- as.logical(warn.nls)))
    warn.nls <- !identical(FALSE, getOption("show.error.messages"))
  Call <- match.call()
  if (!missing(subset)) {
    data <-
      data[eval(asOneSidedFormula(Call[["subset"]])[[2]], data),, drop = FALSE]
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)
  data <- na.action(data)
  if (is.null(grpForm <- getGroupsFormula(model))) {
    if (inherits(data, "groupedData")) {
      if (missing(level))
        level <- length(getGroupsFormula(data, asList = TRUE))
      groups <- getGroups(data, level = level)[drop = TRUE]
      grpForm <- getGroupsFormula(data)
    } else {
      stop("'data' must be a \"groupedData\" object if 'formula' does not include groups")
    }
  } else {
    if (missing(level))
      level <- length(getGroupsFormula(model, asList = TRUE))
    model <- eval(substitute(Y ~ RHS,
			     list(Y  = model[[2]],
				  RHS= getCovariateFormula(model)[[2]])))
    groups <- getGroups(data, form = grpForm, level = level)[drop = TRUE]
  }
  if (is.null(start) && is.null(attr(data, "parameters"))) {
    ## no starting values
    ## checking for old-style selfStart functions
    FUN <- eval(model[[3]][[1]])
    if (is.function(FUN) && class(FUN) != "selfStart" &&
        !is.null(attr(FUN, "initial"))) {
      stop("old-style self-starting model functions\nare no longer supported.\nNew selfStart functions are available.\nUse\n  SSfpl instead of fpl,\n  SSfol instead of first.order.log,\n  SSbiexp instead of biexp,\n  SSlogis instead of logistic.\nIf writing your own selfStart model, see\n  \"help(selfStart)\"\nfor the new form of the \"initial\" attribute.")
    }
  }

  controlvals <- nls.control()
  if(!missing(control)) controlvals[names(control)] <- control
  val <- lapply(split(data, groups),
		function(dat)
                  tryCatch({
                    data <- as.data.frame(dat)
                    if (is.null(start)) {
                      nls(model, data = data, control = controlvals)
                    } else {
                      nls(model, data = data, control = controlvals, start = start)
                    }
                  }, error = function(e) e))
  val <- warnErrList(val, warn = warn.nls)
  if (inherits(data, "groupedData")) {
    ## saving labels and units for plots
    attr(val, "units") <- attr(data, "units")
    attr(val, "labels") <- attr(data, "labels")
    attr(val, "outer") <- attr(data, "outer")
  }

  structure(val, class = c("nlsList", "lmList"),
            call = Call,
            dims = list(N = nrow(data), M = length(val)),
            groups = ordered(groups, levels = names(val)),
            origOrder = match(unique(as.character(groups)), names(val)),
            pool = pool,
            groupsForm = grpForm)
}

###*# Methods for standard generics

coef.summary.nlsList <-
  function(object, ...) object$parameters

formula.nlsList <-
  function(x, ...) eval(attr(x, "call")[["model"]])

summary.nlsList <-
  function(object, ...)
{
  val <- NextMethod("summary") # -> summary.lmList()
  class(val) <- c("summary.nlsList", class(val))
  val
}

update.nlsList <-
    function (object, model., ..., evaluate = TRUE)
{
    call <- attr(object, "call")
    if (is.null(call))
	stop("missing call attribute in \"nlsList\" object")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(model.))
	call$model <- update.formula(formula(object), model.)
    if(length(extras) > 0) {
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

#update.nlsList <-
#  function(object, model, data, start, control, level, subset, na.action,
#	   pool, ...)
#{
#  thisCall <- as.list(match.call())[-(1:2)]
#  if (!missing(model)) {
#    names(thisCall)[match(names(thisCall), "model")] <- "object"
#  }
#  nextCall <- as.list(attr(object, "call")[-1])
#  nextCall[names(thisCall)] <- thisCall
#  do.call("nlsList", nextCall)
#}
