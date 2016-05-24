
# This function is used by contrast.lm instead of predictDesign,
# which only works on rms objects.
predictFrame <- function(object, newdata, env, na.action=na.fail)
{
   tt <- tryCatch(terms(object), error=function(e) terms(formula(object, env=env)))
   Terms <- delete.response(tt)

   if (is.null(object$xlevels))
   {
      mf <- model.frame(object, env=env)
      xlev <- .getXlevels(Terms, mf)
   } else {
      xlev <- object$xlevels
   }

   m <- model.frame(Terms, newdata, na.action=na.action, xlev=xlev)
   if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)

   if (is.null(object$contrasts))
   {
      # this is the way gls does it.  is this appropriate for geese?
      contr <- lapply(m, function(el) if (inherits(el, "factor")) contrasts(el))
      contr <- contr[! unlist(lapply(contr, is.null))]
   } else {
      contr = object$contrasts
   }

   model.matrix(Terms, m, contrasts=contr)
}

# This method recreates the model frame from a gls fit object.
# It is needed since the gls function doesn't provide a way to get the
# model frame (like the method="model.frame" option used by lm).
model.frame.gls <- function (formula, ..., env=parent.frame())
{
   Call <- formula$call
   model <- eval(Call$model, env)
   data <- if ('data' %in% names(Call)) eval(Call$data, env) else sys.frame(sys.parent())
   correlation <- if ('correlation' %in% names(Call)) eval(Call$correlation, env) else NULL
   weights <- if ('weights' %in% names(Call)) eval(Call$weights, env) else NULL
   na.action <- if ('na.action' %in% names(Call)) eval(Call$na.action, env) else na.fail

   if (! is.null(correlation))
      groups <- getGroupsFormula(correlation)
   else
      groups <- NULL

   glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))

   mfArgs <- list(formula = asOneFormula(formula(glsSt), model, groups),
                  data = data, na.action = na.action)

   if ('subset' %in% names(Call))
   {
      mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
   }

   mfArgs$drop.unused.levels <- TRUE

   dataMod <- do.call("model.frame", mfArgs)
   if (! is.null(groups))
   {
      groups <- eval(parse(text = paste("~1", deparse(groups[[2]]), sep = "|")))
      grps <- getGroups(dataMod, groups,
                        level = length(getGroupsFormula(groups, asList = TRUE)))
      ord <- order(grps)
      grps <- grps[ord]
      dataMod <- dataMod[ord, ,drop = FALSE]
   }

   model.frame(model, dataMod)
}

model.frame.lme <- function(formula, ..., env=parent.frame())
{
  Call <- formula$call
  fixed <- eval(Call$fixed, env)
  data <- if ('data' %in% names(Call)) eval(Call$data, env) else sys.frame(sys.parent())
  random <- if ('random' %in% names(Call)) eval(Call$random, env) else pdSymm(eval(as.call(fixed[-2])))
  correlation <- if ('correlation' %in% names(Call)) eval(Call$correlation, env) else NULL
  weights <- if ('weights' %in% names(Call)) eval(Call$weights, env) else NULL
  method <- if ('method' %in% names(Call)) eval(Call$method, env) else "REML"
  na.action <- if ('na.action' %in% names(Call)) eval(Call$na.action, env) else na.fail
  contrasts <- if ('contrasts' %in% names(Call)) eval(Call$contrasts, env) else NULL

  method <- match.arg(method, c("REML", "ML"))
  REML <- method == "REML"
  reSt <- reStruct(random, REML = REML, data = NULL)
  groups <- getGroupsFormula(reSt)
  if (is.null(groups)) {
    if (inherits(data, "groupedData")) {
      groups <- getGroupsFormula(data)
      namGrp <- rev(names(getGroupsFormula(data, asList = TRUE)))
      Q <- length(namGrp)
      if (length(reSt) != Q) { # may need to repeat reSt
        if (length(reSt) != 1) {
          stop("Incompatible lengths for \"random\" and grouping factors")
        }
        randL <- vector("list", Q)
        names(randL) <- rev(namGrp)
        for(i in 1:Q) randL[[i]] <- random
        randL <- as.list(randL)
        reSt <- reStruct(randL, REML = REML, data = NULL)
      } else {
        names(reSt) <- namGrp
      }
    } else {
      ## will assume single group
      groups <- ~ 1
      names(reSt) <- "1"
    }
  }
  ## check if corStruct is present and assign groups to its formula,
  ## if necessary
  if (!is.null(correlation)) {
    if(!is.null(corGrpsForm <- getGroupsFormula(correlation, asList = TRUE))) {
      corGrpsForm <- unlist(lapply(corGrpsForm,
                                   function(el) deparse(el[[2]])))
      corQ <- length(corGrpsForm)
      lmeGrpsForm <- unlist(lapply(splitFormula(groups),
                        function(el) deparse(el[[2]])))
      lmeQ <- length(lmeGrpsForm)
      if (corQ <= lmeQ) {
        if (any(corGrpsForm != lmeGrpsForm[1:corQ])) {
          stop(paste("Incompatible formulas for groups in \"random\"",
                     "and \"correlation\""))
        }
        if (corQ < lmeQ) {
          warning(paste("Cannot use smaller level of grouping for",
                        "\"correlation\" than for \"random\". Replacing",
                        "the former with the latter."))
          attr(correlation, "formula") <-
            eval(parse(text = paste("~",
                    deparse(getCovariateFormula(formula(correlation))[[2]]),
                         "|", deparse(groups[[2]]))))
        }
      } else {
        if (any(lmeGrpsForm != corGrpsForm[1:lmeQ])) {
          stop(paste("Incompatible formulas for groups in \"random\"",
                     "and \"correlation\""))
        }
      }
    } else {
      ## using the same grouping as in random
      attr(correlation, "formula") <-
        eval(parse(text = paste("~",
                     deparse(getCovariateFormula(formula(correlation))[[2]]),
                     "|", deparse(groups[[2]]))))
      corQ <- lmeQ <- 1
    }
  } else {
    corQ <- lmeQ <- 1
  }
  ## create an lme structure containing the random effects model and plug-ins
  lmeSt <- lmeStruct(reStruct = reSt, corStruct = correlation,
                     varStruct = varFunc(weights))

  ## extract a data frame with enough information to evaluate
  ## fixed, groups, reStruct, corStruct, and varStruct
  mfArgs <- list(formula = asOneFormula(formula(lmeSt), fixed, groups),
                 data = data, na.action = na.action)
  if ('subset' %in% names(Call))
  {
     mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
  }
  mfArgs$drop.unused.levels <- TRUE
  dataMix <- do.call("model.frame", mfArgs)

  for(i in names(contrasts))            # handle contrasts statement
      contrasts(dataMix[[i]]) = contrasts[[i]]

  ### [fixed], groups, dataMix, correlation, corQ, lmeQ

  grps <- getGroups(dataMix, groups)
  if (inherits(grps, "factor"))
  {
    ord <- order(grps)  #"order" treats a single named argument peculiarly
  } else {
    ord <- do.call("order", grps)
  }
  if (corQ > lmeQ) {
    ## may have to reorder by the correlation groups
    ord <- do.call("order", getGroups(dataMix, getGroupsFormula(correlation)))
  }
  dataMix <- dataMix[ord, ,drop = FALSE]
  model.frame(fixed, dataMix)
}

# This method recreates the model frame from a geese fit object.
# It is needed since the geese function doesn't provide a way to get the
# model frame.
model.frame.geese <- function(formula, ..., env=parent.frame())
{
   scall <- formula$call
   mnames <- c("", "formula", "data", "offset", "weights", "subset", "na.action", "id", "waves", "corp")
   cnames <- names(scall)
   cnames <- cnames[match(mnames, cnames, 0)]
   mcall <- scall[cnames]
   if (is.null(mcall$id)) mcall$id <- as.name("id")
   mcall[[1]] <- as.name("model.frame")
   m <- eval(mcall, env)
   m
}
