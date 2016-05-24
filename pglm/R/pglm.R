ordinal <- function(link = c('probit', 'logit')){
  link <- match.arg(link)
  list(family = 'ordinal', link = link)
}

negbin <- function(link = c('log'), vlink = c('nb1', 'nb2')){
  link <- match.arg(link)
  vlink <- match.arg(vlink)
  list(family = 'negbin', link = link, vlink = vlink)
}

pglm <-  function(formula, data, subset, na.action,
                  effect = c('individual','time','twoways'),
                  model  = c('random', 'pooling', 'within', 'between'),
                  family, other = NULL, index  = NULL, start = NULL, R = 20, ...){
  dots <- list(...)
  args <- list(model = model, effect = effect)

  if (is.character(family)){
    if (family %in% c("ordprobit", "ordlogit", "tobit")){
      if (family == "ordprobit") family <- list(family = "ordinal", link = "probit")
      if (family == "ordlogit") family <- list(family = "ordinal", link = "logit")
      if (family == "tobit") family <- list(family = "tobit", link = NULL)
    }
    else  family <- get(family, mode = "function")
  }
  if (is.function(family)) family <- family()
  
  link <- thelink <- family$link
  if (family$family == "negbin") vlink <- family$vlink
  family <- family$family
  
  # check and match the arguments
  effect <- match.arg(effect)
  if (!any(is.na(model))) model <- match.arg(model)

  # Check whether data is a pdata.frame and if not create it ; ignore
  # this step if model = pooling and if index = NULL so that ordinary
  # glm models can be fitted
  if (model != "pooling" | !is.null(index)){
    if (inherits(data, "pdata.frame") && !is.null(index))
      warning("the index argument is ignored because data is a pdata.frame")
    if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index)
    # Create a Formula object if necessary
    if (!inherits(formula, "pFormula")) formula <- pFormula(formula)
  }
  
  # eval the model.frame
  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "na.action"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
  mf$data <- data
  data <- eval(mf, parent.frame())
  # return the model.frame or estimate the model
  if (is.na(model)){
    attr(data, "formula") <- formula
    return(data)
  }
  Kw <- NULL
  if (model != "pooling"){
    X <- model.matrix(formula, data, rhs = 1, model = "pooling", effect = effect)
    if (model == "within" && family == "poisson"){
      Xw <- model.matrix(formula, data, rhs = 1, model = "within", effect = effect)
      Kw <- colnames(Xw)
      X <- X[, Kw, drop = FALSE]
    }
    if (ncol(X) == 0) stop("empty model")
    y <- pmodel.response(formula, data, model = "pooling", effect = effect)
    id <- attr(data, "index")[[1]]
  }
  else{
    X <- model.matrix(formula, data)
    y <- model.response(data)
    if (inherits(data, "pdata.frame")) id <- attr(data, "index")[[1]]
    else id <- NULL
  }

  # compute the nodes and the weights for the gaussian quadrature
  if (model == "random" && (! family %in% c("poisson", "negbin", "gaussian")))
      rn <- gauss.quad(R, kind = 'hermite')
  else rn <- NULL

  # compute the starting vgalues
  start <- starting.values(family, link, vlink, rn, model, Kw, X, y, id, cl, start, other)


  # call to maxLik with the relevant arguments
  ml <- cl
  m <- match(c("print.level", "ftol", "tol", "reltol",
               "gradtol", "steptol", "lambdatol", "qrtol",
               "iterlim", "fixed", "activePar", "method"),names(ml),0)
  ml <- ml[c(1, m)]

  argschar <- function(args){
    paste(as.character(names(args)), as.character(args),
          sep="=", collapse=",")
  }

  args <- list(param = "start",
               y = "y", X = "X", id = "id", model = "model", link = "link",
               rn = "rn")
  if (family %in% c("tobit", "gaussian", "poisson")) args$other <- "other"
  if (family == "negbin") args$vlink <- "vlink"

  thefunc <- paste("function(start) lnl.", family,
                   "(", argschar(args), ")", sep = "")
  ml$logLik <- eval(parse(text = thefunc))
  thefunc <- paste("function(start) attr(lnl.", family,
                   "(", argschar(args), "), \"gradient\")", sep = "")
#  ml$grad <- eval(parse(text = thefunc))
  thefunc <- paste("function(start) attr(lnl.", family,
                   "(", argschar(args), "), \"hessian\")", sep = "")
#  ml$hess <- eval(parse(text = thefunc))
  
  ml$start <- start
  ml[[1]] <- as.name('maxLik')
#  if (family == "negbin") ml$activePar <- 14
  result <- eval(ml, parent.frame())
  result[c('call', 'args', 'model')] <- list(cl, args, data)
  result
}

