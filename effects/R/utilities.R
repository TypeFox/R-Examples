# utilities and common functions for effects package
# John Fox, Jangman Hong, and Sanford Weisberg

# 7-25-2013 S. Weisberg modified analyze.model and Analyze.model to ignore
#     default.levels, and use xlevels to set default.  Use grid.pretty by default
# 11-09-2013: fixed error message in Analyze.model(), bug reported by Joris Meys. J. Fox
# 2013-10-15: eliminated functions not needed after effect() methods removed. J. Fox
# 2013-10-29: fixed as.data.frame.*() to handle NA levels. J. Fox
# 2014-03-13: modified Fixup.model.matrix() and Analyze.model() to handle partial residuals; 
#     added is.factor.predictor() and is.numeric.predictor(). J. Fox
# 2014-03-14: error message for non-factor, non-numeric predictor
# 2014-07-08: if no numeric predictor, partial residuals suppressed with warning rather than an error
# 2014-10-09: namespace fixes. J. Fox
# 2015-04-08: added setStrip(), restoreStrip(). J. Fox
# 2015-07-07: fixed matchVarName() so that it handles periods in names properly. J. Fox
# 2015-09-10: added a fix for class = 'array' in Analyze.model.  S. Weisberg
# 2016-02-16: fix Analyze.model(), Fixup.model.matrix() to handle non-focal terms like polynomials correctly; clean up code. J. Fox
# 2016-03-01: correct and improve computation of partial residuals

has.intercept <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

term.names <- function (model, ...) {
  term.names <- gsub(" ", "", labels(terms(model)))
  if (has.intercept(model)) c("(Intercept)", term.names)
  else term.names
}

response.name <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

mfrow <- function(n, max.plots=0){
  # number of rows and columns for array of n plots
  if (max.plots != 0 & n > max.plots)
    stop(paste("number of plots =",n," exceeds maximum =", max.plots))
  rows <- round(sqrt(n))
  cols <- ceiling(n/rows)
  c(rows, cols)
}

expand.model.frame <- function (model, extras, envir = environment(formula(model)),
                                na.expand = FALSE){  # modified version of R base function
  f <- formula(model)
  data <- eval(model$call$data, envir)
  ff <- foo ~ bar + baz
  if (is.call(extras)) 
    gg <- extras
  else gg <- parse(text = paste("~", paste(extras, collapse = "+")))[[1]]
  ff[[2]] <- f[[2]]
  ff[[3]][[2]] <- f[[3]]
  ff[[3]][[3]] <- gg[[2]]
  if (!na.expand) {
    naa <- model$call$na.action
    subset <- model$call$subset
    rval <- if (is.null(data)) eval(call("model.frame", ff, # modified
                                         subset = subset, na.action = naa), envir)           #  lines
    else eval(call("model.frame", ff, data = data,          #
                   subset = subset, na.action = naa), envir)           #
  }
  else {
    subset <- model$call$subset
    rval <- eval(call("model.frame", ff, data = data, subset = subset, 
                      na.action = I), envir)
    oldmf <- model.frame(model)
    keep <- match(rownames(oldmf), rownames(rval))
    rval <- rval[keep, ]
    class(rval) <- "data.frame"
  }
  return(rval)
}

is.relative <- function(term1, term2, factors) {
  all(!(factors[,term1]&(!factors[,term2])))
}

descendants <- function(term, mod, ...){
  names <- term.names(mod)
  if (has.intercept(mod)) names <- names[-1]
  if(length(names)==1) return(NULL)
  which.term <- which(term == names)
  if (length(which.term) == 0){
    factors <- attr(terms(...), "factors")
    rownames(factors) <- gsub(" ", "", rownames(factors))
    colnames(factors) <- gsub(" ", "", colnames(factors))
    (1:length(names))[sapply(names,
                             function(term2) is.relative(term, term2, factors))]
  }
  else {
    factors <- attr(terms(mod), "factors")
    rownames(factors) <- gsub(" ", "", rownames(factors))
    colnames(factors) <- gsub(" ", "", colnames(factors))
    (1:length(names))[-which.term][sapply(names[-which.term],
                                          function(term2) is.relative(term, term2, factors))]
  }
}

is.high.order.term <- function(term, mod,...){
  0 == length(descendants(term, mod, ...))
}

subscripts <- function(index, dims){
  subs <- function(dims, index){
    dim <- length(dims)
    if (dim == 0) return(NULL)
    cum <- c(1,cumprod(dims))[dim]
    i <- index %/% cum
    if (index %% cum != 0) i <- i + 1
    c(i, subs(dims[-dim], index - (i - 1)*cum))
  }
  rev(subs(dims, index))
}

matrix.to.df <- function(matrix, colclasses){
  opt <- options(warn = -1)
  on.exit(options(opt))
  ncol <- ncol(matrix)
  colnames <- colnames(matrix)
  colclasses[sapply(colclasses, function(x) "integer" %in% x)] <- "numeric"
  result <- vector(mode="list", length=ncol)
  names(result) <- colnames
  for (j in 1:ncol){
    result[[j]] <- matrix[, j]
    class <- colclasses[[colnames[j]]]
    result[[colnames[j]]] <- if ("numeric" %in% class) {
      decChar <- getOption('OutDec')
      if (decChar == '.') as.numeric(result[[colnames[j]]])
      else as.numeric(gsub(decChar, '.', matrix[,j]))
    }
    else if ("ordered" %in% class) ordered(result[[colnames[j]]])
    else if ("factor" %in% class) factor(result[[colnames[j]]]) 
    else result[[colnames[j]]]
  }
  as.data.frame(result)
}

# the following function is a modification of code contributed by Steve Taylor

as.data.frame.eff <- function(x, row.names=NULL, optional=TRUE, transform=x$transformation$inverse, ...){
  xx <- x$x
  for (var in names(xx)){
    if (is.factor(xx[[var]])){
      xx[[var]] <- addNA(xx[[var]]) # handle factors with "valid" NA level
    }
  }
  x$x <- xx
  result <- if (is.null(x$se)) data.frame(x$x, fit=transform(x$fit))
  else data.frame(x$x, fit=transform(x$fit), se=x$se, lower=transform(x$lower), upper=transform(x$upper))
  attr(result, "transformation") <- transform
  result
}

as.data.frame.effpoly <- function(x, row.names=NULL, optional=TRUE, ...){
  factors <- sapply(x$variables, function(x) x$is.factor)
  factor.levels <- lapply(x$variables[factors], function(x) x$levels)
  if (!length(factor.levels) == 0){
    factor.names <- names(factor.levels)
    for (fac in factor.names){
      x$x[[fac]] <- factor(x$x[[fac]], levels=factor.levels[[fac]], exclude=NULL)
    }
  }
  result <- data.frame(x$x, x$prob, x$logit)
  if (!is.null(x$confidence.level)) result <- cbind(result,
                                                    x$se.prob, x$se.logit, x$lower.prob, x$upper.prob, x$lower.logit, x$upper.logit)
  result
}

as.data.frame.efflatent <- function(x, row.names=NULL, optional=TRUE, ...){
  xx <- x$x
  for (var in names(xx)){
    if (is.factor(xx$var)){
      xx$var <- addNA(xx$var) # handle factors with "valid" NA level
    }
  }
  x$x <- xx
  if (is.null(x$se)) data.frame(x$x, fit=x$fit)
  else data.frame(x$x, fit=x$fit, se=x$se, lower=x$lower, upper=x$upper)
}

logit2p <- function(logit) 1/(1 + exp(-logit))

p2logit <- function(p) log(p/(1 - p))


lrug <- function(x) {
  if (length(unique(x)) < 0.8 * length(x)) x <- jitter(x)
  grid.segments(x, unit(0, "npc"), x, unit(0.5, "lines"),
                default.units="native")
}

## model.response not generic
model.response.gls <- function(model){
  model.response(model.frame(as.formula(model$call$model), data=eval(model$call$data)))
}

terms.gls <- function(x, ...) terms(formula(x))

## vcov method for eff objects

vcov.eff <- function(object, ...) object$vcov

## [ method for efflist objects

`[.efflist` <- function(x, ...){
  y <- NextMethod("[")
  class(y) <- class(x)
  y
}

### the following functions are for use by Effect() methods

Analyze.model <- function(focal.predictors, mod, xlevels, default.levels=NULL, formula.rhs, 
                          partial.residuals=FALSE, quantiles, x.var=NULL, data=NULL, typical=mean){
  if ((!is.null(mod$nan.action)) && class(mod$na.action) == "exclude") 
    class(mod$na.action) <- "omit"
  all.predictors <- all.vars(formula.rhs)
  check.vars <- !(focal.predictors %in% all.predictors)
  excluded.predictors <- setdiff(all.predictors, focal.predictors)
  number.bad <- sum(check.vars)
  if (any(check.vars)) {
    message <- if (number.bad == 1) paste("the following predictor is not in the model:", 
                                          focal.predictors[check.vars])
    else paste("the following predictors are not in the model:", 
               paste(focal.predictors[check.vars], collapse=", "))
    stop(message)
  }
  X.mod <- model.matrix(mod)
  cnames <- colnames(X.mod)
  factor.cols <- rep(FALSE, length(cnames))
  names(factor.cols) <- cnames
  for (name in all.predictors){
    if (is.factor.predictor(name, mod)) factor.cols[grep(paste("^", name, sep=""), cnames)] <- TRUE
  }
  factor.cols[grep(":", cnames)] <- FALSE   
  X <- na.omit(expand.model.frame(mod, all.predictors))
  bad <- sapply(X[, all.predictors, drop=FALSE], function(x) !(is.factor(x) || is.numeric(x)))
  if (any(bad)){
    message <- if (sum(bad) == 1) paste("the following predictor isn't a factor or numeric:", 
                                        all.predictors[bad])
    else paste("the following predictors aren't factors or numeric:", 
               paste(all.predictors[bad], collapse=", "))
    stop(message)
  }
  x <- list()
  factor.levels <- list()
  if(length(xlevels)==0 & length(default.levels) == 1L) xlevels <- default.levels
  if(is.numeric(xlevels) & length(xlevels) == 1L){
    levs <- xlevels
    for(name in focal.predictors) xlevels[[name]] <- levs
  }
  for (name in focal.predictors){
    levels <- mod$xlevels[[name]]
    if(is.null(levels)) levels <- mod$xlevels[[paste("factor(",name,")",sep="")]]
    fac <- !is.null(levels)
    if (!fac) {    
      levels <- if (is.null(xlevels[[name]])){
        if (partial.residuals){
          quantile(X[, name], quantiles)
        }
        else{
          grid.pretty(range(X[, name]))
        }
      }
      else {
        if(length(xlevels[[name]]) == 1L) { 
          seq(min(X[, name]), max(X[,name]), length=xlevels[[name]])} else
            xlevels[[name]]}
    }
    else factor.levels[[name]] <- levels
    x[[name]] <- list(name=name, is.factor=fac, levels=levels)
  }
  if (partial.residuals){
    numeric.predictors <- sapply(focal.predictors, function(predictor) is.numeric.predictor(predictor, mod))
    if (!any(numeric.predictors)) warning("there are no numeric focal predictors", "\n  partial residuals suppressed")
    else{
      x.var <- which(numeric.predictors)[1]
      x.var.name <- focal.predictors[x.var]
      if (is.null(mod$xlevels[[x.var.name]])){
        x.var.range <- range(X[, focal.predictors[x.var]])
        x[[x.var]][["levels"]] <- seq(from=x.var.range[1], to=x.var.range[2], length=100)
      }
    }
  }
  x.excluded <- list()
  for (name in excluded.predictors){
    levels <- mod$xlevels[[name]]
    fac <- !is.null(levels)
    level <- if (fac) levels[1] else typical(X[, name])
    if (fac) factor.levels[[name]] <- levels
    x.excluded[[name]] <- list(name=name, is.factor=fac,
                               level=level)
  }
  dims <- sapply(x, function(x) length(x$levels))
  len <- prod(dims)
  n.focal <- length(focal.predictors)
  n.excluded <- length(excluded.predictors)
  n.vars <- n.focal + n.excluded
  predict.data <-matrix('', len, n.vars)
  excluded <- sapply(x.excluded, function(x) x$level)
  for (i in 1:len){
    subs <- subscripts(i, dims)
    for (j in 1:n.focal){
      predict.data[i,j] <- x[[j]]$levels[subs[j]]
    }
    if (n.excluded > 0)
      predict.data[i, (n.focal + 1):n.vars] <- excluded
  }
  colnames(predict.data) <- c(sapply(x, function(x) x$name),
                              sapply(x.excluded, function(x) x$name))
  colclasses <- lapply(X, class)
  colclasses[colclasses == "matrix"] <- "numeric"
  colclasses[colclasses == "array"] <- "numeric"
  predict.data <-  matrix.to.df(predict.data, colclasses=colclasses)
  list(predict.data=predict.data, 
       factor.levels=factor.levels, 
       factor.cols=factor.cols, focal.predictors=focal.predictors, n.focal=n.focal,
       excluded.predictors=excluded.predictors, n.excluded=n.excluded,
       x=x, X.mod=X.mod, cnames=cnames, X=X, x.var=x.var)   
}

Fixup.model.matrix <- function(mod, mod.matrix, mod.matrix.all, X.mod,
                               factor.cols, cnames, focal.predictors, excluded.predictors, 
                               typical, given.values){
  attr(mod.matrix, "assign") <- attr(mod.matrix.all, "assign")
  if (length(excluded.predictors) > 0){
    strangers <- Strangers(mod, focal.predictors, excluded.predictors)
    stranger.cols <-  
      apply(outer(strangers, attr(mod.matrix,'assign'), '=='), 2, any)
  }
  else stranger.cols <- rep(FALSE, ncol(mod.matrix))
  if (has.intercept(mod)) stranger.cols[1] <- TRUE
  if (any(stranger.cols)) {
    facs <- factor.cols & stranger.cols
    covs <- (!factor.cols) & stranger.cols
    if (has.intercept(mod)) covs[1] <- FALSE
    if (any(facs)){ 
      mod.matrix[,facs] <-  matrix(apply(as.matrix(X.mod[,facs]), 2, mean), 
                                   nrow=nrow(mod.matrix), ncol=sum(facs), byrow=TRUE)
    }
    if (!is.null(given.values)){
      stranger.names <- cnames[stranger.cols]
      given <- stranger.names %in% names(given.values)
      if (any(given)) {
        mod.matrix[,stranger.names[given]] <- matrix(given.values[stranger.names[given]], nrow=nrow(mod.matrix), 
                                                     ncol=length(stranger.names[given]), byrow=TRUE)
      } 
    }
    for (name in cnames){
      components <- unlist(strsplit(name, ':'))
      components <- components[components %in% cnames]
      if (length(components) > 1) {
        mod.matrix[,name] <- apply(mod.matrix[,components], 1, prod)
      }
    }
  }
  mod.matrix
}

matchVarName <- function(name, expressions){
  scratch <- "zAMIjw4RN3" # randomly generated string
  name <- gsub("\\.", scratch, name)
  expressions <- gsub("\\.", scratch, as.character(expressions))
  a <- !grepl(paste("[.]+", name, sep=""), expressions)
  b <- !grepl(paste(name, "[.]+", sep=""), expressions)
  c <- grepl(paste("\\b", name, "\\b", sep=""), expressions)
  a & b & c
}

Strangers <- function(mod, focal.predictors, excluded.predictors){
  names <- term.names(mod)
  if (has.intercept(mod)) names <- names[-1]
  sel <- apply(sapply(excluded.predictors, matchVarName, expressions=names), 1, any)
  (1:length(sel))[sel]
}

# the following is used by effect.multinom() and Effect.multinom()

eff.mul <- function(x0, B, se, m, p, r, V){
  mu <- exp(x0 %*% B)
  mu <- mu/(1 + sum(mu))
  mu[m] <- 1 - sum(mu)
  logits <- log(mu/(1 - mu))
  if (!se) return(list(p=mu, logits=logits))
  d <- array(0, c(m, m - 1, p))
  exp.x0.B <- as.vector(exp(x0 %*% B))
  sum.exp.x0.B <- sum(exp.x0.B)
  for (j in 1:(m-1)){
    d[m, j,] <- - exp.x0.B[j]*x0
    for (jj in 1:(m-1)){
      d[j, jj,] <- if (jj != j)
        - exp(x0 %*% (B[,jj] + B[,j]))*x0
      else exp.x0.B[j]*(1 + sum.exp.x0.B - exp.x0.B[j])*x0
    }
  }
  d <- d/(1 + sum.exp.x0.B)^2
  V.mu <- rep(0, m)
  for (j in 1:m){
    dd <- as.vector(t(d[j,,]))
    for (s in 1:r){
      for (t in 1:r){
        V.mu[j] <- V.mu[j] + V[s,t]*dd[s]*dd[t]
      }
    }
  }
  V.logits <- V.mu/(mu^2 * (1 - mu)^2)
  list(p=mu, std.err.p=sqrt(V.mu), logits=logits,
       std.error.logits=sqrt(V.logits))
}

# the following are used by effect.polr() and Effect.polr()

eff.polr <- function(x0, b, alpha, V, m, r, se){
  eta0 <- x0 %*% b
  mu <- rep(0, m)
  mu[1] <- 1/(1 + exp(alpha[1] + eta0))
  for (j in 2:(m-1)){
    mu[j] <- exp(eta0)*(exp(alpha[j - 1]) - exp(alpha[j]))/
      ((1 + exp(alpha[j - 1] + eta0))*(1 + exp(alpha[j] + eta0)))
  }
  mu[m] <- 1 - sum(mu)
  logits <- log(mu/(1 - mu))
  if (!se) return(list(p=mu, logits=logits))
  d <- matrix(0, m, r)
  d[1, 1] <- - exp(alpha[1] + eta0)/(1 + exp(alpha[1] + eta0))^2
  d[1, m:r] <- - exp(alpha[1] + eta0)*x0/(1 + exp(alpha[1] + eta0))^2
  for (j in 2:(m-1)){
    d[j, j-1] <- exp(alpha[j-1] + eta0)/(1 + exp(alpha[j-1] + eta0))^2
    d[j, j]   <- - exp(alpha[j] + eta0)/(1 + exp(alpha[j] + eta0))^2
    d[j, m:r] <- exp(eta0)*(exp(alpha[j]) - exp(alpha[j-1]))*
      (exp(alpha[j-1] + alpha[j] + 2*eta0) - 1) * x0 /
      (((1 + exp(alpha[j-1] + eta0))^2)*
         ((1 + exp(alpha[j] + eta0))^2))
  }
  d[m, m-1] <- exp(alpha[m-1] + eta0)/(1 + exp(alpha[m-1] + eta0))^2
  d[m, m:r] <- exp(alpha[m-1] + eta0)*x0/(1 + exp(alpha[m-1] + eta0))^2
  V.mu <- rep(0, m)
  for (j in 1:m){
    dd <- d[j,]
    for (s in 1:r){
      for (t in 1:r){
        V.mu[j] <- V.mu[j] + V[s,t]*dd[s]*dd[t]
      }
    }
  }
  V.logits <- V.mu/(mu^2 * (1 - mu)^2)
  list(p=mu, std.err.p=sqrt(V.mu), logits=logits,
       std.error.logits=sqrt(V.logits))
}

eff.latent <- function(X0, b, V, se){
  eta <- X0 %*% b
  if (!se) return(list(fit=eta))
  var <- diag(X0 %*% V %*% t(X0))
  list(fit=eta, se=sqrt(var))
}

# determine class of a predictor

is.factor.predictor <- function(predictor, model) {
  !is.null(model$xlevels[[predictor]])
}

is.numeric.predictor <- function(predictor, model) {
  is.null(model$xlevels[[predictor]])
}

# manage lattice strips

setStrip <- function(bg=3, fg="black", clip=c("off", "on")){
  clip <- match.arg(clip)
  bg.save <- strip.background <- trellis.par.get("strip.background")
  if (is.numeric(bg) && length(bg) == 1){
    if (bg <= 0) stop("bg should be a positive integer or vector of colors")
    bg <- gray(seq(.95, .5, length=round(bg)))
  }
  strip.background$col <- bg
  fg.save <- strip.shingle <- trellis.par.get("strip.shingle")
  trellis.par.set("strip.background", strip.background)
  if (length(fg) != 1 && length(fg) != length(bg)) 
    stop("lengths of fg and bg incompatible")
  strip.shingle$col <- fg
  trellis.par.set("strip.shingle", strip.shingle)
  clip.save <- .clip <- trellis.par.get("clip")
  .clip$strip <- clip
  trellis.par.set("clip", .clip)
  invisible(list(strip.background=bg.save, strip.shingle=fg.save, clip=clip.save))
}

restoreStrip <- function(saved){
  if (!identical(names(saved), c("strip.background", "strip.shingle", "clip")))
    stop("argument saved does not contain strip parameters")
  trellis.par.set("strip.background", saved$strip.background)
  trellis.par.set("strip.shingle", saved$strip.shingle)
  trellis.par.set("clip", saved$clip)
}