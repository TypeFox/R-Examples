parse.bayesx.input <- function(formula, data, weights = NULL, subset = NULL, offset = NULL, 
  na.action = na.fail, contrasts = NULL, control = bayesx.control(...), ...)
{
  if(missing(data))
    data <- environment(formula)
  if(is.matrix(data))
    data <- as.data.frame(data)
  if(!is.character(weights))
    weights <- deparse(substitute(weights), backtick = TRUE, width.cutoff = 500L)
  if(!is.character(offset))
    offset <- deparse(substitute(offset), backtick = TRUE, width.cutoff = 500L)
  if(!is.character(subset))
    subset <- deparse(substitute(subset), backtick = TRUE, width.cutoff = 500L)
  if(!is.null(weights) && weights == "NULL") weights <- NULL
  if(!is.null(weights) && weights == "ra$weights") weights <- NULL
  if(!is.null(offset) && offset == "NULL") offset <- NULL
  if(!is.null(offset) && offset == "ra$offset") offset <- NULL
  if(!is.null(subset) && subset == "NULL") subset <- NULL
  if(!is.null(subset) && subset == "ra$subset") subset <- NULL
  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  begin <- NULL 
  co.id <- attr(control, "co.id")
  outfile <- control$outfile
  control$oformula <- formula
  control$terms <- terms(formula, specials = c("s", "te", "r", "sx", "t2"), keep.order = TRUE)
  intcpt <- attr(control$terms, "intercept") > 0
  tlsm <- attr(terms(formula), "term.labels")
  formula <- mgcv::interpret.gam(formula)$fake.formula
  control$first <- TRUE
  fc <- as.character(formula)
  tl <- attr(terms(formula), "term.labels")
  if(length(tl) > 0L && any(is.f(tl))) {
    for(k in 1L:length(tl)) {
      if(is.f(tl[k])) {
        tmp <- eval(parse(text = tl[k]))
        if(!class(tmp) %in% paste(c("rsps", "re", "ra", "random"), "smooth.spec", sep = ".")) {
          tl[k] <- tmp$term[1L]
          if(length(tmp$term) > 1L)
            tl <- c(tl, tmp$term[2L])
          if(tmp$by != "NA")
            tl <- c(tl, tmp$by)
        }
        if(class(tmp) %in% "rsps.smooth.spec")
          control$family <- "gaussian_re"
      }
    }
    if(attr(terms(formula), "intercept") < 1L)
      tl <- c(tl, "-1")
    formula <- as.formula(paste(as.character(formula[2L]), "~", paste(tl, collapse = "+")))
  }
  h.variables <- byra <- NULL
  ra.change <- list()
  if(length(tl) > 0L && any(is.rt(tl))) {
    h.random <- list()
    nr <- 0L
    for(k in 1L:length(tl))
      if(is.rt(tl[k])) {
        tmp <- parse.random.bayesx(tl[k], data)
        tl[k] <- tmp$term
        if(tmp$by != "NA")
          byra <- c(byra, tmp$by)
        if(!is.null(tmp$h.random)) {
          nr <- nr + 1L
          h.random[[nr]] <- tmp$h.random
          h.variables <- c(h.variables, tl[k])
        }
      }
    if(length(h.random) < 1L)
      h.random <- NULL
    else {
      control$h.random <- set.hlevel.stuff(h.random, outfile, control)
      control$h.variables <- h.variables
      control$mcmcreg <- TRUE
      control$hlevel <- 1L
      control$max.hlevel <- get.max.hlevel(control$h.random)
      control$h.random <- set.max.hlevel(control$h.random, control$max.hlevel)
    }
    tl <- unique(c(tl, byra))
    if(!intcpt)
      tl <- c("-1", tl)
    formula <- as.formula(paste(as.character(formula[2L]), "~", paste(tl, collapse = "+")))
    control$ra.change <- ra.change
  } else environment(formula) <- parent.frame()
  control$formula <- formula
  if(!is.character(data)) {
    if(!is.null(weights) && is.character(weights)) {
      W <- data[[weights]]
      if(is.null(W))
        W <- eval(parse(text = weights), envir = .GlobalEnv)
      weights <- W
    }
    if(!is.null(offset) && is.character(offset)) {
      O <- data[[offset]]
      if(is.null(O))
        O <- eval(parse(text = offset), envir = .GlobalEnv)
      offset <- O
    }
    if(!is.null(subset) && is.character(subset)) {
      S <- data[[subset]]
      if(is.null(S)) {
        S <- try(eval(parse(text = subset), envir = .GlobalEnv), silent = TRUE)
        if(class(S) == "try-error")
          S <- try(eval(parse(text = subset), envir = data), silent = TRUE)
        if(class(S) == "try-error")
          stop("problems evaluating argument subset!")
      }
      subset <- S
    }
    if(!is.null(control$begin)) {
      begin <- data[[control$begin]]
      if(is.null(begin)) {
        begin <- try(eval(parse(text = control$begin), envir = .GlobalEnv), silent = TRUE)
        if(class(begin) == "try-error")
          begin <- try(eval(parse(text = control$begin), envir = data), silent = TRUE)
        if(class(begin) == "try-error")
          stop("problems evaluating argument begin!")
      }
    }
    ff <- formula
    Yn <- as.character(ff[2L])
    Y <- eval(parse(text = Yn), envir = data)
    if(is.factor(Y)) {
      control$YLevels <- levels(Y)
      Y <- f2int(Y)
      control$nYLevels <- levels(as.factor(Y))
      if(!is.null(control$reference)) {
         aaa <- 1 ## FIXME!
      }
    }
    if(!is.null(control$reference)) {
      if(!is.null(control$YLevels)) {
        if(!is.character(control$reference))
          control$reference <- control$nYLevels[control$reference]
        else
          control$reference <- control$nYLevels[control$YLevels == control$reference]
      }
      if(!(control$reference %in% control$nYLevels))
        stop("argument reference is specified wrong, level not within response levels!")
    }
    ff[2] <- NULL
    only <- only2 <- FALSE
    if(resplit(as.character(ff)) == "~1") {
      ff <- as.formula(paste(Yn, "~1"))
      only2 <- TRUE
    }
    if(resplit(as.character(ff)) == "~-1") {
      data <- as.data.frame(Y)
      names(data) <- Yn 
      only <- TRUE
    } else {
      if(only || only2) {
        if(length(Y) > length(weights))
          weights <- weights[get.unique(Y, 22L)$ind]
        if(length(Y) > length(offset))
          offset <- offset[get.unique(Y, 22L)$ind]
      }
      if(is.function(weights)) weights <- NULL
      if(is.function(subset)) subset <- NULL
      if(is.function(offset)) offset <- NULL
      ff2 <- eval(parse(text = paste("update(ff,", Yn, "~ .)")))
      data[[Yn]] <- Y
      ml <- list(formula = ff2, data = data, weights = if(control$prediction) NULL else weights,
        subset = subset, offset = offset, na.action = na.action, drop.unused.levels = TRUE)
      data <- do.call("model.frame", ml)
      if(control$prediction) {
        if(!is.null(weights))
          data[["(weights)"]] <- weights
        if(!is.null(offset))
          data[["(offset)"]] <- offset
        if(!is.null(subset))
          data <- subset(data, subset)
      }
      if(!is.null(h.variables)) {
        nd <- names(data)
        for(k in nd)
          if(k %in% h.variables) {
            if(is.factor(data[[k]]))
              data[[k]] <- f2int(data[[k]])
          }
      }
      if(any(tlsm2 <- is.sm(tlsm))) {
        for(k in 1L:length(tlsm))
          if(tlsm2[k]) {
            object <- eval(parse(text = tlsm[k]))
            term <- object$term
            for(j in 1L:length(term))
              if(is.factor(data[[term[j]]])) {
                if(class(object) %in% c("mrf.smooth.spec", "gk.smooth.spec", "gs.smooth.spec")) {
                  data[[term[j]]] <- f2int(data[[term[j]]], type = 3L)
                } else {
                  data[[term[j]]] <- f2int(data[[term[j]]])
                }
              }
          }
      }
      if(ncol(data) < 2L && names(data) == Yn)
        only <- TRUE
    }
#    if(!is.null(subset))
#      Y <- Y[subset]
#    if(length(Y) != nrow(data)) {
#      Y <- unique(Y)
#      if(length(Y) != nrow(data))
#        warning(paste("variable lengths differ (found for \'",Yn,"\')", sep = ""))
#    }
#    if(!only) {
#      Y <- as.data.frame(Y)
#      names(Y) <- Yn
#      data <- cbind(Y, data)
#    }
    if(ncol(data) < 2L) {
      control$order <- order(data[, 1L])
      data[, 1L] <- data[order(data[, 1L]), 1L]
    } else {
      control$order <- order(Y)
      data <- data[order(Y), ]
    }
    control <- c(control, list(data = na.action(data), Yn = Yn))
  } else {
    Yn <- as.character(formula[2L])
    data <- path.expand(data)
    control <- c(control, list(data = data, Y = Yn, Yn = Yn, weights = weights, offset = offset))
  }
  attr(control$data, "terms") <- control$formula
  attr(attr(control$data, "terms"), "response") <- 1L
  attr(control$data, "na.action") <- na.action
  control$contrasts <- contrasts
  if(!is.null(begin))
    control$begin.vec <- begin
  attr(control, "co.id") <- co.id
  class(control) <- "bayesx.input"

  return(control)
}


get.max.hlevel <- function(x)
{
  hlevel <- NULL
  for(k in 1L:length(x)) {
    if(!is.null(x[[k]]$hlevel)) {
      hlevel <- max(c(hlevel, x[[k]]$hlevel))
    if(!is.null(x[[k]]$h.random))
      hlevel <- max(c(hlevel, get.max.hlevel(x[[k]]$h.random)))
    }
  }

  return(hlevel)
}


set.max.hlevel <- function(x, max.hlevel)
{
  for(k in 1L:length(x)) {
    if(!is.null(x[[k]]$hlevel))
      x[[k]]$max.hlevel <- max.hlevel
    if(!is.null(x[[k]]$h.random))
      x[[k]]$h.random <- set.max.hlevel(x[[k]]$h.random, max.hlevel)
  }

  return(x)
}
