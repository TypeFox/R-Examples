setMethod("rflexmix", signature(object = "FLXdist", newdata="numeric"), function(object, newdata, ...) {
  newdata <- data.frame(matrix(nrow = as.integer(newdata), ncol = 0))
  rflexmix(object, newdata = newdata, ...)
})

setMethod("rflexmix", signature(object = "FLXdist", newdata="listOrdata.frame"), function(object, newdata, ...) {
  groups <- .FLXgetGrouping(object@formula, newdata)
  object@model <- lapply(object@model, FLXgetModelmatrix, newdata, object@formula, lhs=FALSE)
  group <- if (length(groups$group)) groups$group else factor(seq_len(FLXgetObs(object@model[[1]])))
  object@concomitant <- FLXgetModelmatrix(object@concomitant, data = newdata,
                                                       groups = list(group=group,
                                                         groupfirst = groupFirst(group)))
  rflexmix(new("flexmix", object, group=group, weights = NULL), ...)
})

setMethod("rflexmix", signature(object = "flexmix", newdata="missing"), function(object, newdata, ...) {
  N <- length(object@model)
  object <- undo_weights(object)
  group <- group(object)
  prior <- determinePrior(object@prior, object@concomitant, group)
  class <- apply(prior, 1, function(x) stats::rmultinom(1, size = 1, prob = x))
  class <- if (is.matrix(class)) t(class) else as.matrix(class)
  class <- max.col(class)[group]
  y <- vector("list", N)
  for (i in seq_len(N)) {
    comp <- lapply(object@components, function(x) x[[i]])
    yi <- rFLXM(object@model[[i]], comp, class, group, ...)
    form <- object@model[[i]]@fullformula
    names <- if(length(form) == 3) form[[2]] else paste("y", i, seq_len(ncol(yi)), sep = ".")
    if (ncol(yi) > 1) {
      if (inherits(names, "call")) 
        names <- as.character(names[-1])
      if (length(names) != ncol(yi)) {
        if (length(names) == 1) names <- paste(as.character(names)[1], i, seq_len(ncol(yi)), sep = ".")
        else stop("left hand side not specified correctly")
      }
    }
    else if (inherits(names, "call")) names <- deparse(names)
    colnames(yi) <- as.character(names)
    y[[i]] <- yi
  }
  list(y = y, group=group, class = class)
})

###**********************************************************

determinePrior <- function(prior, concomitant, group) {
  matrix(prior, nrow = length(unique(group)), ncol = length(prior), byrow=TRUE)
}

setGeneric("determinePrior", function(prior, concomitant, group)
           standardGeneric("determinePrior"))

setMethod("determinePrior", signature(concomitant="FLXPmultinom"), function(prior, concomitant, group) {
  exps <- exp(concomitant@x %*% concomitant@coef)
  exps/rowSums(exps)
})

undo_weights <- function(object) {
  if (!is.null(object@weights)) {
    for (i in seq_along(object@model)) {
      object@model[[i]]@x <- apply(object@model[[i]]@x, 2, rep, object@weights)
      object@model[[i]]@y <- apply(object@model[[i]]@y, 2, rep, object@weights)
      object@concomitant@x <- apply(object@concomitant@x, 2, rep, object@weights)
    }
    if (length(object@group) > 0) 
      object@group <- rep(object@group, object@weights)
    object@weights <- NULL
  }
  object
}

###**********************************************************

setMethod("simulate", signature("FLXdist"),
function(object, nsim, seed = NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv, inherits = FALSE))
  }
  ans <- lapply(seq_len(nsim), function(i) rflexmix(object, ...)$y)
  if (all(sapply(ans, ncol) == 1)) ans <- as.data.frame(ans)
  attr(ans, "seed") <- RNGstate
  ans
})
