## prediction method based on refitting with weights
predict.bayesx <- function(object, newdata, model = NULL,
  type = c("response", "link", "terms", "model"), na.action = na.pass,
  digits = 5, ...)
{
  if(!missing(newdata)) {
    object <- get.model(object, model)
    if(length(object) > 1)
      stop("argument model is specified wrong, predictions are only possible one single model!")

    object <- object[[1]]
    mf <- model.frame(object)
    stopifnot(inherits(newdata, c("data.frame", "list", "matrix")))
    newdata <- as.data.frame(newdata)

    ff <- object$formula
    response <- all.vars(ff)[1]
    yf <- is.factor(mf[[response]])
    newdata[[response]] <- if(yf) mf[[response]][1] else 0
    control <- object$bayesx.setup
    control$reference <- NULL
    newdata <- model.frame.bayesx(ff, data = newdata, na.action = na.action)
    newdata[[response]] <- NULL
    nam_nd <- names(newdata)
    nam_mf <- names(mf)
    nam_mf <- nam_mf[nam_mf != response]
    if(!all(nc <- nam_mf %in% nam_nd))
      stop(paste("variables", paste(nam_mf[!nc], collapse = ", "), "are missing in newdata!"))

    nd <- list()
    for(j in nam_mf) {
      nd[[j]] <- unlist(list(mf[[j]], newdata[[j]]))
    }
    nd <- as.data.frame(nd)
    names(nd) <- nam_mf
    cr <- rep(if(yf) mf[[response]][1] else 0, length = nrow(newdata))
    nd[[response]] <- unlist(list(mf[[response]], cr))
    weights <- model.weights(mf)
    if(is.null(weights))
      weights <- rep(1, length = nrow(mf))
    i <- c(rep(FALSE, length(weights)), rep(TRUE, nrow(newdata)))
    nd$weights <- weights <- c(weights, rep(0, length = nrow(newdata)))
    oterms <- names(object$effects)
    object <- update(object, . ~ ., data = nd, weights = weights,
      seed = object$bayesx.setup$setseed, prediction = TRUE, ...)
  } else {
    newdata <- model.frame(object)
    oterms <- names(object$effects)
    i <- rep(TRUE, nrow(newdata))
  }

  type <- match.arg(type)
  if(type == "model") return(object)
  if(type %in% c("link", "response")) {
    pr <- fitted(object)
    if(!is.null(dim(pr))) {
      if(any(j <- grepl(me <- if(type == "link") "eta" else "mu", names(pr)))) {
        pr <- pr[i, j, drop = FALSE]
        colnames(pr) <- gsub(paste(me, ":", sep = ""), "", colnames(pr), fixed = TRUE)
      } else
        pr <- pr[i, 1]
    } else pr <- pr[i]
  }
  if(type == "terms") {
    pr <- fitted(object, term = oterms)
    if(inherits(pr, "data.frame")) pr <- list(pr)
    labels <- NULL
    for(j in seq_along(pr)) {
      if(inherits(pr[[j]], "data.frame")) {
        nt <- names(pr[[j]])[1]
        i <- round(pr[[j]][[nt]], digits) %in% round(newdata[[nt]], digits)
        pr[[j]] <- pr[[j]][i, ]
        rownames(pr[[j]]) <- NULL
      }
      tl <- attr(pr[[j]], "specs")$label
      labels <- c(labels, if(is.null(tl)) "NA" else tl)
    }
    names(pr) <- labels
    if(length(pr) < 2) pr <- pr[[1]]
  }
  
  return(pr)
}

