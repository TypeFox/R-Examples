r <-
function(x, h = NULL, by = NA, xt = NULL, 
  data = NULL, weights = NULL, subset = NULL, 
  offset = NULL, na.action = na.fail, contrasts = NULL, 
  control = bayesx.control(...), ...)
{
  term <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500L)
  call <- match.call()
  is.formula <- FALSE
  if(!any(grepl("~", term)) && is.null(h) && !is.null(data)) {
    term <- paste(term, "~", 1)
  } 
  if(any(grepl("~", term))) {
    tmp <- strsplit(term, "~")[[1L]]
    call$h <- h <- term
    term <- splitme(tmp[1L])
    term <- resplit(term[term != " "])
    call$x <- term
    is.formula <- TRUE
  }
  by.var <- if(!is.character(by)) deparse(substitute(by), backtick = TRUE, width.cutoff = 500L) else by
  ins <- formula <- NULL
  if(by.var == ".") 
    stop("by=. not allowed")
  if(term == ".") 
    stop("r(.) not yet supported.")
  label <- paste("r(", term)
  if(!is.null(h)) {
    ins <- list()
    mlabel <- paste(as.character(call$h), collapse = " ")
    split <- splitme(mlabel)
    if(split[1L] != "~" && !is.formula)
      mlabel <- resplit(c("~", split))
    if(!is.formula)
      formula <- as.formula(paste(term, mlabel))
    else
      formula <- as.formula(mlabel)
    if(length(grep("~", mlabel, fixed = TRUE)))
      label <- paste("r(", mlabel)
    else
      label <- paste(label, ",", mlabel, collapse="")
    mf <- terms.formula(formula, specials=c("sx", "s", "te", "r"))
    mterms <- attr(mf, "term.labels")
    if(length(mterms) > 0L)
      for(k in 1L:length(mterms)) {
        if(is.sm(mterms[k]))
          ins[[k]] <- try(eval(parse(text = mterms[k])), silent = TRUE)
        else {
          ins[[k]] <-list(term = mterms[k], label = mterms[k])
          class(ins[[k]]) <- "lin.smooth.spec"
        }
      }
    }
  if(by.var != "NA")
    label <- paste(label, ",by=", by.var, collapse = "")
  label <- gsub(" ", "", paste(label, ")", sep = ""))
  rval <- list(term = term, label = label, by = by.var, xt = xt, 
    ins = ins, formula = formula, data = data, weights = weights, 
    subset = subset, offset = offset, na.action = na.action, 
    contrasts = contrasts, control = control)
  if(!is.null(control$bs) && control$bs == "rsps") {
    rval$control$bs <- NULL
    class(rval) <- "rsps.smooth.spec"
    if(is.null(rval$ins)) {
      rval$ins <- list()
      rval$formula <- as.formula(paste(rval$term, "~ -1"))
    }
  } else class(rval) <- "ra.smooth.spec"

  return(rval) 
}

