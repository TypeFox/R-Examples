write.prg.setup <- function(response, object, prg.file, data.file, thismodel, terms.specs)
{
  add.terms <- terms.specs$terms
  bt <- NULL
  vars <- terms.specs$vars
  if(!is.null(object$hlevel)) {
    if(("(Intercept)" %in% vars) || ("Ii11iIInterceptIi11iI" %in% vars))
      bt <- c(bt, "const")
  }
  if(!is.null(add.terms)) {
    for(k in 1L:length(add.terms)) {
      st <- eval(parse(text = add.terms[k]))
      st$term <- rmf(st$term)
      st$by <- rmf(st$by)
      stclass <- class(st)
      if(stclass %in% c("tp.smooth.spec", "cs.smooth.spec")) {
        txt <- paste("constructor class \"", stclass, 
          "\" not supported by bayesx(), using bs = \"ps\" in term ", 
          st$label, "!", sep = "")
        warning(txt, call. = FALSE)
        class(st) <- "ps.smooth.spec"
      }
      bttmp <- bayesx.construct(st, object$outfile, prg.file, object$data)
      if(!is.null(object$hvariables))
        for(k in 1L:length(object$hvariables))
          if(bttmp == paste(object$hvariables[k], "(random)", sep = ""))
            bttmp <- paste(object$hvariables[k], "(hrandom)", sep = "")
      bt <- c(bt, bttmp)	
    }
  }
  bt <- c(bt, terms.specs$linvars)
  if("ModelOffset" %in% vars)
    bt <- c(bt, "ModelOffset(offset)")
  bt <- paste(bt, collapse = " + ")
  if("ModelWeights" %in% vars)
    bt <- paste(bt, "weight", "ModelWeights")
  control.values <- object[attr(object, "co.id")]
  if(object$family == "quantreg")
    control.values$quantile <- object$quantile
  hp <- FALSE
  if(!is.null(object$hlevel) && is.null(object$max.hlevel)) {
    object$max.hlevel <- object$hlevel + 1L
    hp <- TRUE
  }
  if(!is.null(object$hlevel) && object$hlevel != object$max.hlevel && !hp)
    control.values[c("iterations", "burnin", "step", "level1", "level2", "hlevel", "maxint")] <- NULL
  if(hp)
    control.values["hlevel"] <- NULL
  if(!is.null(object$hlevel) && object$hlevel != 1L) {
    # control.values[c("iterations", "burnin", "step", "level1", "level2")] <- NULL
    control.values[c("predict", "modeonly", "setseed", "aresp", "bresp", "pred_check", 
      "mse", "mseparam", "centerlinear", "cv", "quantile", "hlevel", "maxint", "oformula")] <- NULL
  }
  if(object$family == "gaussian_re")
    control.values["predict"] <- NULL
  predict <- control.values$predict
  control.values$predict <- NULL
  control.names <- names(control.values)
  equal <- if(is.null(bt) || bt == "") "" else " = "
  if(is.null(object$hlevel))
    fullformula <- paste("b.regress ", response, equal, bt, ",", sep = "")
  if(!is.null(object$hlevel) || object$hmcmc)
    fullformula <- paste("b.hregress ", response, equal, bt, ",", sep = "")
  for(i in 1L:length(control.values)) {
    if(control.names[i] != "hmcmc")
      fullformula <- paste(fullformula, " ", control.names[i], "=", control.values[[i]], sep = "")
  }
  dset <- paste("d", object$hlevel, sep = "")
  if(!is.null(object$hlevel)) {
    hlevel <- object$hlevel
    if(hlevel > 2L)
      hlevel <- 2L
    fullformula <- paste(fullformula, " hlevel=", hlevel, sep = "")
  }
  if(!is.null(predict)) {
    if(is.logical(predict)) {
      if(predict) {
        if(is.null(object$hlevel))
          fullformula <- paste(fullformula, "predict")
        if((!is.null(object$hlevel) || object$hmcmc) && object$hlevel < 2L)
          fullformula <- paste(fullformula, "predict=full")
      }
    } else {
      fullformula <- paste(fullformula, paste("predict", predict, sep = "="))
    }
  }
  fullformula <- paste(fullformula, "using", dset)
  cat("dataset", dset, "\n", file = prg.file, append = TRUE)
  cat(paste(dset, ".infile using ", data.file, "\n\n", sep = ""), file = prg.file, append = TRUE)
  cat(paste("b.outfile = ", object$outfile, "/", object$model.name, thismodel, "\n\n", sep = ""), 
    file = prg.file, append = TRUE)
  cat(fullformula, "\n\n", file = prg.file, append = TRUE)
  if(object$method == "MCMC" && object$first)
    cat("b.getsample\n\n", file = prg.file, append = TRUE)
  bayesx.prg <- paste(paste(readLines(paste(object$outfile, "/", prg.file, sep=""), n = -1L), 
    collapse = " \n"), " \n", sep="")

  ## only for hierarchical
  if(!is.null(object$hlevel) || object$hmcmc) {
    bayesx.prg <- gsub("psplinerw1", "pspline", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("psplinerw2", "pspline", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("random", "hrandom", bayesx.prg, fixed = TRUE)
    bayesx.prg <- gsub("hhrandom", "hrandom", bayesx.prg, fixed = TRUE)
    cat(bayesx.prg, file = prg.file)
  }

  return(bayesx.prg)
}

