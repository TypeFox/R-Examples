anova.DirichletRegModel <- function(object, ..., sorted = FALSE) {

  comp.objs <- list(object, ...)

  if(!all(unlist(lapply(comp.objs, class)) == "DirichletRegModel")) stop('only models fitted using "DirichReg()" can be compared.')
  for(i in seq_along(comp.objs)[-1L]){
    if(!identical(comp.objs[[i-1L]][["Y"]], comp.objs[[i]][["Y"]])) stop("models appear not to be nested.")
  }

  n.mods <- length(comp.objs)
  n.pars <- unlist(lapply(seq_len(n.mods), function(i){ comp.objs[[i]]$npar }))
  calls <- lapply(seq_len(n.mods), function(i) comp.objs[[i]]$call)
  sorting <- if(sorted){ order(n.pars, decreasing=TRUE) } else { seq_len(n.mods) }

  comp.objs <- comp.objs[sorting]

  deviances <- unlist(lapply(seq_len(n.mods), function(i){ -2*comp.objs[[i]]$logLik }))
  dev_diffs <- abs(c(NA, deviances[1L] - deviances[-1L]))
  df <- abs(c(NA, n.pars[1L] - n.pars[-1L]))

  res <- structure(list(
    "Deviance"   = deviances,
    "N. par"     = n.pars,
    "Difference" = dev_diffs,
    "df"         = df,
    "Pr(>Chi)"   = pchisq(dev_diffs, df, lower.tail = FALSE),
    "sorting"    = sorting,
    "n.mods"     = n.mods,
    "calls"      = calls
  ), class = "anova_DirichletRegModel")

  print(res)
  invisible(res)

}



print.anova_DirichletRegModel <- function(x, ...){

  if(interactive()) writeLines("")

  writeLines("Analysis of Deviance Table\n")

  model_numbers <- unlist(lapply(format(x$sorting), function(i){ paste0("Model ", i) }))
  model_names <- unlist(lapply(x$sorting, function(i){ deparse(x$call[[i]], width.cutoff=500L) }))

  for(i in seq_len(x$n.mods)){
    writeLines(strwrap(paste0(model_numbers[[i]], ": ", model_names[[i]]), width=getOption("width"), exdent = 2L))
  }

  writeLines("")

  res <- as.data.frame(x[1:5])
  colnames(res) <- names(x[1:5])
  rownames(res) <- model_numbers

  printCoefmat(res, cs.ind = 1L, tst.ind = 3L, P.values = TRUE, na.print = "")

  writeLines("")

}
