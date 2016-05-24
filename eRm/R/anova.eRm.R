anova.eRm <- function(object, ...){
  models <- c(list(object), list(...))
#browser()
  # exclude LLRAs
  if(any(unlist(lapply(models, function(m){ "llra" %in% class(m) })))) stop("At least one model is an LLRA; comparison to other models not possible.")

  # check if models' data matrices are identical
  for(i in seq_along(models)[-1L]){
    if(!identical(unname(models[[1L]][["X"]]), unname(models[[i]][["X"]]))) stop("Models are not nested.")
  }

  # sort by number of parameters
  models <- models[order(unlist(lapply(models, function(m){ m[["npar"]] })), decreasing = TRUE)]

  # extract information
  calls <- unlist(lapply(models, function(m){ deparse(m[["call"]]) }))
  LLs   <- unlist(lapply(models, function(m){ m[["loglik"]] }))
  npar  <- unlist(lapply(models, function(m){ m[["npar"]] }))
  dev   <- -2*LLs
  LR    <- abs(c(NA, LLs[1L] - LLs[-1L]))
  df    <- c(NA, npar[1L] - npar[-1L])
  p     <- pchisq(LR, df, lower.tail = FALSE)

  return(
    structure(
      list(
        calls      = calls,
        statistics = data.frame(LLs=LLs, dev=dev, npar=npar, LR=LR, df=df, p=p)
      ),
    class="eRm_anova")
  )

}

print.eRm_anova <- function(x, ...){
  if(interactive()) writeLines("")
  writeLines("Analysis of Deviances Table\n")

  for(i in seq_along(x[[1L]])){
    writeLines(strwrap(paste0("Model ", i, ": ", x[[1L]][[i]]), width = getOption("width"), exdent = 4L))
  }
  writeLines("")

  x_print <- as.matrix(x[[2L]])
  rownames(x_print) <- paste0("Model ", seq_along(x[[1L]]))
  colnames(x_print) <- c("cond. LL", "Deviance", "npar", "LR", "df", "p-value")
  printCoefmat(as.matrix(x_print), cs.ind=c(1,2), tst.ind=4, has.Pvalue=TRUE, na.print = "")

  writeLines("")
  message("Note: The models appear to be nested, please check this assumption.")
  if(interactive()) writeLines("")

  invisible(x)
}
