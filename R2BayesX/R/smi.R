smi <-
function(info.file, model.fit)
{
  info <- NULL
  if(file.exists(info.file))
    if(length(info <- readLines(info.file))) {
      if(is.null(model.fit)) {
        model.fit <- eval(parse(text = info[length(info)]))
      } else {
        model.fit2 <- eval(parse(text = info[length(info)]))
        nmf2 <- names(model.fit2)
        nmf <- names(model.fit)
        model.fit[nmf[nmf %in% nmf2]] <- NULL
        model.fit <- c(model.fit, model.fit2)
      }
    }
  if(!is.null(model.fit$family)) {
    if(model.fit$family == " Gaussian" || model.fit$family == "Gaussian")
      model.fit$family <- "gaussian"
  }
  if(model.fit$method == "REML") {
    model.fit$step <- NULL
    model.fit$iterations <- NULL
  }

  return(model.fit)
}

