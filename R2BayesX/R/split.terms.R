split.terms <- function(terms, vars, data, dropvars, intcpt)
{
  if(length(terms) > 0L) {
    sm <- is.sm(terms)
    smt <- terms[sm]
    lv <- rmf(terms[!sm])
    if(length(smt) < 1L)
      smt <- NULL
    if(length(lv) < 1L)
      lv <- NULL
    linvars <- NULL
    for(i in 1L:length(lv)) {
      tmp <- gsub("(offset)", "ModelOffset", lv[i], fixed = TRUE)
      tmp <- gsub("(weights)", "ModelWeights", lv[i], fixed = TRUE)
      if(!is.character(data)) {
        tmp2 <- data[[tmp]]
        if(is.factor(tmp2)) {
          if(intcpt) {
            ff <- eval(parse(text = paste("model.matrix(~ ", tmp,")", sep = "")), envir = data)
            tmp <- rmf(colnames(ff)[2L:ncol(ff)])
          } else {
            ff <- eval(parse(text = paste("model.matrix(~ -1 + ", tmp,")", sep = "")), 
              envir = data)
            tmp <- rmf(colnames(ff))
          }
          tmp <- tmp[tmp %in% vars]
        }
        if(is.matrix(tmp2))
          tmp <- paste(tmp, colnames(tmp2), sep = "")
      }
      linvars <- c(linvars, tmp)
    }
    if(!is.null(linvars) && !is.null(dropvars)) {
      for(drop in dropvars)
        linvars <- linvars[linvars != drop]
    }
  } else smt <- linvars <- NULL
  if(!is.null(smt))
    smt <- na.omit(smt)
  if(!is.null(linvars))
    linvars <- na.omit(linvars)
  if(!is.null(vars))
    vars <- na.omit(vars)

  return(list(terms = smt, linvars = linvars, vars = vars))
}

