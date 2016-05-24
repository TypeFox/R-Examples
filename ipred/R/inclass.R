# $Id: inclass.R,v 1.33 2008/08/04 08:18:41 hothorn Exp $

inclass <- function(formula, data, ...) UseMethod("inclass", data)

inclass.default <- function(formula, data,  ...)
{
  stop(paste("Do not know how to handle objects of class", class(data)))
}

inclass.data.frame <- function(formula, data, pFUN = NULL, cFUN = NULL, ...) 
{
##check formula
  if(missing(formula)
    || (length(formula) != 3)
    || (length(attr(terms(formula[-2], data = data), "term.labels")) < 1))
    stop("formula missing or incorrect")

  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)

##editing formula
###main formula
  if(length(formula[[2]])==3) {
    if(is.function(cFUN)) y.formula <- formula[[2]] else y.formula <- cFUN$formula
    w.formula <- XX~YY
    w.formula[[2]] <- formula[[2]][[3]]
    w.formula[[3]] <- formula[[3]]

    response <-  paste(formula[[2]][[2]])
    w.names <- attr(terms(as.formula(formula[[2]]), data = data), "term.labels")
    x.names <- attr(terms(as.formula(formula), data = data), "term.labels")

    if(x.names[1] == ".") x.names <- colnames(data)[!(colnames(data) %in% c(response, w.names))]
  } else {
    stop(paste("Specified formula has to be of type y~x~w"))   
  }


  if(is.null(w.formula)) stop("no formula for prediction model specified")

  formula.list <- vector(mode = "list", length= length(w.names))
  names(formula.list) <- w.names

  P <- length(pFUN)
  Qi <- length(w.names)
  for(j in 1:Qi) {
    res <- list()
    res$formula <- w.formula
    res$formula[[2]] <- as.name(attr(terms(res$formula[-3], data = data), "term.labels")[j])
    if(res$formula[[3]] == ".") {
      res$formula <- as.formula(paste(res$formula[[2]], "~", paste(x.names, collapse= "+")))
    }
    for(i in 1:P) {  
      if(is.null(pFUN[[i]]$formula)) {
        if(is.null(formula.list[[w.names[j]]]$formula)) formula.list[[w.names[j]]]$formula <- res$formula
        if(is.null(formula.list[[w.names[j]]]$model)) formula.list[[w.names[j]]]$model <- pFUN[[i]]$model
        if(is.null(formula.list[[w.names[j]]]$predict)) formula.list[[w.names[j]]]$predict <- pFUN[[i]]$predict
      } else {
        QQ <- attr(terms(pFUN[[i]]$formula[-3], data = data), "term.labels")
        for(k in QQ) {
          if(w.names[j] == k) {
            res$formula[[3]] <- pFUN[[i]]$formula[[3]]
            if(paste(pFUN[[i]]$formula[[3]]) == ".") {
              res$formula[[3]] <- as.formula(paste(w.names[j], "~", paste(x.names, collapse= "+")))
            }
            formula.list[[w.names[j]]]$formula <- pFUN[[i]]$formula
            formula.list[[w.names[j]]]$model <- pFUN[[i]]$model
            formula.list[[w.names[j]]]$predict <- pFUN[[i]]$predict
          }
        }
      }       

    }
  }
                                                                                         
  if(!is.function(cFUN)) {
   cFUN$formula <- y.formula
   if(is.null(cFUN$training.set)) cFUN$training.set <- "original"
  }

  result <- workhorse.inclass(object = formula.list, data = data, cFUN = cFUN, ...)
  return(result)
}


workhorse.inclass <- function(object, data, cFUN, subset, na.action, ...) 
{
  formula.list <- object
  q <- length(formula.list)

  result <- list()
  namen <- c()

##model fitting
  for(i in 1:q) {
    formula <- formula.list[[i]]$formula
##check necessary?? >
    if(missing(formula)
      || (length(formula) != 3)
      || (length(attr(terms(formula[-2], data = data), "term.labels")) < 1)
      || (length(attr(terms(formula[-3], data = data), "term.labels")) != 1))
      stop("formula missing or incorrect")
## check necessary?? < 
    m <- match.call(expand.dots= FALSE)
    res <- formula.list[[i]]$model(formula = formula, data = data)

    namen <- c(namen, as.character(formula[[2]]))
    result <- c(result, list(res))
  }
  names(result) <- namen

  if(!is.function(cFUN)) {
 ###cFUN can be trained on original intermediate variables or on fitted values or on the subset
    if(!is.null(m$subset) && cFUN$training.set == "subset") dataresp <- data[!subset, ]
    if(cFUN$training.set == "original") dataresp <- data
    if(cFUN$training.set == "fitted") {
      dataresp <- data
      for(i in 1:q){
        if(!is.null(formula.list[[i]]$predict)){
         dataresp[,namen[i]] <- formula.list[[i]]$predict(result[[i]], newdata = data)} else {
         dataresp[,namen[i]] <- predict(result[[i]], newdata = data)
        }
      }
    }
    model.response <- cFUN$model(as.formula(cFUN$formula), data = dataresp, ...)
  } else {
    model.response <- cFUN
  }

###predict specificatiations are not delivered
  result <- list("model.intermediate" = result, "model.response" = model.response, "para.intermediate" = object, "para.response" = cFUN)

  class(result) <- "inclass"
  return(result)
}


print.inclass <- function(x, ...)
{
  x <- x$model.intermediate
  q <- length(x)
  intermediates <- attr(x, "names")
  classes <- c()
  for(i in 1:q) {
    classes <- c(classes, class(x[[i]]))
  }

  text.intermediates <- paste("Indirect classification, with", q, "intermediate variables:")
  if(length(unique(classes)) == 1) { 
    predictive  <- paste("Predictive model per intermediate is", unique(classes))
  } else {
    predictive  <- paste("Predictive model per intermediate is \n", 
                        paste(intermediates, ": ", classes, "\n  ", collapse = ""))
  }
  cat("\n", text.intermediates, "\n", intermediates, "\n", "\n", predictive, "\n")
}


summary.inclass <- function(object, ...)
{
  class(object) <- "summary.inclass"
  object
}


print.summary.inclass <- function(x, ...)
{
  x <- x$model.intermediate
  q <- length(x)
  intermediates <- attr(x, "names")
  classes <- c() 
  for(i in 1:q) {
    classes <- c(classes, class(x[[i]]))
  }

  text.intermediates <- paste("Indirect classification, with", q, "intermediate variables:")
  if(length(unique(classes)) == 1) { 
    predictive  <- paste("Predictive model per intermediate is", unique(classes))
  } else {
    predictive  <- paste("Predictive model per intermediate is", "\n ", 
                        paste(intermediates, ": ", classes, "\n  ", collapse = ""))
  }
  cat("\n", text.intermediates, "\n", intermediates, "\n", "\n", predictive,
        "\n", "\n", "Models:", "\n") 
  print(x)

}


