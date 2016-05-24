DynNom <- function(model, data,
                   clevel = 0.95, covariate = c("slider", "numeric"),
                   ptype = c("st", "1-st")) {

  data <- data.frame(data)

  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if (attr(model$terms, "dataClasses")[[1]] == "logical")
    stop("Error in model syntax: logical form for response not supported")

  if (tail(names(attr(model$terms,"dataClasses")),n=1)=="(weights)") {
    n.terms <- length(attr(model$terms,"dataClasses"))
    attr(model$terms,"dataClasses") <- attr(model$terms,"dataClasses")[1:n.terms - 1]
  }

  if (attr(model, "class")[1] == "lm"|
      attr(model, "class")[1] == "glm") {
    for(i in 1:length(names(attr(model$terms, "dataClasses")))) {
      com1=numeric(length(names(data)))
      for(j in 1:length(names(data))) {
        if (names(attr(model$terms, "dataClasses"))[i]==names(data)[j]) com1[j]=1
      }
      if (sum(com1)==0)
        stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
    }
  }

  if (attr(model, "class")[1] == "coxph.null") {
    stop("Error in model syntax: the model is null")
  }

  if (attr(model, "class")[1] == "coxph") {
    n.strata <- length(attr(model$terms, "specials")$strata)
    dim.terms <- length(names(attr(model$terms, "dataClasses")))

    for (i in 2:dim.terms) {
      if (substr(names(attr(model$terms, "dataClasses"))[i], 1, 6) == "strata") {
        nch <- nchar(names(attr(model$terms, "dataClasses"))[i])
        names(attr(model$terms, "dataClasses"))[i] <- substr(names(attr(model$terms,
                                                                        "dataClasses"))[i], 8, (nch - 1))
      }
    }

    if (!is.null(attr(model$terms, "specials")$tt)) {
      stop("Error in model syntax: coxph models with a time dependent covariate is not supported")
    }

    for(i in 2:length(names(attr(model$terms, "dataClasses")))) {
      com1=numeric(length(names(data)))
      for(j in 1:length(names(data))) {
        if (names(attr(model$terms, "dataClasses"))[i]==names(data)[j]) com1[j]=1
      }
      if (sum(com1)==0)
        stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
    }
  }

  if (attr(model, "class")[1] == "lm") {
    DynNom.lm(model, data, clevel, covariate)
  }
  if (attr(model, "class")[1] == "glm") {
    DynNom.glm(model, data, clevel, covariate)
  }
  if (attr(model, "class")[1] == "coxph") {
    if (attr(model$terms, "dataClasses")[[1]] == "nmatrix.3")
      stop("Error in model syntax: start/stop notation not supported")

    if (attr(model$terms, "dataClasses")[[1]] == "nmatrix.2") {
      DynNom.coxph(model, data, clevel, covariate, ptype)
    }
  }
}
