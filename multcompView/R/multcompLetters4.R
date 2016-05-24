"multcompLetters4" <- 
function (object, comp, ...) {
    #Extract needed data from object
    formula <- terms(object)
    Terms <- colnames(attr(terms(object), "factors"))
    data <- model.frame(object)
    fm <- as.character(formula)
    fm <- fm[-1]
    fms <- list()
    for (i in 1:length(Terms)){
      fms[[i]] <- formula(paste(fm[1], "~", Terms[i]))
    }
    names(fms) <- Terms
    if(is.character(comp) | is.symbol(comp)) {
      comp <- match.fun(comp)
      comp <- comp(object)
    }
    comp <- extract_p(comp)
    ans <- list()
    for(i in 1:length(Terms)){
      ans[[i]] <- list(formula = fms[[i]], p = comp[[i]])
    }
    names(ans) <- Terms
    lapply(ans, function(x) multcompLetters2(x$formula, x$p, data, ...))
}

