ndlVarimp <-
function(object, verbose=TRUE) {

  formula = object$formula
  data = object$data
  response = as.character(formula[2])
  predictors = gsub("[ ]+", " ", paste(deparse(formula[[3]], 
      width.cutoff = 500), collapse = ""))
  n.predictors = length(predictors)
  if (predictors == ".") {
    predictors = colnames(data)
    predictors = predictors[predictors != response]
    formula = as.formula(paste(c(response, paste(predictors, collapse = " + ")), 
            collapse = " ~ "))
  } else {
     predictors = levels(as.factor(unlist(lapply(attr(terms.formula(formula), 
       "term.labels"), function(x) strsplit(x, "[ ]*([\\+]|[\\|]|[:])[ ]*")))))
  }

  if (nlevels(data[,response]) == 2) binary = TRUE
  else binary = FALSE

  if (binary) cvals = rep(0,length(predictors))
  avals = rep(0,length(predictors))
  for (i in 1:length(predictors)) {
    if (verbose) cat(predictors[i],"\n")
    orig = data[,predictors[i]]
    data[,predictors[i]] = sample(data[,predictors[i]])
    tmp.ndl=ndlClassify(formula, data=data)
    data[,predictors[i]] = orig
    if (binary) cvals[i]= summary(tmp.ndl)$statistics$C
    avals[i]= summary(tmp.ndl)$statistics$accuracy
  }
  if (binary) names(cvals) = predictors
  names(avals) = predictors
  
  if (binary) return(list(concordance=cvals, accuracy=avals))
  else return(list(concordance=NA, accuracy=avals))

}
