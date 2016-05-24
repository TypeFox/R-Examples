ndlClassify <- function(formula, data, frequency=NA, variable.value.separator="", ...)
{
  call <- match.call()

#  response = as.character(formula[2])
  response=gsub("[ ]+"," ",paste(deparse(formula[[2]],width.cutoff=500),collapse=""))
  predictors=gsub("[ ]+"," ",paste(deparse(formula[[3]],width.cutoff=500),collapse=""))

  if(predictors == ".")
    { predictors = colnames(data)
      predictors = predictors[predictors!=response]
      if(!is.na(frequency) & is.character(frequency) & length(frequency)==1)
        predictors = predictors[predictors!=frequency]
      formula = as.formula(paste(c(response,paste(predictors,collapse=" + ")),collapse=" ~ "))
    } 
  else
    predictors = levels(as.factor(unlist(lapply(attr(terms.formula(formula),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\|]|[:])[ ]*")))))

  if(is.character(frequency) & length(frequency)==1 & frequency %in% colnames(data))
    data <- data[c(frequency, response, predictors)]
  else
    data <- data[c(response, predictors)]

  cuesOutcomes = ndlCuesOutcomes(formula=formula, data=data, frequency=frequency, variable.value.separator=variable.value.separator, ...)

  weightMatrix = estimateWeights(cuesOutcomes, ...)
  weightMatrix = weightMatrix[order(rownames(weightMatrix)),,drop=FALSE]
  weightMatrix = weightMatrix[,order(colnames(weightMatrix)),drop=FALSE]

  activationMatrix = estimateActivations(cuesOutcomes, weightMatrix, ...)$activationMatrix 

  result <- list(activationMatrix=activationMatrix,  weightMatrix=weightMatrix, cuesOutcomes=cuesOutcomes, frequency=frequency, call=call, formula=formula, data=data)
  class(result) <- "ndlClassify"

  return(result)
}

print.ndlClassify <- function(x, max.print=10, ...)
{
  digits=max(3,getOption("digits")-3)
  if(is.na(max.print))
    max.print=NROW(x$weightMatrix)
#  if(!is.null(x$max.print) & is.numeric(x$max.print))
#       max.print=x$max.print
  cat("\n")
  print(x$call)
  cat("\n")
  print(x$formula)
  cat("\n")
  tabl <- x$weightMatrix[1:min(nrow(x$weightMatrix),max.print),]
  class(tabl) <- "table"
  print(tabl, digits=digits)
  if(nrow(x$weightMatrix)>max.print)
    cat(paste("... [ omitted ",nrow(x$weightMatrix)-max.print," rows ] ...\n",sep=""))
  cat("\n")

  invisible(x)
}
