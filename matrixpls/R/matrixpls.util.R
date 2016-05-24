# =========== Utility functions ===========

estimatesMatrixToVector <- function(est, m, sep, reverse = FALSE){
  
  ind <- which(m != 0)
  ret <- est[ind]
  
  if(reverse){
    names(ret) <- paste(colnames(m)[col(m)[ind]],
                        rownames(m)[row(m)[ind]], sep=sep)
  }
  else{
    names(ret) <- paste(rownames(m)[row(m)[ind]],
                        colnames(m)[col(m)[ind]], sep=sep)
  }
  
  ret
}


#
# Scales the weight matrix so that the resulting composites have unit variance
#

scaleWeights <- function(S, W){
  
  # Calculate the variances of the unscaled composites
  
  var_C_unscaled <- diag(W %*% S %*% t(W))
  
  # Scaling the unscaled weights and return
  
  ret <- diag(x = 1/sqrt(var_C_unscaled)) %*% W
  
  colnames(ret) <- colnames(W)
  rownames(ret) <- rownames(W)
  
  return(ret)
}

lavaanParTableToNativeFormat <- function(partable){
  
  factorLoadings <- partable[partable$op == "=~",]
  regressions <- partable[partable$op == "~",]
  formativeLoadings <- partable[partable$op == "<~",]
  
  # Parse the variables
  latentVariableNames <- unique(c(factorLoadings$lhs, formativeLoadings$lhs))
  observedVariableNames <- setdiff(unique(c(factorLoadings$rhs,formativeLoadings$rhs,
                                            regressions$lhs,regressions$rhs)), latentVariableNames)
  
  # Set up empty model tables
  
  inner <- matrix(0,length(latentVariableNames),length(latentVariableNames))
  colnames(inner)<-rownames(inner)<-latentVariableNames
  
  formative <- matrix(0,length(latentVariableNames),length(observedVariableNames))
  colnames(formative) <- observedVariableNames
  rownames(formative) <- latentVariableNames
  
  reflective <- t(formative)
  
  # Set the relationships in the tables
  
  rows <- match(factorLoadings$rhs, observedVariableNames)
  cols <- match(factorLoadings$lhs,latentVariableNames)
  indices <- rows + (cols-1)*nrow(reflective) 
  
  reflective[indices] <- 1
  
  latentRegressions <- regressions[regressions$rhs %in% latentVariableNames & 
                                     regressions$lhs %in% latentVariableNames,]
  
  rows <- match(regressions$lhs, latentVariableNames)
  cols <- match(regressions$rhs,latentVariableNames)
  indices <- rows + (cols-1)*nrow(inner) 
  
  inner[indices] <- 1
  
  formativeRegressions <- rbind(regressions[regressions$rhs %in% observedVariableNames & 
                                              regressions$lhs %in% latentVariableNames,],
                                formativeLoadings)
  
  
  rows <- match(formativeRegressions$lhs, latentVariableNames)
  cols <- match(formativeRegressions$rhs, observedVariableNames)
  indices <- rows + (cols-1)*nrow(formative) 
  
  formative[indices] <- 1
  return(list(inner=inner, reflective=reflective, formative=formative))
}

parseModelToNativeFormat <- function(model){
  
  if(is.matrixpls.model(model)){
    #Already in native format
    return(model)
  }
  else if(is.character(model)) {
    # Remove all multigroup specifications because we do not support multigroup analyses
    model <- gsub("c\\(.+?\\)","NA",model)
    
    return(lavaanParTableToNativeFormat(lavaan::lavaanify(model)))
  } else if (is.partable(model)) {
    return(lavaanParTableToNativeFormat(model))
  } else if (is.lavaancall(model)) {
    browser()
  } else if (methods::is(model, "lavaan")) {
    return(lavaanParTableToNativeFormat(model@ParTable))
    #	} else if (methods::is(model, "MxModel")) {
    #		browser()
    #	} else if (methods::is(model, "SimSem")) {
    #		browser()
  } else {
    stop("Please specify an appropriate object for the 'model' argument: simsem model template, lavaan script, lavaan parameter table, list of options for the 'lavaan' function, or a list containing three matrices in the matrixpls native model format.")
  }
}

#
# Returns a matrix with composites on rows and observed on columns
#

defaultWeightModelWithModel <- function(model){
  
  nativeModel <- parseModelToNativeFormat(model)
  W.model <- matrix(0,nrow(nativeModel$formative), ncol(nativeModel$formative))
  
  colnames(W.model) <- colnames(nativeModel$formative)
  rownames(W.model) <- rownames(nativeModel$formative)
  
  W.model[nativeModel$formative!=0] <- 1 
  W.model[t(nativeModel$reflective)!=0] <- 1 
  
  return(W.model)
  
}

is.matrixpls.model <- function(model) {
  return (is.list(model) &&
            setequal(names(model),c("inner","reflective","formative")) &&
            all(sapply(model,is.matrix)))
  
}

#
# Functions adapted from simsem
#

is.partable <- function(object) {
  is.list(object) && all(names(object) %in% c("id", "lhs", "op", "rhs", "user", "group", "free", "ustart", "exo", "label", "eq.id", "unco"))
}

is.lavaancall <- function(object) {
  is.list(object) && ("model" %in% names(object))
}
