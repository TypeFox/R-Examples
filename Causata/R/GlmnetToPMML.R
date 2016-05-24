
ToPmml <- function(model.definition, ...){
  UseMethod("ToPmml")
}

ToPmml.GlmnetModelDefinition <- function(model.definition, variable.definition, verbose=FALSE, ...) {
  model <- model.definition$model
  lambda <- model.definition$lambda
  causataData <- model.definition$causata.data
  formula <- model.definition$formula
  coefs <- coef(model, s=lambda)
  
  doc <- newXMLDoc()
  pmml <- newXMLNode(
    "PMML",
    namespaceDefinitions = c("http://www.dmg.org/PMML-4_0"),
    doc=doc)
  
  xmlAttrs(pmml)["version"] <- "4.0"
  InsertPmmlHeader.VariableDefinition(variable.definition, pmml)
  
  newXMLNode("DataDictionary", parent=pmml)
  
  regressionModel <- newXMLNode("RegressionModel", parent=pmml)
  insertRegressionModelAttributes(model, "model", regressionModel)
  
  newXMLNode("MiningSchema", parent=regressionModel)
  newXMLNode("LocalTransformations", parent=regressionModel)
  
  regressionTable <- newXMLNode("RegressionTable", parent=regressionModel, attrs=(c(intercept=0)))
  insertTargetCategory(model, regressionTable)
  
  non.zero.coefs <- coefs[coefs[,1] != 0, ]
  names(non.zero.coefs) <- str_extract(names(non.zero.coefs), "[^`]+")
  
  termsf <- terms.formula(formula)
  terms <- list()
  terms$variables <- attr(termsf, 'variables')
  terms$response  <- attr(termsf, 'response')
  terms$labels    <- attr(termsf, 'term.labels')
  terms$order     <- attr(termsf, 'order')
  terms$dv        <- terms$variables[[terms$response + 1]]
  # I'm adding 1 for the dv here, becuase terms.variables is a quoted expression list(foo, bar)
  # where [[1]] --> list, [[2]] --> foo, [[3]] --> bar
  
  crossTerms <- getCrossTerms(causataData, terms)
  
  if (verbose){
    cat("\n\nGenerating model PMML")
    cat("\n  lambda =", lambda)
    cat("\n  Found", length(non.zero.coefs), "nonzero coefficients for regression model.")
    cat("\n  CrossTerms:", crossTerms)
  }
  
  # loop over all of the nonzero coefficients
  variables.already.written <- c()
  if (verbose) cat("\n  Writing variables:")
  for (name in names(non.zero.coefs)){
    if (name=="(Intercept)") {
      xmlAttrs(regressionTable)["intercept"] = non.zero.coefs[[name]];
      if (verbose) {cat("\n    Writing intercept")}
    } else {
      # get the variable object for this coefficient
      variable <- GetVariable(causataData, r.name=name)
      
      if (name %in% names(crossTerms)) {
        # this is a cross term
        predictorTerm <- newXMLNode("PredictorTerm", parent=regressionTable, attrs=(c(coefficient=non.zero.coefs[[name]])))
        left.term <- crossTerms[[name]][[1]]
        right.term <- crossTerms[[name]][[2]]
        
        left.name <- RVariableToCausata(GetVariable(causataData, left.term)$causata.name)
        right.name <- RVariableToCausata(GetVariable(causataData, right.term)$causata.name) # TODO switch to use crossTerms[[name]]
        
        newXMLNode("FieldRef", parent=predictorTerm, attrs=(c(field=left.name)))
        newXMLNode("FieldRef", parent=predictorTerm, attrs=(c(field=right.name)))
        if (verbose) cat("\n    Writing cross term", left.name, right.name)
      } else {
        # not a cross term
        numericPredictor <- newXMLNode("NumericPredictor", parent=regressionTable)
        xmlAttrs(numericPredictor)["name"] = RVariableToCausata(name)
        xmlAttrs(numericPredictor)["coefficient"] = non.zero.coefs[[name]];
        if (verbose) cat("\n    Writing linear term", RVariableToCausata(name))
      }

      # We are looping over coefficients
      # One factor variable can have multiple coefficients, so if this factor has already been
      # written to the data dictionary etc. then skip the writing steps below
      if (variable$causata.name %in% variables.already.written){
        next
      }

      # This variable has not been written to data dictionary yet, write it
      if (variable$class == "factor") {
        insertDataDictionaryFieldForNominal(variable, pmml)
        insertMiningFieldForNominal(variable, "BLANK", pmml)
        insertFactorTransformation(variable, levels(causataData$df[[variable$name]]), "BLANK", pmml)
      } else {
        insertDataDictionaryFieldForContinuous(variable, pmml)
        insertMiningFieldForContinuous(variable, pmml)
        if ('Discretize' %in% variable$steps) {
          insertDiscretizedTransformation(variable, pmml, causataData, "Discretize")
        }
      }

      # now that the variable is written, update the list of written variables
      variables.already.written <- c(variables.already.written, variable$causata.name)
    }
  }
  return(pmml)
}  
  

# creates a list of all the legal crossed-variable names mapped to the variables which are crossed
# (variable names will be in the pmml variables, e.g. variableNameValue)
getCrossTerms <- function(causataData, terms) {
  
  variables <- causataData$variableList
  crossedVariables <- terms$labels[terms$order == 2] # e.g. "var1:var2", "foo:bar"
  
  if (length(terms$order[terms$order > 2]) > 0) {
    stop(paste("Cannot (yet) handle cross terms of order 3 or more, needed for:", terms$labels[terms$order > 2]))
  }
  
  # for each of the crossed variables, take thier components and produce all the legal names for them.
  # this is only complex for factors, for example two factors being crossed (var1:var2)
  # variable  var1   var2
  # levels    a:b    c:d  
  #           e:f    g:h
  # 
  # when crossed, these are the labels we'd expect to see:
  # var1a:b:var2c:d
  # var1a:b:var2g:h
  # var1e:f:var2c:d
  # var1e:f:var2g:h
  #
  
  crossedLabels <- list()
  
  crossTerm <- function(a, b) { paste(a, b, sep=":") }
  
  for (crossedVariable in crossedVariables) {
    left  <- strsplit(crossedVariable, ":")[[1]][[1]]
    right <- strsplit(crossedVariable, ":")[[1]][[2]]
    
    left.labels  <- getVariableLabels(causataData, left)
    right.labels <- getVariableLabels(causataData, right)
    
    for (left.label in left.labels) {
      for (right.label in right.labels) {
        crossedLabels[[paste(left.label, sep=":", right.label)]] <- c(left.label, right.label)
      }
    }
  }
  
  return(crossedLabels)
}

getVariableLabels <- function(causataData, variableName) {
  if (causataData$variableList[[variableName]]$class == "factor") {
    return( paste(variableName, sep="", levels(causataData$df[[variableName]])) )
  } else {
    return( variableName )
  }
}

insertDataDictionaryFieldForNominal <- function(variable, pmml) {
  dataDictionary <- xmlChildren(pmml)["DataDictionary"][[1]]
  dataField <- newXMLNode("DataField", parent=dataDictionary)
  xmlAttrs(dataField)["name"] <- variable$causata.name
  xmlAttrs(dataField)["optype"] <- "categorical"
  xmlAttrs(dataField)["dataType"] <- "string"
  dataField
}

insertDataDictionaryFieldForContinuous <- function(variable, pmml) {
  dataDictionary <- xmlChildren(pmml)["DataDictionary"][[1]]
  
  dataField <- newXMLNode("DataField", parent=dataDictionary)
  xmlAttrs(dataField)["name"] <- variable$causata.name
  xmlAttrs(dataField)["optype"] <- "continuous"
  xmlAttrs(dataField)["dataType"] <- "double"
  dataField
}

insertMiningFieldForContinuous <- function(variable, pmml) {
  missingValue <- variable$missingReplacement
  if (is.null(missingValue)) {
    missingValue <- 0
  }
  regressionModel <- xmlChildren(pmml)["RegressionModel"][[1]]
  miningSchema <- xmlChildren(regressionModel)["MiningSchema"][[1]]
  
  miningField <- newXMLNode("MiningField", parent=miningSchema)
  
  xmlAttrs(miningField)["name"] <- variable$causata.name
  xmlAttrs(miningField)["optype"] <- "continuous"
  xmlAttrs(miningField)["missingValueReplacement"] <- missingValue
  if (!is.null(variable$outlierLowerLimit) || !is.null(variable$outlierUpperLimit)) {
    xmlAttrs(miningField)["outliers"] <- "asExtremeValues"
    if (!is.null(variable$outlierLowerLimit)) {
      xmlAttrs(miningField)["lowValue"] <- c(variable$outlierLowerLimit)
    }
    if (!is.null(variable$outlierUpperLimit)) {
      xmlAttrs(miningField)["highValue"] <- c(variable$outlierUpperLimit)
    }
  }
}

# This inserts a MiningField to tell us how to handle missing values.
insertMiningFieldForNominal <- function(variable, missingValue, pmml) {
  regressionModel <- xmlChildren(pmml)["RegressionModel"][[1]]
  miningSchema <- xmlChildren(regressionModel)["MiningSchema"][[1]]
  
  miningField <- newXMLNode("MiningField", parent=miningSchema)

  xmlAttrs(miningField)["name"] <- variable$causata.name
  xmlAttrs(miningField)["optype"] <- "categorical"
  xmlAttrs(miningField)["missingValueReplacement"] <- missingValue
}

# This inserts NormDiscrete elements to create N variables from the input variable.
insertFactorTransformation <- function(variable, levels, missingValue, pmml) {
  
  # These go into the data transformations for the model
  regressionModel <- xmlChildren(pmml)["RegressionModel"][[1]]
  localTransformations <- xmlChildren(regressionModel)["LocalTransformations"][[1]]
  
  insert <- function(levelName) {
    derivedField <- newXMLNode("DerivedField", parent=localTransformations)
    xmlAttrs(derivedField)["name"] = paste(variable$causata.name, levelName, sep="")
    xmlAttrs(derivedField)["optype"] = "continuous"
    xmlAttrs(derivedField)["dataType"] = "double"
    
    normDiscrete <- newXMLNode("NormDiscrete", parent=derivedField)
    xmlAttrs(normDiscrete)["field"] <- variable$causata.name
    xmlAttrs(normDiscrete)["method"] <- "indicator"
    xmlAttrs(normDiscrete)["value"] <- levelName
    if (levelName == missingValue) {
      xmlAttrs(normDiscrete)["mapMissingTo"] <- 1
    } else {
      xmlAttrs(normDiscrete)["mapMissingTo"] <- 0
    }
  }
  
  lapply(levels, insert)
}

insertDiscretizedTransformation <- function(variable, pmml, causataData, stepname) {
  # These go into the data transformations for the model
  regressionModel <- xmlChildren(pmml)["RegressionModel"][[1]]
  localTransformations <- xmlChildren(regressionModel)["LocalTransformations"][[1]]
  
  derivedField <- newXMLNode("DerivedField", parent=localTransformations)
  #xmlAttrs(derivedField)["name"] = getBinnedVariableName(variable$name, causataData)
  xmlAttrs(derivedField)["name"] = variable$causata.name
  xmlAttrs(derivedField)["optype"] = "continuous"
  xmlAttrs(derivedField)["dataType"] = "double"
  
  discretize <- newXMLNode("Discretize", parent=derivedField, attrs=list(field=variable$causata.name, dataType="double"))
  binLimits <- variable$binLimits
  
  # set bin values according to the step name
  if (stepname == "Discretize") {
    # Bins are mapped to discrete values
    binValues <- variable$binValues
  } else {
    stop("Invalid step specified in GlmnetToPMML:", stepname)
  }
  # create an index of bins
  binIndex <- seq(1, length(binLimits) - 1)
  
  insert <- function(i.bin) {
    discretizeBin <- newXMLNode("DiscretizeBin", parent=discretize, attrs=list(binValue=binValues[i.bin]))
    interval <- newXMLNode("Interval", parent=discretizeBin)
    xmlAttrs(interval)["closure"] = "openClosed"
    if (i.bin > 1) {
      xmlAttrs(interval)["leftMargin"] = binLimits[i.bin]
    }
    if (i.bin < length(binValues)) {
      xmlAttrs(interval)["rightMargin"] = binLimits[i.bin + 1]
    }
  }
  lapply(binIndex, insert)
}


getBinnedVariableName <- function(name, causataData) {
  if ('BinContinuousP' %in% causataData$variableList[[name]]$steps) {
    return(paste(name, "binned", sep=""))
  }
  return(name)
}

insertRegressionModelAttributes <- function(model, modelName, regressionModel) {
  if (model$name == "Binomial Deviance") {
    addRegressionModelAttributes(regressionModel, modelName, "logisticRegression", "softmax")
  } else {
    if (model$name == "Mean-Squared Error") {
      addRegressionModelAttributes(regressionModel, modelName, "linearRegression", "none")
    } else {
      stop("The model must be either linear regression or logistic regression (built with family=\"binomial\")")
    }
  }
}

addRegressionModelAttributes <- function(regressionModel, modelName, modelType, normalizationMethod) {
  xmlAttrs(regressionModel)["modelName"] = modelName
  xmlAttrs(regressionModel)["functionName"] = "regression"
  xmlAttrs(regressionModel)["modelType"] = modelType
  xmlAttrs(regressionModel)["normalizationMethod"] = normalizationMethod
}

insertTargetCategory <- function(model, regressionTable) {
  if (model$name == "Binomial Deviance") {
    xmlAttrs(regressionTable)["targetCategory"] = "YES"
  }
}

InsertPmmlHeader.VariableDefinition <- function(this, pmml) {
  header <- newXMLNode("Header", attrs=c(copyright="Causata", description=""), parent = pmml)
  insertHeaderExtension("VariableName", this$name, header)
  insertHeaderExtension("DisplayName", this$display.name, header)
  insertHeaderExtension("Description", this$description, header)
  for (label in this$labels) {
    insertHeaderExtension("Label", label, header)
  }
  timestampNode <- newXMLNode("Timestamp", parent=header)
  xmlValue(timestampNode) <- this$timestamp
  return(pmml)
}

insertHeaderExtension <- function(name, value, header) {
  if (!is.null(value)) {
    extension <- newXMLNode("Extension", parent = header);
    xmlAttrs(extension)["name"] = name
    xmlAttrs(extension)["value"] = value
  }
}

RVariableToCausata <- function(vname){
  # Replaces the first double underscore in a string __ with a dollar sign $
  # Used to match variable names in coefficient table with names in other parts of XML
  return(sub('__', '\\$', vname))
}
