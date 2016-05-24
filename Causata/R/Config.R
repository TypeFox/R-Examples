# Functions for dealing with Causata configuration
#

library("RCurl")
library("rjson")
library("glmnet")
library("XML")

CausataConfig <- function(config.server.host, config.server.port, config.username, config.password, protocol="https://", group=NULL) {
  # create a list of arguments supplied in the function call, including defaults
  arglist.fun <- mget(names(formals()),sys.frame(sys.nframe()))
  # arguments with no values will have length zero, remove them
  arglist.fun <- arglist.fun[as.numeric(lapply(as.character(arglist.fun), FUN=str_length)) > 0]
  # create another argument list extracted from causata config file
  arglist.yaml <- LoadCausataConfig(group)
  # combine the arguments from these two sources into a single list
  # if an argument exists in both lists then the function arguments take precendence since they are first
  arglist <- c(arglist.fun, arglist.yaml)

  this <- list(host=arglist$config.server.host, port=arglist$config.server.port, user=arglist$config.username, 
               pass=arglist$config.password, protocol=arglist$protocol)
  class(this) <- "CausataConfig"
  this
}
is.CausataConfig <- function(this) inherits(this, "CausataConfig")

GetModelUploadUrl <- function(this) {
  paste(this$protocol, this$host, ":", this$port, "/status/model", sep="", collapse="")
}

##########################################################################################
# VariableDefinition
##########################################################################################
VariableDefinition <- function(
  name,
  display.name=name,
  description=name,
  labels=list(),
  author=Sys.info()[["user"]],
  timestamp=as.integer(1000*as.numeric(format(Sys.time(),"%H%M%OS3"))),
  archived=FALSE,
  categorizing.attribute="",
  output.value.count=-1,
  data.type="double") {

  if (nchar(categorizing.attribute) > 0 && output.value.count < 1) {
    warning("Cannot create a categorized variable with non-positive output value count")
    return(NULL)
  }
  
  this <- list(
    name=name,
    display.name=display.name,
    description=description,
    labels=labels,
    author=author,
    timestamp=timestamp,
    archived=archived,
    categorizing.attribute=categorizing.attribute,
    output.value.count=output.value.count,
    data.type=data.type)
  class(this) <- "VariableDefinition"
  this
}
is.VariableDefinition <- function(this) inherits(this, "VariableDefinition")

is.numeric.VariableDefinition <- function(this) {
  this$data.type == "double"  ||
  this$data.type == "float"   ||
  this$data.type == "long"    ||
  this$data.type == "integer" ||
  this$data.type == "short"
}

##########################################################################################
# ModelDefinition
##########################################################################################
ModelDefinition <- function(model, ...){
  UseMethod("ModelDefinition")
}

ModelDefinition.cv.glmnet <- function(
  model,
  causata.data,
  formula,
  lambda=model$lambda.1se, ...) {

  this <- list(model=model, causata.data=causata.data, formula=formula, lambda=lambda)
  class(this) <- c("GlmnetModelDefinition", "ModelDefinition")
  return(this)
}
is.ModelDefinition <- function(this) inherits(this, "ModelDefinition")

ToPmml <- function(model.definition, ...){
  UseMethod("ToPmml")
}

GetQuery <- function(this, ...){
  UseMethod("GetQuery")
}

GetQuery.ModelDefinition <- function(this, ...) {
  GetQuery(this$causata.data)
}

predict.GlmnetModelDefinition <- function(object, data, verbose=FALSE, ...) {
  if (verbose) cat("\nGenerating glmnet prediction, lambda =", object$lambda, "\n")
  missing.cols <- c() # default value for missing columns
  
  # remove columns from the data frame that are not used in the formula
  if (verbose) cat('Columns in input data:', ncol(data), '\n')
  data <- data[, names(data) %in% all.vars(object$formula)]
  if (verbose) cat('Columns in input data after removing unused:', ncol(data), '\n')
  
  # set contrasts in data to use all levels in each provided factor
  # don't use the default behavior that omits the first factor level
  mm <- model.matrix(object$formula, data=data, 
    contrasts.arg=lapply(data[, sapply(data, is.factor)], contrasts, contrasts=FALSE) )
  
  if (nrow(mm) != nrow(data)){
    warning("Rows in matrix do not match rows in data, probably caused by missing values in data.")
  }
  
  # ensure that the columns in the new model matrix are aligned with the coefficient matrix
  match.idx <- match(rownames(object$model$glmnet.fit$beta), colnames(mm))
  
  if (any(is.na(match.idx))) {
    # missing columns in model matrix, display warning and add columns of zeros
    missing.cols <- rownames(object$model$glmnet.fit$beta)[is.na(match.idx)]
    msg <- paste("Missing columns (", sum(is.na(match.idx)),
            ") in model matrix for predict.GlmnetModelDefinition, substituting columns of zeros", sep="")
    warning(msg)
    if (verbose) cat(paste(msg, missing.cols, collapse="\n  "))
    # loop for each missing column, append column of zeros to model matrix
    zeros <- rep(0, nrow(mm))
    for (missing.col in missing.cols){
      mm <- cbind(mm, zeros)
      colnames(mm)[ncol(mm)] <- missing.col
    }
  }
  
  # ensure that column order in new matrix matches order in model coefficients
  match.idx <- match(rownames(object$model$glmnet.fit$beta), colnames(mm))
  mm <- mm[, match.idx]
  
  # we added missing columns above, missing columns not allowed now
  stopifnot(!any(is.na(match.idx)))
  # require that columns in new matrix exactly match coefficient names and order
  stopifnot(all(rownames(object$model$glmnet.fit$beta) == colnames(mm)))
  
  predicted <- predict(object$model, newx=mm, type="response", s=object$lambda)
  return(list(
    model.matrix = mm,
    predicted = as.numeric(predicted),
    lambda = object$lambda,
    missing.cols = missing.cols))
}


#
# ToPmml is defined for GlmnetModelDefinition
#


##########################################################################################
# Misc
##########################################################################################

uploadFile <- function(url, user, pass, payload, verbose=FALSE) {
  # Uploads the text in payload to url in a multipart upload.
  #
  content <- paste(
    '--BOUNDARY',
    'Content-Disposition: form-data; name="static-data"; filename="fake.json"',
    '',
    payload,
    '--BOUNDARY--',
    '',
    sep="\r\n", collapse="\r\n")
  
  opts <- curlOptions(
    postfields=content,
    verbose=verbose,
    ssl.verifypeer=FALSE,
    ssl.verifyhost=FALSE,
    followlocation=TRUE,
    userpwd=paste(user, pass, sep=":", collapse=":"),
    httpauth = 1L,
    httpheader=c('Content-Type'='multipart/form-data; boundary=BOUNDARY', Accept='application/json', 'Expect'='')
  )
  
  return(postForm(url, .opts=opts))
}

##########################################################################
# Variable creation and deletion 
##########################################################################

Config.CreatePrimaryVariable <- function(
  causata.config,
  variable.name,
  variable.display.name=variable.name,
  variable.expression) {
  
  # Creates a primay causata variable.
  # A primary variable is one that operates only on events and attributes, rather than on other variables
  #
  # Params:
  # causata.config: A CausataConfig object
  # variable.name: the (system) name of the variable to create
  # variable.display.name: the display name of the variable (defaults to the name)
  # variable.expression: the expression to use for this variable e.g. "COUNT purchase"
  #
  # Returns: 
  #   TRUE if the variable is successfully created or updated
  #   character data otherwise - the response from the configuration server that should provide an error message
  
  json.object <- list(
    "dataSetName" = "current-set",
    "modifications" = list(
      "variables" = list(
        list(
          "name" = variable.name,
          "displayName" = variable.display.name,
          "expression" = variable.expression,
          "type" = "PRIMARY"
        )
      )
    )
  )
  
  url <- paste(causata.config$protocol, causata.config$host, ":", causata.config$port, "/config/static-data", sep="", collapse="")
  response <- uploadFile(url, causata.config$user, causata.config$pass, toJSON(json.object))
  
  if (grepl("Inserted successfully", response, fixed=TRUE) || grepl("Updated successfully", response, fixed=TRUE)) {
    return(TRUE)
  } else {
    return(response)
  }
}

Config.DeleteVariable <- function(causata.config, variable.name) {
  
  # Deletes a causata variable.
  # Params:
  # causata.config: A CausataConfig object
  # variable.name: the (system) name of the variable to delete - same as the name with which it was created.
  #
  # Returns: 
  #   TRUE if the variable is successfully deleted
  #   character data otherwise - the response from the configuration server that should provide an error message
  
  json.object <- list(
    "dataSetName" = "current-set",
    "deletions" = list(
      "variableNames" = list(variable.name)
    )
  )
  
  url <- paste(causata.config$protocol, causata.config$host, ":", causata.config$port, "/config/static-data", sep="", collapse="")
  response <- uploadFile(url, causata.config$user, causata.config$pass, toJSON(json.object))
  
  if (grepl("Deleted successfully", response, fixed=TRUE)) {
    return(TRUE)
  } else {
    return(response)
  }
}


UploadModel <- function(causata.config, model.definition, variable.definition, verbose=FALSE) {
  pmml <- ToPmml(model.definition, variable.definition, verbose)
  url <- GetModelUploadUrl(causata.config)
  if (verbose) {
    cat(paste("Uploading to:", url))
  }
  response <- uploadFile(url, causata.config$user, causata.config$pass, saveXML(pmml), verbose=verbose)

  if (grepl("success", response, fixed=TRUE)) {
    return(TRUE)
  } else {
    return(response)
  }
}


ValidateModel <- function(causata.config, model.definition, test.variable.definition, connection,
          query.function,
          record.error.max, verbose=FALSE, ...) {
  # This funciton re-executes the query from the CausataData object in model.definition with the model
  # variable named in variable.definition included in the new query.
  # When some data has been downloaded (the number of records is defined in the record.count parameter)
  # predict is used to get the 'actual' value from the model in model.definition
  # the values from this are compared to the values from the model variable, and if any pair has a difference
  # of more than record.error.max, then the incorrect record is returned.
  # otherwise TRUE is returned.
  if (verbose) cat("\n\nValidating model")
  if (verbose) cat("\n  Running new query to get validation data")
  query <- query.function(test.variable.definition$name, ...)
  validation.data.df <- GetData(connection, query)
  
  if (verbose) cat("\n  Re-applying data transformations to new data frame")
  transformations <- GetTransforms(model.definition$causata.data)
  validation.data.df <- transformations( validation.data.df )
  
  if (verbose) cat("\n  Getting predictions from native model, lambda =", model.definition$lambda)
  predict.list <- predict(model.definition, validation.data.df, verbose=verbose)
  predictions <- predict.list$predicted
  
  if (verbose) cat("\n  Comparing with predictions from Causata model")
  dv.index <- grep(test.variable.definition$name, names(validation.data.df))
  stopifnot(length(dv.index)==1) # ensure that the dv exists in the data, should be 1 match
  actuals <- validation.data.df[,dv.index]
  
  errors <- abs(predictions - actuals)
  problematic.indices <- errors > record.error.max
  if (max(errors) <= record.error.max) {
    if (verbose) cat("\n  Validation passed")
    result <- TRUE
  } else {
    if (verbose) cat("\n  Found errors")
    result <- FALSE
  }
  return(list(
    result = result,
    validation.data = validation.data.df,
    errors = errors,
    predictions = predictions,
    model.matrix = predict.list$model.matrix,
    actuals = actuals,
    problematic.indices = problematic.indices
    ))
}


WaitForVariable <- function(connection, variable.definition, timeout.ms) {
  var.name <- variable.definition$name
  
  start <- as.integer(1000*as.numeric(format(Sys.time(),"%H%M%OS3")))
  while (!VariableExists(connection, var.name) && as.integer(1000*as.numeric(format(Sys.time(),"%H%M%OS3"))) < (start + timeout.ms)) {
    Sys.sleep(3.0)
  }
  VariableExists(connection, var.name)
}

VariableExists <- function(conn, variable.name) {
  return(variable.name %in% GetMetadata(conn)$variables$NAME)
}

DeleteModel <- function(causata.config, variable.definition) {
  Config.DeleteVariable(causata.config, variable.definition$name)
}

# Uploads the given model and ensures that the Causata model matches the R version of the model.
# If all goes well, TRUE is returned.  otherwise, the customer rows for which the validation failed
# are returned.
#
UploadModelWithValidation <- function(
  causata.config,
  model.definition,
  variable.definition,
  connection,
  query.function,
  record.error.max,
  verbose=FALSE, 
  ...) {
  
  test.variable.name <- paste(sample(letters[1:26],30,replace=TRUE),collapse="")
  
  test.variable.definition <- variable.definition
  test.variable.definition$name <- test.variable.name
  test.variable.definition$display.name <- test.variable.name
  test.variable.definition$description <- "Temporary model for validation"
  test.variable.definition$archived <- TRUE
  
  tryCatch({
    if (UploadModel(causata.config, model.definition, test.variable.definition, verbose) == TRUE) {
      if (WaitForVariable(connection, test.variable.definition, timeout.ms=30000)) {
        validated.list <- ValidateModel(
          causata.config, model.definition, test.variable.definition, connection,
          query.function,
          record.error.max,
          verbose, ...)
      } else {
        print(paste("Model for validation (named", test.variable.definition$name, ") wasn't found after 30 seconds"))
      }
    }
  }, finally={
    DeleteModel(causata.config, test.variable.definition)
  })
  
  if (validated.list$result) {
    # validation succeded, upload the model
    UploadModel(causata.config, model.definition, variable.definition, verbose)
  }

  # return validation data
  return(validated.list)
}
