#' Creates a simnple linear model
#'
#' Returns a simple linear regression model, a linear regression model with a single
#' explanatory variable
#'
#' The function eliminates all pairs for which either the first field or the second field
#' is empty. After the elimination, if the length of the input is less than 2, the function
#' returns the empty sequence. After the elimination, if the standard deviation of the
#' independent variable is 0, the function returns a linear model with intercept = the mean
#' of the dependent variable, coefficients = NaN and r-squared = NaN.
#' After the elimination, if the standard deviation of the dependent variable is 0,
#' the function returns a linear model with r-squared = NaN.
#'
#' @param form an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param mlDf an ml.data.frame object

#' @export
ml.lm <- function(form, mlDf) {

  conn <- mlDf@.conn
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  queryComArgs <- mlDf@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.lm", sep="")
  nPageLength <- mlDf@.nrows
  queryArgs <- c(queryComArgs, 'rs:pageLength'=nPageLength)

  #

  isResponse <- attr(terms(form, keep.order=T, data=data.frame(x=1)), "response")
  vars <- all.vars(form)
  if (length(vars) != 2) {
    stop("Can only use two variables!")
  }
  # below will return a tabel with all fields on the right to the  ~ as columnnames and
  # all variables as rownames
  tab1 <- attr(terms(form), "factors")
  # our independent variable is in the columnname
  independent <- colnames(tab1)
  # dependent
  dependent <- vars[isResponse]
  if (!(independent %in% mlDf@.col.name) && !(dependent %in% mlDf@.col.name)) {
    stop("Both variables must be part of a ml.data.frame")
  }
  # need to very that they are number fields...

  fields <- "{"
  # check if dependent or independent is existing fields
  # or new, if new we ned to use the expersion
  if (is.null(mlDf@.col.defs[[dependent]])) {
    fieldDefDep <- paste('rfmlResult[\'', dependent,'\']',sep='')
    i <- match(dependent, mlDf@.col.name)
    fieldOrgNameDep <- mlDf@.col.org_name[i]
    fieldFormatDep <- mlDf@.col.format[i]
  } else {
    fieldDefDep <- mlDf@.col.defs[[dependent]]
    fieldOrgNameDep <- dependent
    fieldFormatDep <- ''
  }
  if (is.null(mlDf@.col.defs[[independent]])) {
    fieldDefInd <- paste('rfmlResult[\'', independent,'\']',sep='')
    i <- match(independent, mlDf@.col.name)
    fieldOrgNameInd <- mlDf@.col.org_name[i]
    fieldFormatInd <- mlDf@.col.format[i]
  } else {
    fieldDefInd <- mlDf@.col.defs[[independent]]
    fieldOrgNameInd <- independent
    fieldFormatInd <- ''
  }

  fields <- paste(fields, '"',dependent , '":{"fieldDef":"',fieldDefDep,'","orgField":"', fieldOrgNameDep, '","orgFormat":"', fieldFormatDep  ,'"},"',
                  independent, '":{"fieldDef":"',fieldDefInd,'","orgField":"', fieldOrgNameInd, '","orgFormat":"', fieldFormatInd ,'"}' ,sep='')

  fields <- paste(fields, '}', sep='')
  #message(fields)
  queryArgs <- c(queryArgs, 'rs:fields'=fields)

  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())

  rContent <- content(response) #, as = "text""
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  res <- list(intercept=rContent$`math:linear-model`$intercept, coefficients=rContent$`math:linear-model`$coefficients, rsquared=rContent$`math:linear-model`$rsquared)
  class(res) = c("mlLm")
  res
}
#' Prints information for a simnple linear model returned by \link{ml.lm}
#'
#' @param x a ml.lm result
#' @param ... not used
#' @method print mlLm
#' @export
print.mlLm <- function(x, ...) {
  cat("intercept: ", x$intercept)
  cat("\ncoefficients: ", x$coefficients)
  cat("\nr-squared :", x$rsquared, "\n")
}
