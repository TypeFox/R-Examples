#' Fit Linear Model and return its coefficients.
#' 
#' Outputs coefficients of the linear model fitted to Aster table according
#' to the formula expression containing column names. The zeroth coefficient corresponds 
#' to the slope intercept. R formula expression with column names for response and 
#' predictor variables is exactly as in \code{\link{lm}} function (though less 
#' features supported).
#' 
#' Models for \code{computeLm} are specified symbolically. A typical model has the form 
#' \code{response ~ terms} where response is the (numeric) column and terms is a series of
#' column terms which specifies a linear predictor for response. A terms specification of 
#' the form \code{first + second} indicates all the terms in first together with all the 
#' terms in second with duplicates removed. A specification of the form \code{first:second} 
#' and \code{first*second} (interactions) are not supported yet.
#' 
#' @param channel connection object as returned by \code{\link{odbcConnect}}
#' @param tableName Aster table name
#' @param formula an object of class "formula" (or one that can be coerced to that class): 
#'   a symbolic description of the model to be fitted. The details of model 
#'   specification are given under `Details`.
#' @param tableInfo pre-built table summary with data types 
#' @param categories vector with column names containing categorical data. Optional if the column is of
#'   character type as it is automatically treated as categorical predictors. But if numerical 
#'   column contains categorical data then then it has to be specified for a model to view it
#'   as categorical. Apply extra care not to have columns with too many values (approximaltely > 10) 
#'   as categorical because each value results in dummy predictor variable added to the model.
#' @param sampleSize function always computes regression model coefficent on all data in the table.
#'   But it computes predictions and returns an object of \code{\link{class}} "lm" based on sample
#'   of data. The sample size is in an absolute value for number of rows in the sample. 
#'   Be careful not overestimating the size as all results are loaded into memory. 
#'   Special value \code{"all"} or \code{"ALL"} will include all data in computation.
#' @param where specifies criteria to satisfy by the table rows before applying
#'   computation. The creteria are expressed in the form of SQL predicates (inside
#'   \code{WHERE} clause).
#' @param test logical: if TRUE show what would be done, only (similar to parameter \code{test} in \link{RODBC} 
#'   functions like \link{sqlQuery} and \link{sqlSave}).  
#' @return
#' \code{computeLm} returns an object of \code{\link{class}} \code{"toalm", "lm"}.
#' 
#' The function \code{summary} .....
#' 
#' For backward compatibility 
#' Outputs data frame containing 3 columns:
#' \describe{
#'   \item{coefficient_name}{name of predictor table column, zeroth coefficient name is "0"}
#'   \item{coefficient_index}{index of predictor table column starting with 0}
#'   \item{value}{coefficient value}
#' }
#' @export
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' # batting average explained by rbi, bb, so 
#' lm1 = computeLm(channel=conn, tableName="batting_enh", formula= ba ~ rbi + bb + so)
#' summary(lm1)
#' 
#' # with category predictor league and explicit sample size
#' lm2 = computeLm(channel=conn, tableName="batting_enh", formula= ba ~ rbi + bb + so + lgid,
#'                 , sampleSize=10000, where="lgid in ('AL','NL') and ab > 30") 
#' summary(lm2)
#' }
#' 
computeLm <- function(channel, tableName, formula, tableInfo = NULL, categories = NULL,
                      sampleSize = 1000,
                      where = NULL, test = FALSE) {
  
  if (missing(channel)) {
    stop("Must provide connection.")
  }
  
  if (missing(tableName) || missing(formula)) {
    stop("Must provide table and expression.")
  }
  
  if (missing(tableInfo) && test) {
    stop("Must provide tableInfo when test==TRUE.")
  }
    
  cl <- match.call()
  ft = stats::terms(formula)
  #fvars = all.vars(formula)
  vars = as.character(attr(ft, "variables"))[-1]
  if (length(vars) < 2) {
    stop("No predictors found in formula.")
  }
  
  responseIdx = attr(ft, "response")
  if (responseIdx < 1) {
    stop("No response variable found in formula.")
  }
  
  responseVar = vars[[responseIdx]]
  
  predictors = vars[-responseIdx]
  
  
  # check columns and their data types to make dummy variables for categorical predictors
  if (missing(tableInfo)) {
    summary = sqlColumns(channel, tableName)
  }else {
    summary = includeExcludeColumns(tableInfo, include=vars, except=NULL)
  }
  
  num_cols = getNumericColumns(summary, names.only=TRUE)
  char_cols = getCharacterColumns(summary, names.only=TRUE)
  if (!all(vars %in% c(num_cols, char_cols))) {
    stop(paste0("Columns ", paste(vars[!vars %in% c(num_cols, char_cols)], collapse=", "), 
                " are not found in table ", tableName, ". Make sure all predictors are character or numeric types."))
  }
  
  if (!responseVar %in% num_cols) {
    stop("Response variable's column is not of numeric type.")
  }
  
  char_predictors = predictors[predictors %in% char_cols]
  num_predictors = predictors[predictors %in% num_cols]
  
  # categories defines numerical columns that are actually categorical
  if (!is.null(categories) && length(categories) > 0) {
    nums_to_move = intersect(num_predictors, categories)
    num_predictors = setdiff(num_predictors, nums_to_move)
    char_predictors = union(char_predictors, nums_to_move)
  }
  
  predictorColumns = num_predictors
  predictorNames = num_predictors
  predictorNamesSQL = predictorNames
  
  # implement categorical variables
  xlevels = NULL
  if (length(char_predictors) > 0) {
    
    xlevels = list()
    
    for(name in char_predictors) {
      col_values = getColumnValues(channel, tableName, name, where, test)
      if (length(col_values) < 2) {
        stop(paste0("Categorical column '", name, "' has 1 or no values. Consider removing or replacing it."))
      }
      numValuesBefore = length(unique(col_values))
      col_values = gsub("[^[:alnum:]_]", "", col_values)
      if (length(unique(col_values)) < numValuesBefore) {
        stop(paste0("Categorical column '", name, "' values are not valid strings. Consider replacing them with alpha-numeric versions."))
      }
      
      col_values = col_values[order(col_values)]
      xlevels[[name]] = col_values
      col_values = col_values[-1]
      
      # case when lgid = 'NL' then 1 else 0 end as "lgid_NL"
      for (value in col_values) {
        categoryExpr = paste0("CASE WHEN ", name, " = '", value, "' THEN 1 ELSE 0 END")
        categoryName = paste0(tolower(name), value)
        categoryNameSQL = paste0("\"", tolower(name), value, "\"")
        
        predictorColumns = c(predictorColumns, categoryExpr)
        predictorNames = c(predictorNames, categoryName)
        predictorNamesSQL = c(predictorNamesSQL, categoryNameSQL)
      }
    }
  }
  
  predictorList = paste0(predictorColumns, rep(" x", length(predictorColumns)), 
                         as.character(seq(1, length(predictorColumns))), collapse=", ")
  selectList = paste0(predictorList, ", ", responseVar, " y ")

  where_clause = makeWhereClause(where)
  
  sql = paste0(
        "SELECT * 
           FROM linreg(
                  ON linregmatrix(
                    ON (SELECT ", selectList, " FROM ", tableName, where_clause, ")
                  )
                  PARTITION BY 1
           )"
  )
  
  if (test) 
    return (sql)
  else {
    result = toaSqlQuery(channel, sql, stringsAsFactors=FALSE)
  }
  
  # handle empty data set: currently Aster ODBC driver doesn't return SQL/MR error:
  # ERROR: SQL-MR function LINREG failed: The input data results in a singular matrix and hence there is no solution
  # so we handle it when no results returned
  if (nrow(result) == 0) 
    stop(paste0("Please re-run computeLm with test=TRUE and check that generated sql indeed triggers following error:
ERROR: SQL-MR function LINREG failed: The input data results in a singular matrix and hence there is no solution
Then refer to Aster Analytics Foundation Guide, Linear Regression Function, in particular:
If two or more input columns are co-linear, or very closely correlated, then no solution to
linear regression exists, so the function will fail. Looking at correlations between columns
using Aster Database correlation (stats correlation) function can help uncover sources of colinearity.
Removing co-linear columns should resolve the issue.
This inconvinience will be addressed in one of future releases of toaster."))
  
  z = createLm(channel, tableName, stats::as.formula(formula), cl, result$value, ft, xlevels, 
               predictors, predictorColumns, predictorNames, predictorNamesSQL, 
               responseVar, sampleSize, where)
  
  # for previous version 0.2.5 support (not backward compatible)
  z$old.result = cbind(data.frame(coefficient_name = c("0", predictorNames)), result)
  
  return(z)
}

createLm <- function(channel, tableName, formula, cl, coefficients, ft, xlevels, 
                     predictors, predictorColumns, predictorNames, predictorNamesSQL, 
                     response, sampleSize, where) {
  
  z <- structure(list(coefficients = coefficients,
                      rank = length(predictors) + 1,
                      call = cl,
                      terms = ft,
                      contrasts = NULL,
                      xlevels = xlevels

                      ),
                 class = c("toalm", "lm"))
  names(z$coefficients) = c("(Intercept)", predictorNames)
  
  if (!is.null(sampleSize) && sampleSize >= 30) {
    fit = predictLm(channel, tableName, predictors, predictorColumns, predictorNamesSQL, response,
                    coefficients, sampleSize, where)
    
    rownumbers = fit$`__row_number__`
    
    z$residuals = fit[, response] - fit$`__yt__`
    names(z$residuals) = rownumbers
    
    z$fitted.values = fit$`__yt__`
    names(z$fitted.values) = rownumbers
    
    z$df.residual = length(rownumbers) - z$rank
    
    z$qr = qr(cbind(1, fit[, predictorNames]))
    
    # for now model.frame is compatible with models without categorical predictors only
    if (is.null(xlevels) || length(xlevels)==0)
      z$model = stats::model.frame(formula=formula, data=fit)
  }else
    warning("No sampling performed if sample size is NULL or < 30.")
  
  return (z)
}

predictLm <- function(channel, tableName, predictors, predictorColumns, predictorNamesSQL, response,
                      coefficients, sampleSize, where) {
  
  
  # omit rows with NULLs for predictors or response (na.omit action)
  whereNotNull = paste0(predictors, " IS NOT NULL", collapse=" AND ")
  whereNotNull = paste0(whereNotNull, " AND ", response, " IS NOT NULL ")
  if (is.null(where)) 
    where = whereNotNull
  else 
    where = paste0(where, " AND ( ", whereNotNull, " ) ")
  
  # generate sample SQL
  if (tolower(sampleSize) == "all") {
    selectSample = paste0("SELECT * FROM ", tableName,  makeWhereClause(where))
  }else
    selectSample = computeSample(NULL, tableName=tableName, sampleSize=sampleSize, where=where, test=TRUE)
  
  a0 = coefficients[[1]]
  coefficients = coefficients[-1]
  
  predictExpr = paste0(" ( ", coefficients, " * ( ", predictorColumns, " ) ) ", collapse = " + ")
  predictExpr = paste0(" ( ", a0, " + ", predictExpr, " ) __yt__ ")
  selectList = paste0(predictorColumns, " ", predictorNamesSQL, collapse = ", ")
  selectList = paste0(selectList, ", ", predictExpr, ", ", response, ", row_number() over (order by 1) as __row_number__ ")
  
  sql = paste0(
    "SELECT ", selectList, " FROM ( ", selectSample, " ) t "
  )
  
  result = toaSqlQuery(channel, sql, stringsAsFactors=FALSE)
  
  # handle hidden errors by analyzing number of rows
  if (nrow(result) == 0) {
    stop("No rows to compute lm characteristics. One of the reason could be insufficient number of rows for sample - try using sampleSize='ALL'.")
  }
  
  # sometimes numeric columns get converted to character 
  # (guessing due to large number of digits after decimal point)
  for (i in 1:length(result)) {
    if (class(result[,i]) == "character") 
      result[,i] = as.numeric(result[,i])
  }
  
  return(result)
}

