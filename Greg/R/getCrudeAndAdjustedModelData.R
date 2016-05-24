#' This function helps with printing regression models
#' 
#' This function is used for getting the adjusted and unadjusted values
#' for a regression model. It takes a full model and walks through each
#' variable, removes in the regression all variables except one then
#' reruns that variable to get the unadjusted value. This functions not
#' intended for direct use, it's better to use \code{\link{printCrudeAndAdjustedModel}}
#' that utilizes this function.
#' 
#' This function saves a lot of time creating tables since it compiles a fully
#' unadjusted list of all your used covariates.
#' 
#' If the model is an exponential poisson/logit/cox regression model then it automatically
#' reports the exp() values instead of the original values
#' 
#' The function skips by default all spline variables since this becomes very complicated
#' and there is no simple \deqn{\beta}{beta} to display. For the same reason it skips
#' any interaction variables since it's probably better to display these as a contrast table. 
#' 
#' Note that the rms regression has a separate function that uses the rms:::summaryrms function
#' that returns a matrix that is then pruned.
#' 
#' @param model The regression model
#' @param level The confidence interval level
#' @param remove_interaction_vars Removes the interaction terms as they in 
#'  the raw state are difficult to understand
#' @param remove_strata Strata should most likely not be removed in the crude
#'  version. If you want to force the removal of stratas you can specify the
#'  \code{remove_strata = TRUE}
#' @param remove_cluster Cluster information should most likely also retain
#'  just as the \code{remove_strata} option. Clusters are sometimes used in
#'  cox regression models, \code{\link[survival]{cluster}}
#' @param var_select A vector with regular expressions for choosing what variables
#'  to return (the same format as for the \code{order} argument in 
#'  \code{\link{printCrudeAndAdjustedModel}} call). It can be useful when working with 
#'  large datasets only to report a subsection of all tested variables. This
#'  makes the function both run faster and the data presentation more consice. 
#' @param ... Not used
#' @return Returns a matrix with the columns: 
#'   \code{c("Crude", "2.5 \%", "97.5 \%", "Adjusted", "2.5 \%", "97.5 \%")}.
#'   The row order is not changed from the original model. The percentages can vary depending
#'   on the set level.
#'
#' @seealso \code{\link{printCrudeAndAdjustedModel}}
#' 
#' @example inst/examples/getCrudeAndAdjustedModelData_example.R
#' 
#' @importFrom stringr str_split
#' @import stats
#' 
#' @rdname getCrudeAndAdjustedModelData
#' @export
getCrudeAndAdjustedModelData <- function(model, level=.95, 
                                         remove_interaction_vars = TRUE, 
                                         remove_strata = FALSE,
                                         remove_cluster = FALSE,
                                         var_select,
                                         ...)
  UseMethod("getCrudeAndAdjustedModelData")

#' @export
#' @keywords internal
getCrudeAndAdjustedModelData.default <- function(model, level=.95, 
                                                 remove_interaction_vars = TRUE, 
                                                 remove_strata = FALSE,
                                                 remove_cluster = FALSE, 
                                                 var_select,
                                                 ...){
  
  var_names <- prGetModelVariables(model, 
                                   remove_interaction_vars = remove_interaction_vars,
                                   add_intercept = TRUE)
  if (length(var_names) == 0)
    stop("You have no variables that can be displayed as adjusted/unadjusted.",
      " They are most likely all are part of an interaction, spline, I(),",
      " strata, or some other function.")
  
  # Get the adjusted variables
  adjusted <- prCaDefaultGetCoefAndCI(model, level = level)
  
  # Map rows to variables
  var_rows <- prMapVariable2Name(var_names = var_names, 
                                 available_names = rownames(adjusted),
                                 data = prGetModelData(model))
  
  if (!missing(var_select)){
    greps <- 
      prCaSelectAndOrderVars(names = names(var_rows),
                             order = var_select, 
                             ok2skip = TRUE)
    if (length(greps) == 0)
      stop("The variables that you have tried to select:",
           " '", paste(var_select, collapse="', '"), "'",
           " do not seem to exist - there is no match",
           " for the names: '", paste(names(var_rows), collapse="', '"), "'")
    
    var_rows <- var_rows[sort(greps)]
    var_names <- local({
      tmp <- var_names[var_names %in% names(var_rows)]
      copyAllNewAttributes(from = var_names, to = tmp)
    })
  }
  
  keep <- unlist(lapply(var_rows,
                        function(x) x$location), 
                 use.names = FALSE)
  
  if (length(keep) == 0)
    stop("Error when trying to extract the variable names",
      " from the adjusted values. These names: ", paste(var_names, collapse=", "),
      "\n seem not to exist within the rownames of: ", paste(rownames(adjusted), collapse=", "))
  
  # Sort in order to keep the order
  adjusted <- adjusted[sort(keep), ,drop=FALSE]
  
  unadjusted <- c()
  for(variable in var_names){
    if (!grepl("intercept", variable, 
               ignore.case = TRUE)){
      
      # We should keep any strata information when running the models
      # TODO: Add the nlmn | options
      vars_4_frml <- variable
      if (!is.null(attr(var_names, "strata")) &&
            !remove_strata)
        vars_4_frml <- c(vars_4_frml, attr(var_names, "strata"))

      if (!is.null(attr(var_names, "cluster")) &&
            !remove_cluster)
        vars_4_frml <- c(vars_4_frml, attr(var_names, "cluster"))
      
      frml_4_single_var <- 
        paste(".~", paste(vars_4_frml,collapse="+"))
      
      # Run the same model but with only one variable
      model_only1 <- prEnvModelCall(model, update, frml_4_single_var)
      
      # Get the coefficients processed with some advanced
      # round part()
      new_vars <- prCaDefaultGetCoefAndCI(model_only1, 
                                          level = level,
                                          skip_intercept = TRUE)
      
      # Add them to the previous
      unadjusted <- rbind(unadjusted, new_vars)
    }else{
      # Run the same model but without any variables
      model_only1 <- prEnvModelCall(model, update, ".~1")
      
      # Get the coefficients
      new_vars <- prCaDefaultGetCoefAndCI(model_only1, 
                                          level = level,
                                          skip_intercept = FALSE)
      
      # Add
      unadjusted <- rbind(new_vars, unadjusted)
      
      # Change name back to the original
      rownames(unadjusted)[1] <- variable[grepl("intercept", variable, 
                                                ignore.case = TRUE)]
    }
  }

  if (any(rownames(adjusted) != rownames(unadjusted)))
    stop("The rownames of the adjusted don't match:", 
         "\n\t a:", rownames(adjusted),
         "\n\tUa: ", rownames(unadjusted))
  
  # If just one variable it's not a proper matrix
  if (is.null(dim(adjusted))){
    both <- matrix(c(unadjusted, adjusted), nrow=1)
  }else{
    both <- cbind(unadjusted, adjusted)
  }
  
  levels_str <- c(sprintf("%.1f %%", 100*(1-level)/2),
    sprintf("%.1f %%", 100*(level + (1-level)/2)))
  
  colnames(both) <- c("Crude", levels_str,
                      "Adjusted", levels_str)
  
  attr(both, "model") <- model
  class(both) <- c("getCrudeAndAdjustedModelData", class(both))
  return(both)
}

#' @rdname getCrudeAndAdjustedModelData
#' @export
#' @importFrom Gmisc copyAllNewAttributes
#' 
#' @keywords internal
`[.getCrudeAndAdjustedModelData` <- function(x, i, j, ...){
  ret <- NextMethod()
  attr2skip <- c("dimnames", "dim")
  copyAllNewAttributes(x, ret, attr2skip = attr2skip)
}
