#####################################################
# predicted_results
# Calculates results for a model and a set of parameters exactly as anneal or
# likeli would do it, thus allowing the user to find out what the "predicted
# results" would be.
#
# Arguments:
# model = model function for which to calculate predicted results.
# Arguments to this function will be provided from par and source_data.
# par = list of parameters.  If this were simulated annealing,
# these would be the varying parameters.  Doesn't make too
# much sense in this context but this makes it match anneal.
# Each par element name matches an argument in a function
# (either model or pdf). All values listed in par must be
# numeric vectors.
# var = list of other variables and data needed by the model
# and pdf functions, any type as needed.
# source_data = data frame with dependent variable and
# associated independent variables
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
######################################################
predicted_results <- function(model, par, var, source_data, ...) {

  ##
  ## Error checking
  ##
  if (!is.function(model)) {
    stop("predicted_results: model is not a function.\n")
  }
  if (!is.list(par)) {
    stop("predicted_results: par is not a list.\n")
  }
  if (!is.list(var)) {
    stop("predicted_results: var is not a list.\n")
  }
  if (!is.data.frame(source_data)) {
    stop("predicted_results: source_data is not a data frame.\n")
  }
  if (length(source_data) < 1) {
    stop("predicted_results: source_data contains no data.\n")
  }
  if (length(source_data[[1]]) < 1) {
    stop("predicted_results: source_data contains no data.\n")
  }
  if (!is.function(pdf)) {
    stop("predicted_results: pdf is not a function.\n")
  }

  # Here's where we'll put results of all pre-evaluations
  eval_results <- NULL

  # This is where to put the predicted
  predicted <- NULL
  
  # Process the model call

  datasets<-NULL
  datasets[[1]]<-list(value=source_data, varname="source_data")
  
  # Put var and par together
  var <- c(par,var)

  # Get the list of arguments for model
  model_call<-analyze_function(model, list(value=var, varname="var"), NULL, datasets, ...)

  ##
  ## Perform the calculations
  ##

  # First:  any sub-functions needed by the model
  if (!is.null(model_call$pre_eval)) {
    level<-1
    pre_evals<-NULL
    pre_evals[[level]]<-model_call$pre_eval
    curr_index<-1
    max_index<-NULL
    max_index[[level]]<-length(model_call$pre_eval)

    while (level > 0) {
      if (is.null(pre_evals[[level]][[curr_index[[level]]]]$pre_eval)) {
        eval_results[[pre_evals[[level]][[curr_index[[level]]]]$parname]]<-do.call(pre_evals[[level]][[curr_index[[level]]]]$fun, pre_evals[[level]][[curr_index[[level]]]]$call)
        curr_index[[level]]<-curr_index[[level]]+1
      
        if (curr_index[[level]] > max_index[[level]]) {
          # We've evaluated all the pre-evals in this level - back up
          # one level
          level <- level - 1
        }
      }  
      else {
        # The current level has pre-evals - before evaluating the current level,
        # go down a level to the first of its pre-evals
        level <- level + 1
        curr_index[[level]] <- 1
        max_index[[level]] <- length(pre_evals[[level]][[curr_index[[level]]]]$pre_eval)
        pre_evals[[level]] <- pre_evals[[level]][[curr_index[[level]]]]$pre_eval
      }        
    }    
  }
  eval_results <- as.list(eval_results)

  # Second:  the model predicted values
  predicted <- do.call(model,model_call$call)

  if (any(is.infinite(predicted)) && any(is.nan(predicted))) {
    stop("predicted_results: the model produced math errors.\nLikelihood cannot be calculated.\n")
  }

  predicted
}
