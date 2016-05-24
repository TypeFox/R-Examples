#####################################################
# likeli
# Calculates likelihood.
# Arguments:
# model = model function for which to calculate likelihood.
# Arguments to this function will be provided from par
# and source_data.
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
# pdf = probability density function.  To calculate
# negative log likelihood, use a function such as
# dnorm that can calculate log probability.
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
######################################################

likeli <- function(model, par, var, source_data, pdf, ...) {

  ##
  ## Error checking
  ##
  if (!is.function(model)) {
    stop("likeli: model is not a function.\n")
  }
  if (!is.list(par)) {
    stop("likeli: par is not a list.\n")
  }
  if (!is.list(var)) {
    stop("likeli: var is not a list.\n")
  }
  if (!is.data.frame(source_data)) {
    stop("likeli: source_data is not a data frame.\n")
  }
  if (length(source_data) < 1) {
    stop("likeli: source_data contains no data.\n")
  }
  if (length(source_data[[1]]) < 1) {
    stop("likeli: source_data contains no data.\n")
  }
  if (!is.function(pdf)) {
    stop("likeli: pdf is not a function.\n")
  }
  # Make sure all values in par are numeric
  for (i in 1:length(par)) {
    if (!is.numeric(par[[i]])) {
      stop("likeli: All values in par must be numeric.\n")
    }
    if (!is.vector(par[[i]])) {
      stop("likeli: All values in par must be vectors.\n")
    }
  }

  # Here's where we'll put results of all pre-evaluations
  eval_results <- NULL

  # This is where to put the predicted
  predicted <- NULL
  
  # Put together var and par
  var <- c(par,var)

  # Process the PDF call
  datasets<-NULL
  datasets[[1]]<-list(value=source_data, varname="source_data")
  pdf_call<-analyze_function(pdf, list(value=var, varname="var"), NULL, datasets, ...)

  # Process the model call

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
    stop("likeli: the model produced math errors.\nLikelihood cannot be calculated.\n")
  }

  # Third: any sub-functions needed by the pdf
  if (!is.null(pdf_call$pre_eval)) {
    level<-1
    pre_evals<-NULL
    pre_evals[[level]]<-pdf_call$pre_eval
    curr_index<-1
    max_index<-NULL
    max_index[[level]]<-length(pdf_call$pre_eval)

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

  # Fourth: sum over the PDF
  lh<-sum (do.call(pdf,pdf_call$call))

  lh
}
