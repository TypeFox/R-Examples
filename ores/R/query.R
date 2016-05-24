ores_query <- function(path, ...){
  
  url <- paste0("http://ores.wmflabs.org/", path)
  ua <- httr::user_agent("ORES R Client - https://github.com/Ironholds/ores")
  result <- httr::GET(url, ua, ...)
  httr::stop_for_status(result)
  return(httr::content(result))
}

#'@title List Supported Projects
#'@description \code{\link{list_wikis}} lists Wikimedia
#'projects that support some or all of the ORES models.
#'
#'@inheritParams list_models
#'
#'@export
list_wikis <- function(...){
  
  result <- ores_query("scores/", ...)
  return(unlist(result$contexts))
}

#'@title List Model Information
#'@description \code{\link{list_models}} lists information about
#'the models for a particular wiki, including what models are available,
#'how they have been trained, information about the model's accuracy and
#'ROC, and the model's version.
#'
#'@param project a Wikimedia project. Supported projects can be obtained with
#'\code{\link{list_wikis}}.
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@examples
#'# Get model information for the English-language Wikipedia
#'model_data <- list_models("enwiki")
#'
#'@seealso \code{\link{list_wikis}} for retrieving the list of supported projects,
#'and \code{\link{check_reverted}} and similar for actual checking
#'against models.
#'
#'@export
list_models <- function(project, ...){
  result <- ores_query(paste0("scores/", project, "/", ...))
  return(result$models)
}

#'@title Check Revert Probabilities
#'@description \code{check_reverted} identifies
#'whether or not an edit is considered likely, by
#'the ORES models, to be reverted.
#'
#'@param edits a vector of edit IDs, as integers.
#'
#'@return A data.frame of five columns; \code{edit}, the
#'edit ID, \code{project}, the project, \code{prediction},
#'whether the model predicts that the edit will be reverted,
#'\code{false_prob}, the probability that the model's prediction
#'is wrong, and \code{true_prob}, the probability that the model's
#'prediction is correct. In the event of an error (due to the edit
#'not being available) NAs will be returned in that row.
#'
#'@examples
#'# A simple, single-diff example
#'revert_data <- check_reverted("enwiki", 34854345)
#'
#'@seealso
#'\code{\link{check_goodfaith}} to identify if a set of edits were made
#'in good faith, and \code{\link{check_damaging}} to check if a set of edits
#'were damaging.
#'
#'@inheritParams list_models
#'@export
check_reverted <- function(project, edits, ...){

  result <- lapply(edits, function(id, project, ...){
    
    data <- ores_query(path = paste0("scores/", project, "/reverted/", id),
                       ...)
    if("error" %in% names(data[[1]])){
      return(data.frame(edit = names(data),
                        project = project,
                        prediction = NA,
                        false_prob = NA,
                        true_prob = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(edit = names(data),
                      project = project,
                      prediction = data[[1]]$prediction,
                      false_prob = data[[1]]$probability$false,
                      true_prob = data[[1]]$probability$true,
                      stringsAsFactors = FALSE))
    
  }, project = project)
  
  return(do.call("rbind", result))
}

#'@title Check Good-Faith Probability
#'@description \code{check_goodfaith} identifies whether
#'or not an edit was made in 'good faith' - whether it was well-intentioned,
#'even if it is not a high-quality contribution.
#'
#'@return A data.frame of five columns; \code{edit}, the
#'edit ID, \code{project}, the project, \code{prediction},
#'whether the model predicts that the edit was made in good faith,
#'\code{false_prob}, the probability that the model's prediction
#'is wrong, and \code{true_prob}, the probability that the model's
#'prediction is correct. In the event of an error (due to the edit
#'not being available) NAs will be returned in that row.
#'
#'@examples
#'# A simple, single-diff example
#'goodfaith_data <- check_goodfaith("enwiki", 34854345)
#'
#'@seealso
#'\code{\link{check_reverted}} to identify if a set of edits are likely
#'to be reverted, and \code{\link{check_damaging}} to check if a set of edits
#'were damaging.
#'
#'@inheritParams check_reverted
#'@export
check_goodfaith <- function(project, edits, ...){
  
  result <- lapply(edits, function(id, project, ...){
    
    data <- ores_query(path = paste0("scores/", project, "/goodfaith/", id),
                       ...)
    if("error" %in% names(data[[1]])){
      return(data.frame(edit = names(data),
                        project = project,
                        prediction = NA,
                        false_prob = NA,
                        true_prob = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(edit = names(data),
                      project = project,
                      prediction = data[[1]]$prediction,
                      false_prob = data[[1]]$probability$false,
                      true_prob = data[[1]]$probability$true,
                      stringsAsFactors = FALSE))
    
  }, project = project)
  
  return(do.call("rbind", result))
}

#'@title Check Damaging Probability
#'@description \code{check_damaging} identifies whether
#'or not an edit was damaging - the type that caused actual
#'harm to an article.
#'
#'@return A data.frame of five columns; \code{edit}, the
#'edit ID, \code{project}, the project, \code{prediction},
#'whether the model predicts that the edit was damaging,
#'\code{false_prob}, the probability that the model's prediction
#'is wrong, and \code{true_prob}, the probability that the model's
#'prediction is correct. In the event of an error (due to the edit
#'not being available) NAs will be returned in that row.
#'
#'@examples
#'# A simple, single-diff example
#'damaging_data <- check_damaging("enwiki", 34854345)
#'
#'@seealso
#'\code{\link{check_goodfaith}} to identify if a set of edits were made
#'in good faith, and \code{\link{check_reverted}} to check if a set of edits
#'are likely to be reverted.
#'
#'@inheritParams check_reverted
#'@export
check_damaging <- function(project, edits, ...){
  
  result <- lapply(edits, function(id, project, ...){
    
    data <- ores_query(path = paste0("scores/", project, "/damaging/", id),
                       ...)
    if("error" %in% names(data[[1]])){
      return(data.frame(edit = names(data),
                        project = project,
                        prediction = NA,
                        false_prob = NA,
                        true_prob = NA,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(edit = names(data),
                      project = project,
                      prediction = data[[1]]$prediction,
                      false_prob = data[[1]]$probability$false,
                      true_prob = data[[1]]$probability$true,
                      stringsAsFactors = FALSE))
    
  }, project = project)
  
  return(do.call("rbind", result))
}