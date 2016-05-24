# Internal functions for the cvAUC package

.process_input <- function(predictions, labels, label.ordering = NULL, folds = NULL, ids = NULL, confidence = NULL){

  # If 'folds' is specified as a list or vector, 'predictions'
  # and 'labels' must be specified as a vector and will be
  # converted into list format programmatically.

  .vec_to_list <- function(idxs, vec){
    return(vec[idxs])
  }

  if (!is.null(folds)){

    if (class(predictions) == "list" | class(labels) == "list"){
      stop("If folds is specified, then predictions and labels must both be vectors.")
    }
    if (length(predictions) != length(labels)) {
      stop("predictions and labels must be equal length")
    }
    
    if (is.vector(folds) && !is.list(folds)){  #If folds is a vector, convert into a list.
      if (length(folds) != length(labels)) {
        stop("folds vector must be the same length as the predictions/labels vectors.")
      } else {
        fids <- as.list(unique(folds))
        folds <- lapply(fids, function(fid, folds){which(folds == fid)}, folds)
      }
    } else if (!is.list(folds)){
      stop("If specifying the folds argument, folds must be a list\n of vectors of indices that correspond to each CV fold or a vector of fold numbers\n the same size as the predictions/labels vectors.")
    } else if (length(unlist(folds)) != length(labels)) {
      stop("Number of observations in the folds argument does not equal number of predictions/labels.")
    }

    #Convert predictions and labels vector arguments into lists.
    predictions <- sapply(folds, .vec_to_list, vec = predictions)
    labels <- sapply(folds, .vec_to_list, vec = labels)
    if (length(labels) > length(unlist(labels))){
      stop("Number of folds cannot exceed the number of observations.")
    }
  }
  # Might add more checking for AUC function; Something like this:
  #else {
  #  if (length(unique(unlist(labels))) != 2) {
  #    stop("AUC only implemented for binary response")
  #  }
  #  if (length(predictions) != length(labels)) {
  #    stop("predictions and labels must be equal length")
  #  }
  #}

  #Process predictions/labels arguments using ROCR prediction function
  pred <- prediction(predictions = predictions, labels = labels, 
    label.ordering = label.ordering)
  predictions <- pred@predictions
  labels <- pred@labels

  #For ci.pooled.cvAUC only, process ids argument
  if (!is.null(ids)){  
    
    if (is.list(ids)){
      if (length(unlist(ids)) != length(unlist(labels))){
        stop("ids must contain same number of observations as predictions/labels.")
      }
    } else if (is.vector(ids)){  #Convert ids vector to list
      if (is.null(folds)){  #Single fold
        ids <- list(ids)
      } else {  #CV folds
        ids <- sapply(folds, .vec_to_list, vec=ids)
      }
    } else if (is.matrix(ids) | is.data.frame(ids)){
      ids <- as.list(data.frame(ids))
    } else {
      stop("Format of ids is invalid.")
    }
    
    #Check to make sure common ids are in the same fold
    if (length(ids) > 1){
      n_ids <- sum(sapply(ids, function(i){length(unique(i))}))
      if (length(unique(unlist(ids))) != n_ids){
        warning("Observations with the same id are currently spread across multiple folds.\nAll observations with the same id must be in the same fold to avoid bias.")
      }
    }
  }

  #Check confidence value
  if (!is.null(confidence)){
    if (is.numeric(confidence) && length(confidence)==1){
      if (confidence <= 0 | confidence >= 1) {
        stop("confidence value must fall within (0,1)")
      }
    }
  }
  
  return(list(predictions = predictions, labels = labels, folds = folds, ids = ids))
}
