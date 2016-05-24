#' Change the name of elements that an experiment contains
#'
#' This function change the name of problems, methods or parameter values that
#' an existing experiment object contains.
#'
#' @export
#' @param e Input experiment
#' @param elements A list of arrays of strings containing the new names. The old
#' name will be specified as the name of the element in such array, and the name
#' for the parameter, method or problem will be given by the name of the 
#' corresponding object in the list.
#' If a name is not present in the set of parameter names or parameter values,
#' it will be ignored.
#' @param name The name of the new experiment. If NULL, the previous name will
#' be used.
#' @return A modified exreport experiment object with some changes on the name of
#' the elements.
#'
#' @examples
#' # We load the wekaExperiment problem as an experiment and then change the name
#' # of one value for the parameter discretization and for one method.
#' 
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' expRename(experiment, list(featureSelection = c("no"="false"),
#'                            method=c("RandomForest"="RndForest")))
#' 
expRename  <- function(e, elements=list(), name = NULL){
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  if (!is.list(elements))
    stop(.callErrorMessage("wrongParameterError", "elements", "non-empty list"))
  # Check that all elements have a proper name
  if (length(elements)>0 && is.null(names(elements)))
    stop(.callErrorMessage("noNamesError"))
  oldMethods    <- levels(e$data[[e$method]])
  oldProblems   <- levels(e$data[[e$problem]])
  oldParameters <- list()
  for(param in e$parameters){
    oldParameters[[param]] <- levels(e$data[[param]])
  }

  # Copy the experiment
  result <- e
  
  # Now we apply all the renames
  for(var in names(elements)){
    # If the variable name does not exist, we continue
    if(!(var %in% names(e$data)))
      next
    l <- levels(e$data[[var]])
    for(val in names(elements[[var]])){
      idx <- which(l==val)
      l[idx] <- elements[[var]][[val]]
    }
    levels(result$data[[var]]) <- l
  }
  
  if(!is.null(name)){
    result$tags$title = name
    result$tags$context = name
  }

  # Append this operation in the historic
  varNames <- names(elements)
  oldElemNames <- lapply(elements, FUN=names)
  newElemNames <- lapply(elements, FUN=identity)
  renames <- c()
  for(i in varNames){
    renames <- c(renames,paste0(i,": [",paste(oldElemNames[[i]],newElemNames[[i]],sep="->",collapse=", "), "]"))
  }
  result$historic <- c(result$historic, 
                       list(paste0("From experiment renamed to ",result$name, ", discrete values from method, problem or parameters columns have been renamed: ",
                                   paste(renames,collapse = "; "))))

  result
}