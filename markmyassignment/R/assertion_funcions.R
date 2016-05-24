#' @title
#' Check that tasks exist in assignment.
#' 
#' @description
#' Checks that tasks that are inputted in mark_my_assignment exists. If not, a warning is produced.
#' 
#' @param tasks
#' The \code{task} vector from \code{mark_my_assignment}.
#' @param path
#' Path to assignment file. Passed to \code{read_assignment_yml}.
#' @return
#' A warning message or nothing.
#' 
#' @keywords internal
#' 
check_existance_tasks <- function(tasks, path = NULL){
  res <- read_assignment_yml(path = path)
  if(!all(tasks %in% names(res$tasks))){
    warning(paste("The following tasks do not exist:", paste(
      tasks[!(tasks %in% names(res$tasks))], collapse = ", ")))
  }  
}



#' @title
#' Checks the input arguments for mark_my_assignment and mark_my_file
#' 
#' @description
#' Checks that the input arguments are of the correct type.
#' 
#' @param tasks
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @param mark_file
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @param lab_file
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @param force_get_tests
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @param quiet
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @param reporter
#' Arguments from \code{mark_my_assignment} and \code{mark_my_file}. 
#' @return
#' If all inputs are OK, nothing. Otherwise the functions stop.
#' 
#' @keywords internal
assert_function_arguments_in_API <- function(
  tasks, mark_file, lab_file = NULL, force_get_tests, quiet, reporter){
  
  if(!is.null(tasks) & !is.character(tasks))
    stop("Tasks must be a character vector or NULL.")
  
  if(!is.null(mark_file))
    if(!file.exists(mark_file))
      stop("Mark file does not exist.")
  
  if(!is.null(lab_file))
    if(!file.exists(lab_file))
      if(!url_success(lab_file))
        stop("Lab file could not be found.")
  
  if(!is.logical(force_get_tests) | !is.logical(quiet))
    stop("force_get_tests and quiet must be logical.")
  
  if(!is.null(reporter) & !is.character(reporter))
    stop("reporter must be character.")
  
}

