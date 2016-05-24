#' @title
#' Mark assignment file
#' 
#' @description
#' Mark a specific assignment file
#' 
#' @param tasks
#'   Which task should be corrected (if more than one). Default is all. 
#'   To see the different task, see \code{\link{show_tasks}}.
#' @param mark_file
#'   Path to the file to mark.
#' @param assignment_path
#'   Assignment file to set before marking the assignment (url or local path).
#' @param force_get_tests
#'   Force download of test files before marking of assignments. Default is FALSE.
#' @param quiet
#'   Should test be run without output?
#' @param reporter to use. Default is the 'summary' or specified in assignment yml file.
#' 
#' @examples
#' assignment_path <- 
#'   paste0(system.file(package = "markmyassignment"), "/extdata/example_assignment01.yml")
#' file_path <- paste0(system.file(package = "markmyassignment"), "/extdata/example_lab_file.R")
#' mark_my_file(mark_file = file_path, assignment_path = assignment_path)
#' 
#' @export
mark_my_file <- function(tasks = NULL, mark_file=file.choose(), assignment_path, force_get_tests = FALSE, quiet = FALSE, reporter){
  
  if(!missing(assignment_path)) {
    old_warn_opt <- options(warn = 2)
    set_assgn_result <- try(suppressMessages(set_assignment(assignment_path)), silent = TRUE)
    options(warn = old_warn_opt$warn)
    if(is(set_assgn_result, "try-error"))
      stop(set_assgn_result[1])
  }
  if(missing(reporter)) reporter <- get_mark_my_reporter()
  
  assert_function_arguments_in_API(
    tasks = tasks, mark_file = mark_file, lab_file = assignment_path,
    force_get_tests = force_get_tests, quiet = quiet, reporter = reporter)
  
  
  get_tests(tasks = tasks, force_get_tests = force_get_tests)
  test_results <- run_test_suite("mark_my_file", tasks, mark_file, quiet, reporter = reporter)
  test_results_df <- as.data.frame(test_results)
  if(!any(test_results_df$error) & sum(test_results_df$failed) == 0 & is.null(tasks) & !quiet) cheer()
  check_existance_tasks(tasks = tasks)
  return(invisible(test_results))
}

