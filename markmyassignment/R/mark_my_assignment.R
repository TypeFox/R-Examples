#' @title
#' Mark assignment in global environment
#' 
#' @description
#' Mark assignment in global environment.
#' 
#' @param tasks
#'   Which task should be corrected (if more than one). Default is all. 
#'   To see the different task, see \code{\link{show_tasks}}.
#' @param mark_file
#'   Argument is deprecated, use mark_my_file instead.
#' @param force_get_tests
#'   Force download of test files before marking of assignments. Default is FALSE.
#' @param quiet
#'   Should test be run without output?
#' @param reporter to use. Default is the 'summary' or specified in assignment yml file.
#' 
#' @examples
#' assignment_path <- 
#'  paste0(system.file(package = "markmyassignment"), "/extdata/example_assignment01.yml")
#' set_assignment(assignment_path)
#' source(paste0(system.file(package = "markmyassignment"), "/extdata/example_lab_file.R"))
#' mark_my_assignment()
#' 
#' @export
mark_my_assignment <- function(tasks = NULL, mark_file = NULL, force_get_tests = FALSE, quiet = FALSE, reporter = NULL){
  assert_function_arguments_in_API(
    tasks = tasks, mark_file = mark_file, force_get_tests = force_get_tests,
    quiet = quiet, reporter = reporter)
  if(!is.null(mark_file)){
    .Deprecated("mark_my_file", old = "mark_file")
  }
  get_tests(tasks = tasks, force_get_tests = force_get_tests)
  if(is.null(reporter)) reporter <- get_mark_my_reporter()
  test_results <- run_test_suite("mark_my_assignment", tasks, mark_file, quiet, reporter = reporter)
  test_results_df <- as.data.frame(test_results) 
  if(!any(test_results_df$error) & sum(test_results_df$failed) == 0 & is.null(tasks) & !quiet) cheer()
  check_existance_tasks(tasks = tasks)
  return(invisible(test_results))
}

#' @title
#' Mark assignments in a directory
#' 
#' @description
#' Marks assignments in a directory. Stores the results.
#' 
#' @param directory
#'   Directory with assignments files.
#' @param lab_file
#'   Assignment file to set before marking the assignment (url or local path).
#' @param tasks
#'   Which task should be corrected (if more than one). 
#'   Default is all. To see the different task, see \code{\link{show_tasks}}.
#' @param force_get_tests
#'   Force download of test files before marking of assignments. Default is FALSE.
#'   
#' @keywords internal
#'   
#' @export
mark_my_dir <- function(directory, lab_file, tasks = NULL, force_get_tests = FALSE){
  file_names <- dir(directory, pattern = "\\.[Rr]")
  if(length(file_names) == 0) stop("No files to mark.")
  files_to_mark <- paste0(directory, "/", file_names)
  res_mark <- vector(mode = "list", length = length(files_to_mark))
  names(res_mark) <- file_names
  
  for(i in seq_along(files_to_mark)){
    res_mark_temp <- try(
      mark_my_file(tasks = tasks, 
                   mark_file = files_to_mark[i],
                   assignment_path = lab_file,
                   force_get_tests = force_get_tests, 
                   quiet = TRUE), silent=TRUE)
    force_get_tests <- FALSE
    if(class(res_mark_temp) == "try-error") {
      message(res_mark_temp[1])
      message(file_names[i], " could not be marked.")
      res_mark[[i]] <- as.character(res_mark_temp[1])
    } else {
      res_mark[[i]] <- res_mark_temp
      print(paste(file_names[i], "was marked."))
    }
  }
  return(res_mark)
}



#' @title
#' Get test files
#' 
#' @description
#' Downloads the test files for the current assignment and save them to 
#' temp directory.
#' 
#' @param tasks
#'   Which task should be downloaded. Default is "all".
#' @param force_get_tests
#'   Force download/get test (ignore cached tests).
#' 
#' @keywords internal
#' 
get_tests <- function(tasks = NULL, force_get_tests = FALSE){
  assignment <- read_assignment_yml()
  dir.create(path = mark_my_test_dir(), recursive = TRUE, showWarnings = FALSE)
  
  tasks_to_get <- names(assignment$tasks)
  if(!is.null(tasks)) tasks_to_get <- tasks_to_get[tasks_to_get %in% tasks]
  if(!force_get_tests) tasks_to_get <- tasks_to_get[!tasks_to_get %in% cached_tasks()]
  
  for(task in tasks_to_get) {
    for(i in seq_along(assignment$tasks[[task]]$url)){
      dest <- paste0(mark_my_test_dir(), "/test-", task, "-", i, ".R")
      path <- path_type(assignment$tasks[[task]]$url[i]) 
      get_file(path = path, dest = dest)
    }
  }
  
  if(force_get_tests | !"00mandatory" %in% cached_tasks()){
    for(i in seq_along(assignment$mandatory$url)){
      dest <- paste0(mark_my_test_dir(), "/test-00mandatory-", i, ".R")
      path <- path_type(assignment$mandatory$url[i]) 
      get_file(path = path, dest = dest)
    }
  }
  return(invisible(TRUE))
}

#' @title
#' Cached tasks
#' 
#' @description
#'   Checks which assignments that are cached (ie already downloaded to temp dir).
#' 
#' @return
#'   character vector with cached assignments.
#' 
#' @keywords internal
#' 
cached_tasks <- function(){    
  files <- dir(mark_my_test_dir())
  unique(unlist(lapply(strsplit(files, split = "-"), FUN=function(X) X[2])))
}


#' @title
#'   Run test suite
#' 
#' @description
#'   Runs test on the tasks. Always run mandatory tests.
#' 
#' @param caller
#'   Either "mark_my_assignment" or "mark_my_file"
#' @param tasks
#'   Which task should be tested
#' @param mark_file
#'   Run tests on a R-file. Default is NULL means global environment.
#' @param quiet
#'   Should the output be supressed (only returning test results)
#' @param reporter
#'   Reporter to use. Standard is student.
#'
#' @return
#'   test_suite results
#'   
#' @keywords internal
#'   
run_test_suite <- function(caller, tasks = NULL, mark_file = NULL, quiet = FALSE, reporter = "summary"){
  
  test_directory <- mark_my_test_dir()
  
  if(caller == "mark_my_assignment" & is.null(mark_file))
    mark_my_env <- new.env(parent = .GlobalEnv)
  else
    mark_my_env <- new.env(parent = parent.env(env = .GlobalEnv))
  
  if(!is.null(mark_file)){
    mark_file <- delete_circular_calls(mark_file)
    tf_path <- tempfile(pattern = "mark_file", fileext = ".txt")
    writeLines(text = mark_file, con = tf_path)
    source(file = tf_path, local = mark_my_env)
    unlink(x = tf_path)
  } 
  
  if(quiet) reporter <- "silent"
  
  if(is.null(tasks)) tasks <- "all" else tasks <- paste(c("00mandatory", paste(tasks, collapse="|")), collapse="|")
  
  if(tasks == "all") tasks <- NULL
  test_res <- test_dir(path = test_directory, 
                       filter = tasks, 
                       reporter = reporter, env = mark_my_env)
  test_res
}


#' @title
#'  Functions to create directories
#'  
#' @description
#'  Functions to create directories
#'  
#' @name directories
#' 
#' @keywords internal
#' 
mark_my_base_dir <- function() paste0(tempdir(), "/markmyassignment")

#' @rdname directories
#' @param no assignment number
#' @keywords internal
mark_my_assignment_dir <- function(no = 1) paste0(mark_my_base_dir(), "/assignment", no)

#' @rdname directories
#' @param ... to send to \code{\link{mark_my_assignment_dir}}
#' @keywords internal
mark_my_test_dir <- function(...) paste0(mark_my_assignment_dir(...), "/tests")



#' @title
#'  Cheer when all tasks pass
#'  
#' @description
#' Cheer when all tasks pass
#' 
#' @keywords internal
#' 
cheer <- function() {
  cat(sample(x = c("Yay! All done!",
                   "Good work!",
                   "You're a coding rockstar!",
                   "Keep up the good work!",
                   "Everything's correct!"), 1))
}

#' @title
#'  Get reporter from yml file
#'  
#' @description
#'  Get reporter from yml file
#'  
#'  Default reporter is 'summary'. 
#'  
#' @keywords internal
#'  
get_mark_my_reporter <-function(){
  assign_yml <- read_assignment_yml()
  if("reporter" %in% names(assign_yml)){
    reporter <- assign_yml$reporter
    output <- capture_output(
      check_reporter <- try(test_file(path = file.path(system.file(package = "markmyassignment"), "extdata/test_reporter_file.R"), reporter = reporter), silent = TRUE)
    )
    if(inherits(check_reporter, what = "try-error")) {
      warning("Reporter '", reporter, "' not found. Summary reporter is used.")
      reporter <- "summary"
      }
  } else {
    reporter <- "summary"
  }
  reporter
}

#' @title
#'  Checks and deletes circular calls 
#'  
#' @description
#'  Checks and deletes circular calls 
#'  
#' @param mark_file File to check
#'  
#' @return
#'  Character vector of the possibly changed mark file
#'  
#' @keywords internal
#' 
delete_circular_calls <- function(mark_file){
  txt_in <- txt_out <- as.character(parse(mark_file))
  forbidden <- c(
    "mark_my_assignment", "mark_my_dir", "set_assignment", "mark_my_file",
    "install.packages", "utils::install.packages",
    "devtools::install_github", "install_github", "data", "system")
  regex <- paste("(^|;| )", forbidden, "\\(.*\\)", sep = "")
  for(pattern in regex){
    txt_out <- gsub(pattern = pattern, replacement = "", x = txt_out)
  }
  if(!identical(txt_in, txt_out))
    message("The following statements were ignored when running mark_my_file:\n",
            paste(txt_in[txt_in != txt_out], collapse = "\n"))
  
  #   indices <- grep(pattern = "^data\\(.*\\)", x = txt_out, value = F)
  #   txt_out[indices] <- gsub(pattern = "^data\\(", replacement = "markmyassignment:::data_mma\\(", x = txt[indices])
  
  return(txt_out)
}


