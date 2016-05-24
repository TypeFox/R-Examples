#' @title
#' Set assignment to mark
#' 
#' @description
#' Sets the assignment to mark. Behind the scenes it download the test suite for the assignment.
#' 
#' @param path
#' Path to the yml file
#' @param auth_token 
#' Authorization token (for github). Not implemented.
#' 
#' @examples
#' assignment_path <- path <- 
#'   paste0(system.file(package = "markmyassignment"), "/extdata/example_assignment01.yml")
#' set_assignment(assignment_path)
#' 
#' @export
set_assignment <- function(path, auth_token = NULL){
  path <- path_type(path)
  if(inherits(path, what = "path_error")) stop("Assignment path/url does not work.")
  temp_folder_check_create()
  temp_file <- tempfile()
  on.exit(unlink(temp_file))
  
  if(file.exists(mark_my_assignment_dir())) unlink(mark_my_assignment_dir(), recursive = TRUE, force = TRUE)
  dir.create(mark_my_assignment_dir(), recursive = TRUE, showWarnings = FALSE)
  dest <- paste0(mark_my_assignment_dir(), "/assignment1.yml")
  get_file(path, temp_file)  
  assignment_yml_ok(path = temp_file)
  file.copy(from = temp_file, to = dest, overwrite = TRUE)
  assignment <- read_assignment_yml()
  if("packages" %in% names(assignment))
    check_installed_packages(assignment$packages)
  message("Assignment set:\n", assignment$name, " : ", assignment$description)
  invisible(dest)
}

#' @title
#'  Check and create folder if missing.
#' 
#' @description
#'   Checks if markmyassignment folder exist in R temp directory.
#'   If not, the folder is created.
#'   
#' @keywords internal
#' 
temp_folder_check_create <- function() {
  if(!"markmyassignment" %in% dir(tempdir())){
    dir.create(path = mark_my_base_dir())
  }
}

#' @title
#'   Check the yml file to be a correct assignment yml.
#' 
#' @description
#'   Check the yml file that it is correct and can be used. Otherwise warn.
#' 
#' @param path \code{path} object from \code{\link{path_type}}
#' 
#' @return 
#'   boolean 
#'   
#' @keywords internal
#' 
assignment_yml_ok <- function(path = NULL){
  assignment <- try(read_assignment_yml(path), silent = TRUE)
  if(inherits(assignment, "try-error")) return(FALSE)  
  check_assignment_file(assignment)
}

#' @title
#' Check assignment yml file that it is a correct assignment file
#' 
#' @description
#' Check assignment yml file that it is a correct assignment file
#' @param assignment list to test.
#' 
#' @return 
#'   boolean 
#'   
#' @keywords internal
#' 
check_assignment_file <- function(assignment){
  # The yml contain at most 6 slots.
  check <- all(names(assignment) %in% c("name", "description", "reporter", "tasks", "mandatory", "packages"))
  if(!check) stop("Assignment file contain erroneous parts (except name, desc., reporter, tasks and mandatory.")
  # The name and description is of length 1
  check <- all(unlist(lapply(assignment[c("name", "description")], length)) == 1)
  if(!check) stop("Name/description is not of length 1 in assignment file.")
  # Check that all url exists/works
  urls <- try(as.list(unlist(lapply(assignment[["tasks"]], FUN = function(X) return(X$url)))), silent = TRUE)
  if(inherits(urls, "try-error")) stop("All tasks in assignments file do not contain urls")
  check <- !any(unlist(lapply(lapply(urls, path_type), class)) == "path_error")
  if(!check) stop("Not all tasks in assignments have working urls.")
  
  if("mandatory" %in% names(assignment)) {
    # Check mandatory urls
    urls <- try(as.list(assignment[["mandatory"]]$url), silent = TRUE)
    if(inherits(urls, "try-error")) stop("Mandatory urls are missing.")
    check <- !any(unlist(lapply(lapply(urls, path_type), class)) == "path_error")
    if(!check) stop("Mandatory urls are not working.")
  }
  
  TRUE
}

#' @title
#' Get the path type.
#' 
#' @description
#' Check the path type. 
#' 
#' @param path Character element of url or local search path.
#' 
#' @return path type as character element c("path_local", "path_http", "path_error")
#' 
#' @keywords internal
#' 
path_type <- function(path){
  if(file.exists(path)){
    class(path) <- c("path_local", "character")
  } else {
    try_http <- try(identical(httr::status_code(httr::HEAD(path)), 200L), silent = TRUE)
    if (!is(try_http, "try-error") && try_http){
      class(path) <- c("path_http", "character")
    } else {
      class(path) <- c("path_error", "character")
    }
  }  
  path
}

#' @title
#' Get the file from the path
#' 
#' @description
#' Get/download the file from the path.
#' 
#' @param path
#'   Path object
#' @param dest
#'   Destination for the file
#' @param ...
#'   Further arguments to send to \code{httr::GET()}.
#' 
#' @keywords internal
#' 
get_file <- function(path, dest, ...){
  stopifnot(!inherits(path, "path_error"))
  UseMethod("get_file")
}

get_file.path_local <- function(path, dest, ...){
  file.copy(from = path, to = dest, overwrite = TRUE)
}

get_file.path_http <- function(path, dest, ...){
  request <- httr::GET(path, ...)
  httr::stop_for_status(request)
  writeBin(httr::content(request, "raw"), dest)
}

#' @title
#' Load assignment information
#' 
#' @description
#' Check if there exist an assignmentfile and then load it.
#' 
#' @param path \code{path object} from \code{\link{path_type}}
#' 
#' @return assignment object
#' 
#' @keywords internal
#' 
read_assignment_yml <- function(path = NULL){
  if(is.null(path)){
    assignment_file <- paste0(mark_my_assignment_dir(), "/assignment1.yml")
  } else {
    assignment_file <- path
  }
  if(file.exists(assignment_file)){
    res <- suppressWarnings(yaml::yaml.load_file(assignment_file))
    if(is.list(res)) return(res) else stop("Not a correct .yml file")
  } else {
    stop("No assignment has been set. Please use set_assignment().", call. = FALSE)
  }
}


#' @title
#'   Get the name of the tasks in the assignment.
#' 
#' @description
#'   Get the name of the tasks in the assignment.
#'   
#' @examples
#' # We first set the assignment
#' assignment_path <- 
#'  paste0(system.file(package = "markmyassignment"), "/extdata/example_assignment01.yml")
#' set_assignment(assignment_path)
#'  
#' show_tasks()
#' 
#' @export
show_tasks <- function(){
  assignment <- read_assignment_yml()
  names(assignment$task)
}


#' @title
#'  Check whether required packages are installed and loaded.
#' 
#' @description
#'   Checks if the packages listed in assignment file are loaded and installed.
#'   If not, a warning message is printed.
#' @param packages
#'   Packages to check
#'   
#' @keywords internal
#'   
check_installed_packages <- function(packages) {
  
  if(all(paste("package:", packages, sep="") %in% search())){
    # All packages are loaded and installed
  }else{
    if(all(packages %in% rownames(utils::installed.packages()))){
      warning("The following packages need to be loaded:\n",
              paste(packages[!paste("package:", packages, sep="") %in% search()], collapse = ", "))
    }
    else{
      warning("The following packages need to be installed and then loaded:\n",
              paste(packages[!packages %in% rownames(utils::installed.packages())], collapse=", "))
    }
  }
}

