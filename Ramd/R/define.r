#' Load a bunch of dependencies by filename
#' 
#' This is useful for reducing pollution in the global namespace,
#' and not loading multiple files twice unnecessarily.
#'
#' @export
#' @param ... see examples.
#' @examples
#' \dontrun{
#' helper_fn <- define('some/dir/helper_fn')
#' define(c('some/dir/helper_fn', 'some/other_dir/library_fn'), function(helper_fn, library_fn) { ... }
#' helper_fns <<- define('some/dir/helper_fn1', 'some/otherdir/helper_fn2')
#' helper_fns[[1]]('do something'); helper_fns[[2]]('do something else')
#' }
define <- (function() {

  number_of_required_arguments <- function(fn) {
    function_has_variable_number_of_arguments <- '...' %in% names(formals(fn))
    if (function_has_variable_number_of_arguments) return(NA_real_)
    function_arguments <- formals(fn)
    required_arguments <- sapply(function_arguments, class) == 'name'
    sum(required_arguments)
  }

  process_function_with_no_dependencies <- function(fn) {
    number_of_arguments <- number_of_required_arguments(fn)
    if (number_of_arguments == 0) fn()
    else if (number_of_arguments == 1) fn(define)
    else stop("Ramd::define only processes functions with <= 1 ",
              "arguments if no dependencies are given, but the ",
              "passed function has ", number_of_required_arguments, 
              " required arguments")
  }

  flatten <- function(lists) {    
    atomic_vector <- unlist(c(lists))
    delimited_string <- paste(atomic_vector, collapse = ' ')
    strsplit(delimited_string, '[^a-zA-Z0-9.-_`:\\\\\\/]+')[[1]]
  }

  parse_dependencies <- function(arguments) {
    cd <- current_directory()
    if (any(sapply(arguments, class) != 'character'))
      stop("Ramd::define only accepts atomic character vectors for ",
           "specifying dependencies")
    dependencies <- unlist(c(arguments))
    if ('Ramd.no_flatten' %in% names(.Options) &&
        getOption('Ramd.no_flatten')) dependencies
    else flatten(dependencies)
  }

  fetch_dependencies <- function(arguments) {
    dependency_names <- parse_dependencies(arguments)
    dependencies <- lapply(dependency_names, load_dependency)
    names(dependencies) <- dependency_names
    dependencies
  }

  verify_number_of_required_arguments_matches_number_of_dependencies <-
    function(fn, number_of_dependencies) {
      num_of_required_arguments <- number_of_required_arguments(fn)
      if (is.na(num_of_required_arguments)) return(TRUE)
      if (num_of_required_arguments != number_of_dependencies)
        stop("Ramd::define was not able to load dependencies because ",
             number_of_dependencies, " dependenc",
             # Pluralization, for fun!
             if (number_of_dependencies == 1) 'y was' else 'ies were',
             " passed in but the given function has ",
             num_of_required_arguments, " required argument",
             if (num_of_required_arguments == 1) '' else 's')
      TRUE
    }

  function(...) {
    arguments <- list(...)
    if ('packages' %in% names(arguments)) {
      if (length(arguments) == 1)
        stop("Ramd::define does more than just load packages, ",
             "please provide some dependencies or a function. ",
             "To just load packages, use Ramd::packages")
      packages(arguments$packages)
      arguments <- arguments[names(arguments) != 'packages']
    }

    fn <- tail(arguments, 1)[[1]]
    valid_function <- is.function(fn)
    if (valid_function) {
      dependencies <- head(arguments, -1)
      if (length(dependencies) == 0)
        return (process_function_with_no_dependencies(fn))
    } else dependencies <- arguments

    if (valid_function)
      verify_number_of_required_arguments_matches_number_of_dependencies(
        fn, length(unlist(dependencies)))

    dependencies <- fetch_dependencies(dependencies)
    if (valid_function) do.call(fn, unname(dependencies))
    else dependencies
  }

})()
