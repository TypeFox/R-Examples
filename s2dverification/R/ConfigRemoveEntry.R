ConfigRemoveEntry <- function(configuration, dataset_type, dataset_name = NULL, var_name = NULL, position = NULL) {
  table_name <- dataset_type
  if (!is.null(dataset_name) && !is.null(var_name)) {
    if (dataset_name == ".*") {
      if (var_name == ".*") {
        level <- 1
      } else {
        level <- 3
      }
    } else {
      if (var_name == ".*") {
        level <- 2
      } else {
        level <- 4
      }
    }

    position <- which(unlist(lapply(configuration[[table_name]][[level]], "[", 1)) == dataset_name &
                         unlist(lapply(configuration[[table_name]][[level]], "[", 2)) == var_name)[1]
    if (is.na(position)) {
      stop("No entry found that matches 'dataset_name' and 'var_name'.")
    }
  } else {
    if (is.null(position)) {
      stop("At least ('dataset_name', 'var_name') or 'position' must be specified.")
    }
    
    all_entries <- length(unlist(configuration[[table_name]], recursive = FALSE))
    if (position < 1 || position > all_entries) {
      stop("'position' must be in the range [1, # of table entries]")
    }

    found <- FALSE
    level <- 1
    index_of_first <- 1
    while (!found && level < 5) {
      if (position <= (index_of_first + length(configuration[[table_name]][[level]]) - 1)) {
        found <- TRUE
      } else {
        index_of_first <- index_of_first + length(configuration[[table_name]][[level]])
        level <- level + 1
      }
    }
    position <- position - index_of_first + 1
  }

  configuration[[table_name]][[level]][[position]] <- NULL

  configuration
}
