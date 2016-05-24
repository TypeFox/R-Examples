ConfigAddEntry <- function(configuration, dataset_type, position = 'last', dataset_name = ".*", var_name = ".*", main_path = "*", file_path = "*", nc_var_name = "*", suffix = "*", varmin = "*", varmax = "*") {
  table_name <- dataset_type
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

  index_of_first <- 0
  index_of_last <- 0
  for (i in 1:level) {
    index_of_first <- index_of_first + ifelse(i == 1, 1, length(configuration[[table_name]][[i - 1]]))
    index_of_last <- index_of_last + length(configuration[[table_name]][[i]])
  }

  if (position == 'last') {
    position <- index_of_last - index_of_first + 1 + 1
  } else if (position == 'first') {
    position <- 1
  } else {
    if (position < index_of_first || position > index_of_last + 1) {
      stop("'position' must be in the range [index of first table entry in corresponding level, index of last table entry in corresponding level + 1]")
    }
    position <- position - index_of_first + 1
  }

  if (dataset_type == 'experiments' || dataset_type == 'observations') {
    configuration[[table_name]][[level]] <- append(configuration[[table_name]][[level]], list(c(dataset_name, var_name, main_path, file_path, nc_var_name, suffix, varmin, varmax)), after = position - 1)
  } else {
    stop("'dataset_type' must be one of 'experiments' or 'observations'")
  }

  configuration
}
