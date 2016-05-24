ConfigShowTable <- function(configuration, dataset_type, line_numbers = NULL) {
  table_name <- dataset_type
  header <- paste("| Matches in", gsub("_", " ", table_name), "|\n")
  cat(paste0(paste(rep("-", nchar(header) - 1), collapse = ''), "\n"))
  cat(header)
  cat(paste0(paste(rep("-", nchar(header) - 1), collapse = ''), "\n"))
  cat("#dataset_name, var_name[, main_path[, file_path[, nc_var_name[, suffix[, var_min[, var_max]]]]]]\n")

  if (is.null(line_numbers)) {
    line_numbers <- 1:length(unlist(configuration[[table_name]], recursive = FALSE))
  }
  line_number <- 1

  level <- 1
  invisible(lapply(configuration[[table_name]], 
    function(x) {
      cat(paste("# Level", level, "#\n"))
      lapply(x, 
        function(y) {
          cat(paste(line_numbers[line_number], ": ", paste(unlist(y), collapse = ', '), "\n", sep = ''))
          line_number <<- line_number + 1
        }
      )
      level <<- level + 1 
    }
  ))
}
