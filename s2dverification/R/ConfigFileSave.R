ConfigFileSave <- function(configuration, file_path, confirm = TRUE) {
  continue <- TRUE
  if (file.exists(file_path)) {
    if (confirm) {
      while (continue != 'y' && continue != 'n') {
        continue <- readline(paste0("WARNING: The configuration file '", file_path, "' already exists. It will be replaced. Continue? (y/n)\n"))
      }
      continue <- ifelse(continue == 'y', TRUE, FALSE)
    }
  }
  if (continue) {
    file_conn <- file(file_path)
    file_text <- c(
"# s2dverification configuration file",
"#",
"# Check ?ConfigFileOpen after loading s2dverification for detailed ",
"# documentation on this configuration file.",
""
                  )

    file_text <- c(file_text,
      paste(rep("#", nchar("definitions") + 2), collapse = ''),
      paste0("!!definitions"),
      paste(rep("#", nchar("definitions") + 2), collapse = '')
                  )
    defaults <- configuration$definitions[grep("^DEFAULT_", names(configuration$definitions))]
    definitions <- configuration$definitions[-grep("^DEFAULT_", names(configuration$definitions))]
    file_text <- c(file_text, as.vector(paste(names(defaults), unlist(defaults), sep = " = ")))
    file_text <- c(file_text, as.vector(paste(names(definitions), unlist(definitions), sep = " = ")))
    file_text <- c(file_text, "")

    table_names <- c("experiments", "observations")
    for (table_name in table_names) {
      if (table_name == "experiments") {
        dataset_type <- 'exp'
      } else {
        dataset_type <- 'obs'
      }
      file_text <- c(file_text,
"",
  paste(rep("#", nchar(table_name) + 11), collapse = ''),
  paste0("!!table of ", gsub("_", " ", table_name)),
  paste(rep("#", nchar(table_name) + 11), collapse = ''),
  paste0("#", dataset_type, "_name, var_name[, ", dataset_type, "_main_path[, ", dataset_type, "_file_path[, nc_var_name[, suffix[, var_min[, var_max]]]]]]")
                  )
      # Some heavy entry processing still to do here, to put asterisks, empty spaces, double quotes, and reduce options
      file_text <- c(file_text, unlist(lapply(configuration[[table_name]], function (x) lapply(x, function (y) paste(unlist(y), collapse = ", ")))))
    }
  
    writeLines(file_text, file_conn)
    close(file_conn)
  }

  invisible(continue)
}
