ConfigFileOpen <- function(file_path, silent = FALSE, stop = FALSE) {
  if (!silent) {
    cat(paste("* Reading configuration file:", file_path, "\n"))
  }
  # Read the data from the configuration file.
  ## Remove comments, tabulations, spaces, empty lines, ...
  all_lines <- readLines(file_path)
  all_lines <- gsub("\t", "", all_lines)
  all_lines <- gsub(" ", "", all_lines)
  all_lines <- all_lines[-grep("^#", all_lines)]
  all_lines <- all_lines[-grep("^$", all_lines)]
  ## Detect key lines
  key_positions <- grep("^!!", all_lines)

  ## Check that the format of the configuration file is right.
  if (length(key_positions) != 3) {
    stop('Error: The configuration file is corrupted or outdated: the key lines do not match the expected pattern.')
  }

  ## Start parsing the configuration.
  # The variables that are used in the configuration filed are kept in 
  # 'definitions', an associative array (key-value array or dictionary).
  definitions <- list()
  ## Parse the variables definitions in the whole configuration file
  if (key_positions[1] + 1 < key_positions[2]) {
    all_definitions <- all_lines[(key_positions[1] + 1):(key_positions[2] - 1)]
  } else {
    all_definitions <- c()
  }
  if (length(grep("=", all_definitions)) == length(all_definitions)) {
    for (definition in all_definitions) {
      if (length(which(strsplit(definition, "")[[1]] == "=")) == 1) {
        var_name <- strsplit(definition, "=")[[1]][1]
        tmp_value <- strsplit(definition, "=")[[1]][2]
        var_value <- ifelse(is.na(tmp_value), "", tmp_value)
        if ((length(which(strsplit(var_value, "")[[1]] == "$")) %% 2) == 0) {
          definitions[[var_name]] <- var_value
        } else {
          stop('Error: The configuration file is corrupted: there are incorrect variable definition lines in the definition zone. A closing "$" symbol may be missing.')
        }
      } else {
        stop('Error: The configuration file is corrupted: there are incorrect definition lines in the definition zone.')
      }
    }
  } else {
    stop('Error: The configuration file is corrupted: there are malformed definition lines in the definition zone.')
  }
  mandatory_definitions <- c("DEFAULT_EXP_MAIN_PATH", "DEFAULT_EXP_FILE_PATH", 
                             "DEFAULT_NC_VAR_NAME", "DEFAULT_SUFFIX", "DEFAULT_VAR_MIN", 
                             "DEFAULT_VAR_MAX", "DEFAULT_OBS_MAIN_PATH", 
                             "DEFAULT_OBS_FILE_PATH", "DEFAULT_DIM_NAME_LONGITUDES",
                             "DEFAULT_DIM_NAME_LATITUDES", "DEFAULT_DIM_NAME_MEMBERS")
  if (any(!(mandatory_definitions %in% names(definitions)))) {
    cat("* WARNING: Some of the mandatory variables below are not defined in the configuration file. You can add them with ConfigFileOpen(), ConfigEditDefinition() and ConfigFileSave() or by editing the configuration file by hand, as specified in ?ConfigFileOpen.")
    if (stop) {
      stop(paste(mandatory_definitions, collapse = ', '))
    } else {
      cat(paste(mandatory_definitions, collapse = ', '))
    }
  }

  # Parse the entries in the tables
  ## These are the indices of the key positions in the vector of key positions
  tables_key_positions <- c(2, 3)
  current_table <- 1
  for (table_key_position in tables_key_positions) {
    datasets <- list(c(), c(), c(), c())

    if (table_key_position == 2) {
      id <- 'EXP'
    } else {
      id <- 'OBS'
    }
    default_values <- c(paste0("$DEFAULT_", id, "_MAIN_PATH$"), paste0("$DEFAULT_", id, "_FILE_PATH$"), "$DEFAULT_NC_VAR_NAME$", '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
    previous_values <- c(".*", ".*", default_values)
    table_lines <- c()
    table_end <- ifelse(table_key_position == max(tables_key_positions), length(all_lines), key_positions[table_key_position + 1] - 1)
    if ((key_positions[table_key_position] + 1) <= table_end) {
      table_lines <- all_lines[(key_positions[table_key_position] + 1):table_end]
      table_lines <- strsplit(table_lines, ",")
    }

    current_line <- 1
    for (entry in table_lines) {
      if (entry[1] == '"') {
        entry[1] <- previous_values[1]
      }
      if ((length(entry) > 1)) {
        if (entry[2] == '"') {
          entry[2] <- previous_values[2]
        }
      } else {
        stop('Error: The variable column must be defined in all the entries in the tables in the configuration file.')
      }
      if (length(entry) > length(default_values) + 2) {
        stop(paste0("Error: More elements than expected in the entry ", current_line, " in the configuration file."))
      }
      for (value_position in 1:length(default_values)) {
        if ((length(entry) > value_position + 1)) {
          if (entry[value_position + 2] == '"') {
            entry[value_position + 2] <- previous_values[value_position + 2]
          }
        } else {
          entry[value_position + 2] <- '*'
        }
      }
      if (entry[1] == '.*') {
        if (entry[2] == '.*') {
          datasets[[1]] <- c(datasets[[1]], list(entry))
        } else {
          datasets[[3]] <- c(datasets[[3]], list(entry))
        }
      } else {
        if (entry[2] == '.*') {
          datasets[[2]] <- c(datasets[[2]], list(entry))
        } else {
          datasets[[4]] <- c(datasets[[4]], list(entry))
        }
      }
      current_line <- current_line + 1
      previous_values <- entry
    }

    if (current_table == 1) {
      exps <- datasets
    } else if (current_table == 2) {
      obs <- datasets
    }

    current_table <- current_table + 1
  }
  
  if (!silent) {
    cat("* Config file read successfully.\n")
  }

  invisible(list(definitions = definitions, 
                 experiments = exps, 
                 observations = obs))
}
