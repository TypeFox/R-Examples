ConfigApplyMatchingEntries <- function(configuration, var, exp = NULL, obs = NULL, show_entries = FALSE, show_result = TRUE) {
  ## Function to tell if a regexpr() match is a complete match to a specified name
  isFullMatch <- function(x, name) {
    ifelse(x > 0 && attributes(x)$match.length == nchar(name), TRUE, FALSE)
  }

  var_entries_in_exps <- c()
  if (length(unlist(configuration$experiments, recursive = FALSE)) > 0) {
    var_entries_in_exps <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$experiments, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }
  var_entries_in_obs <- c()
  if (length(unlist(configuration$observations, recursive = FALSE)) > 0) {
    var_entries_in_obs <- which(unlist(lapply(lapply(as.list(unlist(lapply(configuration$observations, lapply, "[[", 2))), regexpr, var), isFullMatch, var) > 0))
  }

  exp_info <- list()
  jmod <- 1
  for (mod in exp) {
    mod_var_matching_entries <- mod_var_matching_indices <- mod_var_matching_entries_levels <- c()
    
    if (length(unlist(configuration$experiments, recursive = FALSE)) > 0) {
      mod_entries_in_exps <- which(unlist(lapply(lapply(unlist(lapply(configuration$experiments, lapply, "[[", 1), recursive = FALSE), regexpr, mod), isFullMatch, mod)))
      if (length(mod_entries_in_exps) > 0) {
        mod_var_matching_indices <- intersect(var_entries_in_exps, mod_entries_in_exps)
        mod_var_matching_entries <- unlist(configuration$experiments, recursive = FALSE)[mod_var_matching_indices]
        exps_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$experiments[[x]]))})
        mod_var_matching_entries_levels <- unlist(exps_levels)[intersect(var_entries_in_exps, mod_entries_in_exps)]
      }
    }

    if (length(mod_var_matching_entries) == 0) {
      stop(paste('Error: There are no matching entries in the configuration file for the experiment', mod, 'and the variable', var, 
                 '. Please check the configuration file.)'))
    } else {
      if (show_entries) {
        header <- paste0("# Matching entries for experiment '", exp[jmod], "' and variable '", var, "' #\n")
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(header)
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        ConfigShowTable(list(experiments = list(mod_var_matching_entries)), 'experiments', mod_var_matching_indices)
        cat("\n")
      }
      result <- .ConfigGetDatasetInfo(mod_var_matching_entries, 'experiments')
      if (show_result) {
        cat(paste0("The result of applying the matching entries to experiment name '", exp[jmod], "' and variable name '", var, "' is:\n"))
        configuration$definitions[["VAR_NAME"]] <- var
        configuration$definitions[["EXP_NAME"]] <- exp[jmod]
        fields <- c("MAIN_PATH: ", "FILE_PATH: ", "NC_VAR_NAME: ", "SUFFIX: ", "VAR_MIN: ", "VAR_MAX: ")
        values <- lapply(result, lapply, function (x) .ConfigReplaceVariablesInString(x, configuration$definitions, TRUE))
        lapply(paste0(fields, unlist(values), "\n"), cat)
        cat("\n")
      }
      exp_info <- c(exp_info, list(result))
    }

    jmod <- jmod + 1
  }

  obs_info <- list()
  jobs <- 1
  for (ref in obs) {
    ref_var_matching_entries <- ref_var_matching_indices <- ref_var_matching_entries_levels <- c()
    
    if (length(unlist(configuration$observations, recursive = FALSE)) > 0) {
      ref_entries_in_obs <- which(unlist(lapply(lapply(unlist(lapply(configuration$observations, lapply, "[[", 1), recursive = FALSE), regexpr, ref), isFullMatch, ref)))
      if (length(ref_entries_in_obs) > 0) {
        ref_var_matching_indices <- intersect(var_entries_in_obs, ref_entries_in_obs)
        ref_var_matching_entries <- unlist(configuration$observations, recursive = FALSE)[ref_var_matching_indices]
        obs_levels <- lapply(as.list(1:4), f <- function(x) {x <- array(x, length(configuration$observations[[x]]))})
        ref_var_matching_entries_levels <- unlist(obs_levels)[intersect(var_entries_in_obs, ref_entries_in_obs)]
      }
    }

    if (length(ref_var_matching_entries) == 0) {
      stop(paste('Error: There are no matching entries in the configuration file for the observation', ref, 'and the variable', var, 
                 '. Please check the configuration file.)'))
    } else {
      if (show_entries) {
        header <- paste0("# Matching entries for observation '", obs[jobs], "' and variable '", var, "' #\n")
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        cat(header)
        cat(paste0(paste(rep("#", nchar(header) - 1), collapse = ''), "\n"))
        ConfigShowTable(list(observations = list(ref_var_matching_entries)), 'observations', ref_var_matching_indices)
        cat("\n")
      }
      result <- .ConfigGetDatasetInfo(ref_var_matching_entries, 'observations')
      if (show_result) {
        cat(paste0("The result of applying the matching entries to observation name '", obs[jobs], "' and variable name '", var, "' is:\n"))
        configuration$definitions[['VAR_NAME']] <- var
        configuration$definitions[["OBS_NAME"]] <- obs[jobs]
        fields <- c("MAIN_PATH: ", "FILE_PATH: ", "NC_VAR_NAME: ", "SUFFIX: ", "VAR_MIN: ", "VAR_MAX: ")
        values <- lapply(result, lapply, function (x) .ConfigReplaceVariablesInString(x, configuration$definitions, TRUE))
        lapply(paste0(fields, unlist(values), "\n"), cat)
        cat("\n")
      }
      obs_info <- c(obs_info, list(result))
    }

    jobs <- jobs + 1
  }

  invisible(list(exp_info = exp_info, obs_info = obs_info))
}
