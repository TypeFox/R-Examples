#' Search the working directory and the run path for an executable
#'
#' This function tries to find a binary of a given name by looking
#' in the current working directory and the directories listed in
#' the $PATH environment variable. If an environment variable with
#' name equal to to binaries name in upper case is given, it also
#' tries to use this as path of the binary.
#'
#' @param name The name of the executable to look for
#' @param envir_var the name of the environment variable to use
#' @return The complete path of the executable is found, or "NULL" if not.
#' @keywords internal
search_executable <- function(name, envir_var = NULL) {
  # See if an environment variable is given
  exe <- NULL

  if (!is.null(envir_var)) {
    exe_path <- Sys.getenv(envir_var)
    if (exe_path != "") {
      if (!file.exists(exe_path)) {
        warning("Can not find ", name, " executable at '", exe_path, "'")
      } else {
        return(exe_path)
      }
    }
  }

  # Try to find it in the PATH folders and the Working directory
  if (Sys.info()[["sysname"]] == "Windows") {
    run_path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
  } else {
    run_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  }

  candidates <- do.call(c, lapply(name, function(x) {
    file.path(c(getwd(), run_path), x)
  }))

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      exe <- candidate
      break
    }
  }

  exe
}
