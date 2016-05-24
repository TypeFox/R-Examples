print_model_part <- function(model, part, part_name) {
  part_list <- model[[part]]
  cat(part_name, ": ", sep = "")
  if (length(part_list) == 0) cat("None\n")
  else if (length(part_list) == 1) part_list[[1]]$print()
  else {
    cat("\n")
    for (x in part_list) {
      cat("* ")
      x$print()
    }
  }
}

#' @export
print.coalmodel <- function(x, ...) {
  print_model_part(x, "features", "Features")
  cat("\n")

  print_model_part(x, "parameter", "Parameter")
  cat("\n")

  print_model_part(x, "loci", "Loci")
  cat("\n")

  print_model_part(x, "sum_stats", "Summary Statistics")
  cat("\n")

  suppressWarnings(simprog <- select_simprog(x))
  cat("Simulator: ")
  if (is.null(simprog)) {
    cat("No suitable program available\n")
  } else {
    cat(simprog$get_name(), "\n")
    cmd <- try(cmd <- get_cmd(x), silent = TRUE)
    if ("try-error" %in% class(cmd)) {
      cat("Failed to generate command: \n  ")
      cat(cmd)
    } else {
      cat("Command:", paste(cmd, collapse = "\n"), "\n")
    }
  }
}
