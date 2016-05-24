ConfigEditDefinition <- function(configuration, name, value, confirm = TRUE) {
  continue <- TRUE
  if (name %in% names(configuration$definitions)) {
    if (confirm) {
      while (continue != 'y' && continue != 'n') {
        continue <- readline("WARNING: The definition already exists. It will be replaced. Continue? (y/n)\n")
      }
      continue <- ifelse(continue == 'y', TRUE, FALSE)
    }
  }
  if (continue) {
    configuration$definitions[[name]] <- value
  }
  
  configuration
}
