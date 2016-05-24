ConfigRemoveDefinition <- function(configuration, name) {
  configuration$definitions[[name]] <- NULL

  configuration
}
