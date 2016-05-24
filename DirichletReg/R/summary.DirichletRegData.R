summary.DirichletRegData <- function(object, ...){

  if(interactive()) writeLines("")
  writeLines(paste0("This object contains compositional data with ", attr(object, "dims"), " dimensions."))

  writeLines(paste0("Number of observations: ", attr(object, "obs"),
    " of which ", attr(object, "valid_obs"), " ( ", round(100*attr(object, "valid_obs")/attr(object, "obs"),2), "% ) are valid."))

  if(attr(object, "normalized") || attr(object, "transformed")){
    cat("\nNote: The data were ")
    if(attr(object, "normalized")) cat("normalized")
    if(attr(object, "normalized") && attr(object, "transformed")) cat(" and ")
    if(attr(object, "transformed")) cat("transformed")
    cat(".\n")
  }

  if(interactive()) writeLines("")

}
