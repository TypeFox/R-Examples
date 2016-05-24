#' @title Shiny Demonstrations
#' @description Demonstrate package functionality via Shiny apps
#' @param ex Example to run.
#' @details Demonstrations available:
#'   \code{"gvf"} Gradually-varied flow.
#'
#' @examples
#' \dontrun{
#' # get list of available demos
#' demo_shiny()
#' # run the gradually-varied flow demo
#' demo_shiny("gvf")
#' }
#'
#' @export
demo_shiny = function(ex){
  # locate all the shiny app examples that exist
  app.dir = system.file("shiny-examples", package = "rivr")
  valid.examples = list.files(app.dir)
  if(missing(ex))
    return(
      paste0("Available shiny examples: '", 
             paste(valid.examples, collapse = "', '"), "'")
    )
  if(!(ex %in% valid.examples))
    stop("example '", ex, "' not found")  
  if(!requireNamespace("shiny", quietly = TRUE))
    stop("package 'shiny' is not installed.")
  shiny::runApp(file.path(app.dir, ex), display.mode = "normal")
  invisible(NULL)
}
