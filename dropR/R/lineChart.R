#' Show Line Chart Output on a HTML canvas
#' 
#' create a HTML SVG Line Chart Output using NVD3.js
#' This function should be called from ui.R
#' in a shiny web application. 
#' 
#' @param inputId input identifier for the output function, i.e. name of the list element in shiny.
#' @param width defaults to 100\%. 
#' @param height defaults to 400px.
#' @export
lineChartOutput <- function(inputId, width="100%", height="400px") {
  style <- sprintf("width: %s; height: %s;",
                   shiny::validateCssUnit(width), shiny::validateCssUnit(height))
  
  tags = NULL
  
  shiny::tagList(
    # Include CSS/JS dependencies. Use "singleton" to make sure that even
    # if multiple lineChartOutputs are used in the same page, we'll still
    # only include these chunks once.
#     shiny::singleton(tags$head(
#       shiny::tags$script(src="d3/d3.v3.min.js"),
#       shiny::tags$script(src="nvd3/nv.d3.min.js"),
#       shiny::tags$link(rel="stylesheet", type="text/css", href="nvd3/nv.d3.min.css"),
#       shiny::tags$script(src="linechart-binding.js")
#     )
#     ),
    
    shiny::tags$head(
      shiny::tags$script(src="d3/d3.v3.min.js"),
      shiny::tags$script(src="nvd3/nv.d3.min.js"),
      shiny::tags$link(rel="stylesheet", type="text/css", href="nvd3/nv.d3.min.css"),
      shiny::tags$script(src="linechart-binding.js")
      ),
    
    
    shiny::div(id=inputId, class="nvd3-linechart", style=style,
               shiny::tag("svg", list()))
  )
}

