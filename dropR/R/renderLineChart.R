#' Render a line Chart in shiny 
#' 
#' This function renders a line chart and should be called from server.R in a shiny application.
#'  
#' @param expr an expression to rendered
#' @param env environment, defaults to parent.frame().
#' @param quoted logical, defaults to FALSE.
#' @export
renderLineChart <- function(expr, env=parent.frame(), quoted=FALSE) {
  # This piece of boilerplate converts the expression `expr` into a
  # function called `func`. It's needed for the RStudio IDE's built-in
  # debugger to work properly on the expression.
  # get rid of the NOTE for CRAN compliance, see here:
  # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  func = NULL
  shiny::installExprFunction(expr, "func", env, quoted)
  
  function() {
    dataframe <- func()
    
    mapply(function(col, name) {
      
      values <- mapply(function(val, i) {
        list(x = i, y = val)
      }, col, 1:nrow(dataframe), SIMPLIFY=FALSE, USE.NAMES=FALSE)
      
      list(key = name, values = values)
      
    }, dataframe, names(dataframe), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }
}

