dataTableDependency <- function() {
  list(
    htmlDependency(
      "datatables", "1.10.2",
      c(file = system.file("www/datatables", package = "shinybootstrap2")),
      script = "js/jquery.dataTables.min.js"
    ),
    htmlDependency(
      "datatables-bootstrap", "1.10.2",
      c(file = system.file("www/datatables", package = "shinybootstrap2")),
      stylesheet = c("css/dataTables.bootstrap.css", "css/dataTables.extra.css"),
      script = "js/dataTables.bootstrap.js"
    )
  )
}

#' Create a table output element
#'
#' Render a \code{\link{renderDataTable}} within an application page.
#' \code{renderDataTable} uses the DataTables Javascript library to create an
#' interactive table with more features.
#'
#' @param outputId output variable to read the table from
#' @return A table output element that can be included in a panel
#'
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   library(shiny)
#'   shinyApp(
#'     ui = bootstrapPage(
#'       dataTableOutput('table')
#'     ),
#'     server = function(input, output) {
#'       output$table <- renderDataTable(iris)
#'     }
#'   )
#' }
#' @export
dataTableOutput <- function(outputId) {
  attachDependencies(
    div(id = outputId, class="shiny-datatable-output"),
    dataTableDependency()
  )
}
