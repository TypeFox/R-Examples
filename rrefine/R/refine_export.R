#' export data from OpenRefine
#'
#' @param project.name name of project to be exported
#' @param project.id unique identifier for project to be exported
#' @param format file format of project to be exported, default is 'csv'
#' @param col.names logical indicator for whether column names should be included
#' @param encoding character encoding for exported data, default is 'UTF-8'
#' @export
#' @examples
#' \dontrun{
#' refine_export("purple_rain")
#' refine_export(project.id = 1901017388690, format = "csv")
#' refine_export(project.name = "Untitled", project.id = 1901017388888)
#' }
#'

refine_export <- function(project.name = NULL, project.id = NULL, format = "csv", col.names = TRUE, encoding = "UTF-8") {

    refine_check()

    project.id <- refine_id(project.name, project.id)

    cont <- httr::content(
        httr::POST(
            paste0(refine_path(), "/",
                   "command/core/export-rows/",
                   project.id, ".", format),
            body = c(engine = list(facets = "", mode="row-based"),
                     project = project.id, format = format),
            encode = "form"),
        type = "text/csv",
        as = "text",
        encoding = encoding)

    readr::read_csv(cont, col_names = col.names, col_types = NULL)

}
