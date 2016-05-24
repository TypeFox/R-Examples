#' delete project from OpenRefine
#'
#' @param project.name name of project to be deleted
#' @param project.id unique identifier for open refine project to be deleted
#' @export
#' @examples
#' \dontrun{
#' refine_delete(project.name = "foo")
#' refine_delete(project.id = 1901017388690)
#' refine_delete(project.name = "Untitled", project.id = 1901018888332)
#' }
#'
refine_delete <- function(project.name = NULL, project.id = NULL) {

    refine_check()

    project.id <- refine_id(project.name, project.id)

    x <- readline(prompt = "are you sure you want to delete this project? (Y/N): ")

    if(x == "Y") {
        httr::POST(
            paste0(refine_path(), "/", "command/core/delete-project"),
            body = list(project = project.id),
            encode = "form")
        message("project deleted")
    } else
        stop("aborting delete process")

}
