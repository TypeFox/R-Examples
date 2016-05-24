#' helper function to configure and call path to OpenRefine
#'
#' @param host host path for your OpenRefine application
#' @param port port number for your OpenRefine application
#' @return path to be executed

refine_path <- function(host = "127.0.0.1", port ="3333") {

    paste0("http://", host, ":", port)

}

#' get all project metadata from OpenRefine
#'
#' @examples
#' \dontrun{
#' refine_metadata()
#' }
#'

refine_metadata <- function() {

    httr::content(httr::GET(paste0(refine_path(),
        "/",
        "command/core/get-all-project-metadata"))
    )

}

#' helper function to get OpenRefine project.id by project.name
#'
#' @param project.name name of project to be exported or deleted
#' @param project.id unique identifier for project to be exported or deleted
#' @return unique id of project to be exported or deleted

refine_id <- function(project.name, project.id) {

    if (is.null(project.name) & is.null(project.id)){

        stop("you must supply either a project name or project id")

    }

    if (!is.null(project.id)) {

        project.id <- project.id

    } else {

        resp <- refine_metadata()

        name <- NULL

        id <- names(rlist::list.mapv(
        rlist::list.filter(resp[["projects"]],
            name == project.name),
        name))

        if (length(id) == 1) {

            project.id <- id

        } else if (length(id) == 0) {

            stop(paste0("there are no projects found named '",
            project.name,
            "' ... try adjusting the project.name argument or use project.id"))

        } else {

            stop(paste0("there are mulitple projects named '",
                project.name,
                "' ... try addding a project.id"))

        }
    }
}

#' helper function to check if rrefine can connect to OpenRefine
#'
#' @return error message if rrefine is unable to connect to OpenRefine, otherwise is invisible

refine_check <- function() {

    tryCatch(
        expr = httr::GET(refine_path()),
        error = function(e)
            cat("rrefine is unable to connect to OpenRefine ... make sure OpenRefine is running"))

    invisible()
}
