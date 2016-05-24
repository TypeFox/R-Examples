#' Check which simulator can simulate a model
#'
#' This function checks which of the available simulators can
#' simulate a given model. It also states the problems for the ones that
#' are incompatible with the model.
#'
#' @param model The model which is checked
#' @export
#' @seealso Do view the priority of the simulators: \code{\link{list_simulators}}
#' @examples
#' model <- coal_model(10, 1) +
#'   feat_mutation(5, fixed = TRUE)
#' check_model(model)
check_model <- function(model) {
  force(model)
  for (simprog_name in ls(simulators)) {
    cat(simprog_name, ": ")
    simprog <- get_simulator(simprog_name)

    cmd <- try(simprog$get_cmd(model), silent = TRUE)

    if (inherits(cmd, "try-error")) {
      cat(cmd, "\n")
    } else {
      cat("OK\n\n")
    }
  }
}
