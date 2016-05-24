stat_file_class <- R6Class("stat_file", inherit = sumstat_class, #nolint
  private = list(
    folder = NULL,
    req_files = TRUE
  ),
  public = list(
    initialize = function(folder) {
      dir.create(folder, showWarnings = FALSE)
      private$folder <- folder
      super$initialize("file", identity)
    },
    calculate = function(seg_sites, trees, files, model) {
      if (is.list(files)) files <- unlist(files)
      if (!all(file.copy(files, private$folder, overwrite = FALSE)))
        stop("Failed to copy simulated files. Look for warnings.")

      file.path(private$folder, basename(files))
    }
  )
)


#' Summary Statistic: Files
#'
#' This "summmary statistic" returns files with the raw results of
#' the simulation. Multiple files are returned in case coala needs
#' multiple calls to simulators to simulate the model. These files
#' do not contain any post processing of the results done by coala,
#' e.g. \code{\link{feat_unphased}} and
#' \code{\link{feat_ignore_singletons}}.
#'
#' @param folder The path to a folder. The files will be created there.
#' @return A character vector containing the files in order in which they
#'         where created.
#' @export
#' @template summary_statistics
#' @examples
#' folder <- tempfile("coala-test")
#' model <- coal_model(10, 1) +
#'   feat_mutation(5) +
#'   sumstat_file(folder)
#' simulate(model)$file
#'
#' unlink(folder, recursive = TRUE)  # Clean up
sumstat_file <- function(folder) {
  stat_file_class$new(folder)
}
