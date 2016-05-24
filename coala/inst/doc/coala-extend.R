## ------------------------------------------------------------------------
library(R6)
library(coala)
stat_segsites_class <- R6Class("stat_segsites", inherit = sumstat_class,
  private = list(req_segsites = TRUE),
  public = list(
    calculate = function(segsites, trees, files, model) segsites
  )
)

sumstat_seg_sites <- function(name = "seg_sites", transformation = identity) {
  stat_segsites_class$new(name, transformation)
}

## ------------------------------------------------------------------------
stat_file_class <- R6Class("stat_file", inherit = sumstat_class,
  private = list(folder = NULL, req_files = TRUE),
  public = list(
    initialize = function(folder) {
      dir.create(folder, showWarnings = FALSE)
      private$folder <- folder
      super$initialize("file", identity)
    },
    calculate = function(seg_sites, trees, files, model) {
      file.copy(files, private$folder, overwrite = FALSE)
      file.path(private$folder, basename(files))
    }
  )
)

