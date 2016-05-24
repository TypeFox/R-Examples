#' Generate command line arguments for features
#'
#' These functions are exported only for technical reasons
#' (because they are S3 methods) and are not intended for
#' users.
#'
#' @param feature The feature for which the argument is generated
#' @param model The complete model for which the argument is generated
#' @keywords internal
conv_to_ms_arg <- function(feature, model) UseMethod("conv_to_ms_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.default <- function(feature, model) {
  stop("Unknown feature", call. = FALSE)
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
ms_class <- R6Class("ms", inherit = simulator_class,
  private = list(
    name = "ms",
    priority = 100,
    binary = NULL
  ),
  public = list(
    initialize = function(priority = 300) {
      if (!requireNamespace("phyclust", quietly = TRUE)) {
        stop("Please install package 'phyclust' to use ms", call. = FALSE)
      }
      assert_that(is.numeric(priority) && length(priority) == 1)
      private$priority <- priority
    },
    create_cmd_tempalte = function(model) {
      cmd <- read_cache(model, "ms_cmd")
      if (is.null(cmd)) {
        cmd <- paste(vapply(model$features, conv_to_ms_arg,
                            FUN.VALUE = character(1), model),
                     collapse = "")
        cmd <- paste0("c('", cmd, "')")
        cache(model, "ms_cmd", cmd)
      }

      cmd
    },
    simulate = function(model, parameters=numeric()) {
      stopifnot(length(parameters) == 0 | all(is.numeric(parameters)))

      # Generate the simulation commands
      cmd_template <- self$create_cmd_tempalte(model)
      sample_size <- sum(get_sample_size(model, for_sim = TRUE))

      sim_cmds <- lapply(1:get_locus_group_number(model), function(group) {
        fill_cmd_template(cmd_template, model, parameters, group)
      })

      wd <- getwd()
      setwd(tempdir())

      # Do the actual simulation
      files <- lapply(sim_cmds, function(sim_cmd) {
        vapply(1:nrow(sim_cmd), function(j) {
          file <- tempfile("ms")
          ret <- phyclust::ms(sample_size,
                              sim_cmd[j, "locus_number"],
                              sim_cmd[j, "command"],
                              temp.file = file)

          if (!file.exists(file)) stop("ms simulation failed")
          file
        }, character(1))
      })

      setwd(wd)

      # Parse the output and calculate summary statistics
      if (requires_segsites(model) || requires_trees(model)) {
        output <- parse_ms_output(files, #nolint
                                  get_sample_size(model, for_sim = TRUE),
                                  get_locus_number(model))
      } else {
        output <- list(seg_sites = NULL, trees = NULL)
      }

      cmds <- lapply(sim_cmds, function(cmd) {
        paste("ms", sample_size, cmd[ , 1], cmd[ , 2])
      })

      sum_stats <- calc_sumstats_from_sim(output$segsites, output$trees,
                                          files, model, parameters, cmds, self)

      # Clean Up
      unlink(unlist(files))
      sum_stats
    },
    get_cmd = function(model) {
      template <- self$create_cmd_tempalte(model)
      cmd <- fill_cmd_template(template, model, NULL, 1, eval_pars = FALSE)
      paste("ms",
            sum(get_sample_size(model, TRUE)),
            cmd[1, "locus_number"],
            cmd[1, "command"])
    },
    get_info = function() c(name = "ms",
                            version = paste0("phyclust_",
                                             packageVersion("phyclust")))
  )
)

has_ms <- function() !is.null(simulators[["ms"]])


#' Simulator: ms
#'
#' This function adds the simulator 'ms' to the list of available simulators.
#' In order to use 'ms', you need to install the CRAN package \pkg{phyclust}.
#' By default, 'scrm' will be preferred over 'ms'. Raise the priority of 'ms'
#' to change this behavior.
#'
#' @references
#' Richard R. Hudson.
#' Generating samples under a Wright-Fisher neutral model of genetic variation.
#' Bioinformatics (2002) 18 (2): 337-338
#' doi:10.1093/bioinformatics/18.2.337
#'
#' @name simulator_ms
#' @param priority The priority for this simulator. If multiple simulators
#'   can simulate a model, the one with the highest priority will be used.
#' @export
#' @examples
#' # To prefer ms to scrm:
#' \dontrun{activate_ms(priority = 500)}
#' @family simulators
activate_ms <- function(priority = 300) {
  register_simulator(ms_class$new(priority))
  reset_cache()
  invisible(NULL)
}
