conv_to_scrm_arg <- function(feature, model) UseMethod("conv_to_scrm_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.default <- function(feature, model) {
  stop("Unknown feature when generating scrm command")
}


scrm_create_cmd_template <- function(model) {
  cmd <- read_cache(model, "scrm_cmd")
  if (is.null(cmd)) {
    cmd <- paste(vapply(model$features, conv_to_scrm_arg,
                        FUN.VALUE = character(1), model),
                 collapse = "")
    cmd <- paste0("c('", cmd, "')")
    cache(model, "scrm_cmd", cmd)
  }
  cmd
}


#' @importFrom scrm scrm
#' @include simulator_class.R
#' @include simulator_ms.R
scrm_class <- R6Class('Scrm', inherit = simulator_class, #nolint
  private = list(
    name = "scrm",
    version = paste0("scrm_", packageDescription("scrm", fields = "Version"))
  ),
  public = list(
    simulate = function(model, parameters) {
      cmd_template <- scrm_create_cmd_template(model)
      sample_size <- sum(get_sample_size(model, for_sim = TRUE))

      sim_cmds <- lapply(1:get_locus_group_number(model), function(group) {
        fill_cmd_template(cmd_template, model, parameters, group)
      })

      locus_number <- sum(get_locus_number(model))
      seg_sites <- list()
      if (requires_segsites(model)) length(seg_sites) <- locus_number

      trees <- list()
      if (requires_trees(model)) length(trees) <- locus_number

      sim_number <- sum(sapply(sim_cmds, nrow))
      if (requires_files(model)) {
        files <- sapply(1:sim_number, function(x) tempfile("scrm"))
      } else {
        files <- rep("", sim_number)
      }

      cl <- 1
      j <- 0
      for (cmds in sim_cmds) {
        for (i in 1:nrow(cmds)) {
          j <- j + 1
          stats <- scrm(paste(sample_size, cmds[i, 1], cmds[i, 2]), files[j])

          if (requires_segsites(model)) {
            seg_sites[cl:(cl + cmds[i, 1] - 1)] <- stats$seg_sites[] #nolint
          }

          if (requires_trees(model)) {
            split_tree <- lapply(stats$trees, function(x) {
              strsplit(x, split = "\n", fixed = TRUE)[[1]]
            })
            trees[cl:(cl + cmds[i, 1] - 1)] <- split_tree[]  #nolint
          }
          cl <- cl + cmds[i, 1]
        }
      }

      seg_sites <- lapply(seg_sites, function(x) {
        if (!is_segsites(x)) {
          x <- create_segsites(x, as.numeric(colnames(x)), check = FALSE)
        }
        x
      })

      cmds <- lapply(sim_cmds, function(cmd) {
        paste("scrm", sample_size, cmd[ , 1], cmd[ , 2])
      })

      stats <- calc_sumstats_from_sim(seg_sites, trees, files,
                                      model, parameters, cmds, self)
      unlink(files)

      stats
    },
    get_cmd = function(model) {
      template <- scrm_create_cmd_template(model)
      cmd <- fill_cmd_template(template, model, NULL, 1, eval_pars = FALSE)
      paste("scrm",
            sum(get_sample_size(model, TRUE)), cmd[1, 1], cmd[1, 2])
    },
    get_info = function() c(name = "scrm", version = private$version)
  )
)

#' Simulator: scrm
#'
#' This function adds the simulator 'scrm' to the list of available simulators.
#' It is provided via the CRAN package \pkg{scrm} and should be always installed
#' alongside with \pkg{coala}. It should be activated automatically, and this
#' function is only needed to change it \code{priority}.
#'
#' @references
#' Paul R. Staab, Sha Zhu, Dirk Metzler and Gerton Lunter (2015).
#' "scrm: efficiently simulating long sequences using the approximated
#' coalescent with recombination."
#' Bioinformatics, 31(10), pp. 1680-1682.
#' http://dx.doi.org/10.1093/bioinformatics/btu861
#'
#' @name simulator_scrm
#' @inheritParams simulator_ms
#' @family simulators
#' @export
#' @examples
#' # Change scrm's priority
#' model <- coal_model(10, 1) + feat_mutation(5)
#' model # scrm is used by default
#' activate_scrm(250)
#' model # Now ms is used instead (if installed)
#' activate_scrm(550)
#' model # Now scrm is used again
activate_scrm <- function(priority = 400) {
  register_simulator(scrm_class$new(priority))
  reset_cache()
  invisible(NULL)
}
