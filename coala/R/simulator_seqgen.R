sg_mutation_models <- c("HKY", "GTR")


generate_tree_model <- function(model) {
  tree_model <- read_cache(model, "tree_model")

  if (is.null(tree_model)) {
    tree_model <- model

    # Features
    tree_model_features <- !vapply(model$features, function(x) {
      any(c("seg_sites_feat",
            "mutation",
            "outgroup") %in% class(x))
    }, logical(1)) #nolint
    if (all(tree_model_features)) stop("seq-gen not required")
    tree_model$features <- model$features[tree_model_features]

    # Summary Stastics
    tree_model$sum_stats <- create_sumstat_container()
    tree_model <- tree_model + sumstat_sg_trees()

    cache(model, "tree_model", tree_model)
  }

  tree_model
}


conv_to_seqgen_arg <- function(feature, model) UseMethod("conv_to_seqgen_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.default <- function(feature, model) {
  stop("Unknown feature when generating seqgen command")
}


sg_generate_opts <- function(model, parameters, locus,
                             seeds, for_cmd = FALSE) {
  locus_lengths <- get_locus_length(model, group = locus, total = FALSE)

  if (length(locus_lengths) == 5) {
    locus_lengths <- locus_lengths[c(1, 3, 5)]
  }

  cmd <- sg_generate_opt_cmd(model)
  #print(cmd)

  # Fill the parameters in the template
  sapply(seq(along = locus_lengths), function(i) {
    par_envir <- create_par_env(model, parameters, locus = locus,
                                locus_length = locus_lengths[i],
                                seed = seeds[i], for_cmd = for_cmd)
    paste(eval(parse(text = cmd[[i]]), envir = par_envir), collapse = " ")
  })
}


sg_generate_opt_cmd <- function(model) {
  cmd <- read_cache(model, "seqgen_cmd")

  if (is.null(cmd)) {
    if (has_trios(model)) is_outer <- c(TRUE, FALSE, TRUE)
    else is_outer <- FALSE

    cmd <- lapply(is_outer, function(outer) {
      cmd <- paste(vapply(model$features, conv_to_seqgen_arg,
                          FUN.VALUE = character(1), model = model),
                   collapse = "")
      cmd <- paste0("c('", cmd, "')")
    })

    cache(model, "seqgen_cmd", cmd)
  }
  cmd
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
seqgen_class <- R6Class("seqgen", inherit = simulator_class,
  private = list(
    name = "seqgen",
    binary = NULL,
    priority = 100
  ),
  public = list(
    initialize = function(binary = NULL, priority = 100) {
      # Try to automatically find a jar file and java if not given
      if (is.null(binary)) {
        binary <- search_executable(c("seqgen", "seq-gen",
                                      "seqgen.exe", "seq-gen.exe"), "SEQGEN")
      }
      if (is.null(binary)) stop("No binary file for seqgen found.")
      if (!file.exists(binary)) stop("seqgen binary (", binary,
                                     ") does not exist.")
      message("Using '", binary, "' as seqgen binary")
      assert_that(is.character(binary) && length(binary) == 1)
      private$binary <- binary

      super$initialize(priority)
    },
    call = function(args) {
      suppressWarnings(results <- system2(private$binary, args, stdout = TRUE))
      results
    },
    simulate = function(model, parameters = numeric()) {
      # Simulate the ancestral trees
      tree_model <- generate_tree_model(model)
      tree_sim_data <- simulate.coalmodel(tree_model, pars = parameters)
      trees <- tree_sim_data$trees
      assert_that(!is.null(trees))

      cmd_store <- new.env()
      cmd_store$list <- list()
      length(cmd_store$list) <- length(trees)

      # Call seq-gen for each locus group
      seg_sites <- lapply(1:length(trees), function(locus_group) {

        seqgen_args <-
          sg_generate_opts(model, parameters, locus_group,
                           sample_seed(length(trees[[locus_group]])))

        if (length(seqgen_args) == 1) {
          locus_length <- get_locus_length_matrix(model)[locus_group, 3]
        } else {
          locus_length <- get_locus_length_matrix(model)[locus_group,
                                                         c(1, 3, 5)]
        }

        cmd_store$list[[locus_group]] <- rep(NA, length(seqgen_args))

        # Call seq-gen of each trio locus
        group_loci <- lapply(seq(along = seqgen_args), function(trio_locus) {
          seqgen_file <- tempfile("seqgen")
          cmd <- paste(seqgen_args[trio_locus],
                       trees[[locus_group]][trio_locus])

          cmd_store$list[[locus_group]][trio_locus] <- paste("seq-gen", cmd)
          sim_output <- self$call(cmd)

          parse_seqgen_output(sim_output,
                              individuals = sum(get_sample_size(model, TRUE)),
                              locus_length = locus_length[trio_locus],
                              locus_number = get_locus_number(model,
                                                              locus_group,
                                                              TRUE),
                              outgroup_size = get_outgroup_size(model, TRUE),
                              calc_segsites = requires_segsites(model))
        })

        if (length(group_loci) == 1) return(group_loci[[1]])

        assert_that(length(group_loci) == 3)
        create_locus_trio(group_loci[[1]], group_loci[[2]], group_loci[[3]])
      })
      unlink(unlist(trees))

      # Generate the summary statistics
      seg_sites <- unlist(seg_sites, recursive = FALSE)

      if (requires_segsites(model)) {
        sequences <- NULL
      } else {
        sequences <- seg_sites
        seg_sites <- NULL
      }

      # Return the commands
      cmds <- list(trees = tree_sim_data$cmds,
                   seqgen = cmd_store$list)

      calc_sumstats_from_sim(seg_sites, NULL, sequences, model,
                             parameters, cmds, self)
    },
    get_cmd = function(model) {
      c(trees = get_cmd(generate_tree_model(model)),
        sequence = paste("seqgen", sg_generate_opts(model, NULL, 1, 0, TRUE),
                         collapse = " "))
    },
    get_info = function() c(name = "seqgen", binary = private$binary)
  )
)

has_seqgen <- function() !is.null(simulators[["seqgen"]])


#' Simulator: seq-gen
#'
#' This allows you to use seq-gen to simulate finite sites mutation models.
#' When using seq-gen, coala will simulate ancestral tress using the other
#' simulators, and call seq-gen to simulate finite sites mutations using the
#' trees. Seq-gen has a low priority, but will always be used when finite
#' sites mutation models are used.
#'
#' @section Installation:
#' You need to download the program from
#' \url{http://tree.bio.ed.ac.uk/software/seqgen/}
#' and compile the binary prior to invoking this function.
#' On Debian-based systems, you can alternatively install the package
#' 'seg-gen'.
#'
#' @references
#' Andrew Rambaut and Nicholas C. Grassly.
#' Seq-Gen: an application for the Monte Carlo simulation of DNA sequence
#' evolution along phylogenetic trees.
#' Comput Appl Biosci (1997) 13 (3): 235-238
#' doi:10.1093/bioinformatics/13.3.235
#'
#' @param binary The path of the seqgen binary that will be used
#'  for simulations. If none is provided, coala will look for a
#'  binary called 'seqgen' or 'seq-gen' using the PATH variable.
#' @inheritParams simulator_ms
#' @name simulator_seqgen
#' @family simulators
#' @export
#' @examples
#' \dontrun{activate_seqgen("./bin/seqgen")}
activate_seqgen <- function(binary = NULL, priority = 100) {
  register_simulator(seqgen_class$new(binary, priority))
  reset_cache()
  invisible(NULL)
}
