#' @importFrom assertthat assert_that is.number
stat_omega_class <- R6Class("stat_omega", inherit = sumstat_class,
  private = list(
    req_segsites = TRUE,
    binary = NULL,
    min_win = NULL,
    max_win = NULL,
    grid = NULL,
    create_empty_result = function(locus, locus_length) {
      data.frame(locus = locus, pos = locus_length / 2, omega = 0)
    }
  ),
  public = list(
    initialize = function(name, min_win, max_win,
                          grid, binary, transformation) {
      assert_that(is.number(min_win))
      private$min_win <- min_win
      assert_that(is.number(max_win))
      assert_that(min_win < max_win)
      private$max_win <- max_win
      assert_that(is.number(grid))
      private$grid <- grid

      if (identical(binary, "automatic")) {
        binary <- search_executable("OmegaPlus", envir_var = "OMEGAPLUS")
        if (is.null(binary)) stop("No binary for OmegaPlus found.")
      } else {
        assert_that(length(binary) == 1)
        assert_that(is.character(binary))
        assert_that(file.exists(binary))
      }
      private$binary <- binary
      super$initialize(name, transformation)
    },
    check = function(model) {
      if (has_trios(model)) {
        stop("OmegaPlus can not be calculated from locus trios")
      }
      if (any(get_locus_length(model, total = TRUE) < self$get_grid())) {
        stop("Grid value in stat_omega can not be larger than the locus length")
      }
      invisible(TRUE)
    },
    calculate = function(seg_sites, trees, files, model) {
      cur_wd <- getwd()

      tmp_dir <- tempfile("omegaprime")
      dir.create(tmp_dir)
      setwd(tmp_dir)
      grid <- self$get_grid()

      op_list <- lapply(seq(along = seg_sites), function(i) {
        locus_length <- get_locus_length(model, locus = i)

        # Return 0 if there are few SNPs
        if (ncol(seg_sites[[i]]) <= 2) {
          return(private$create_empty_result(i, locus_length))
        }

        # Adjust number of grid points if needed
        pos_rel <- get_positions(seg_sites[[i]])
        max_grid <- floor(diff(pos_rel[c(1, length(pos_rel))] * locus_length))
        if (grid > max_grid) grid <- max_grid - 2
        if (grid < 1) return(private$create_empty_result(i, locus_length))

        # Create an input file
        tmp_file <- tempfile("omega")
        cat(c("ms 10 1 -t 5",
              "3579 27011 59243",
              "",
              "//",
              conv_to_ms_output(seg_sites[[i]])),
              "", sep = "\n", file = tmp_file)

        # Execute OmegaPlus
        system2(private$binary,
                args = c("-name", i,
                         "-minwin", self$get_min_win(),
                         "-maxwin", self$get_max_win(),
                         "-grid", grid,
                         "-length", locus_length,
                         "-input", tmp_file),
                stdout = TRUE)
        unlink(tmp_file)

        # Parse the results
        self$parse_report(tmp_dir, grid, i)
      })

      unlink(tmp_dir, recursive = TRUE)
      setwd(cur_wd)

      do.call(rbind, op_list)
    },
    parse_report = function(dir, n_grid, locus) {
      op_file <- file.path(dir, paste0("OmegaPlus_Report.", locus))
      if (!file.exists(op_file)) stop("Calculation of omega failed.")
      values <- read.delim(op_file, header = FALSE, comment.char = "/")
      colnames(values) <- c("pos", "omega")
      assert_that(nrow(values) %% n_grid == 0)
      data.frame(locus = locus, values)
    },
    get_grid = function() private$grid,
    get_min_win = function() private$min_win,
    get_max_win = function() private$max_win
  )
)


#' Summary Statistic: Omega
#'
#' Calculates the Omega Statistic introduced by
#' Kim & Nielsen (2004) from the simulated data. The statistic is sensitive for
#' hard selective sweeps. To calculate
#' the statistic, coala relies on the command line program
#' \href{http://sco.h-its.org/exelixis/web/software/omegaplus/index.html}{OmegaPlus},
#' which needs to be downloaded and compiled manually in order to use the
#' statistic.
#'
#' @references
#' Linkage disequilibrium as a signature of selective sweeps.
#' Y. Kim and R. Nielsen (2004). Genetics, 167, 1513-1524.
#'
#' OmegaPlus: a scalable tool for rapid detection of selective
#' sweeps in whole-genome datasets.
#' N. Alachiotis, A. Stamatakis and P. Pavlidis (2012).
#' Bioinformatics Vol. 28 no. 17 2012, pages 2274-2275
#' doi:10.1093/bioinformatics/bts419
#'
#' @inheritParams sumstat_four_gamete
#' @param min_win The minimum distance from the grid point that a SNP must have
#'   to be included in the calculation of omega.
#' @param max_win The maximum distance from the grid point that a SNP must have
#'   to be included in the calculation of omega.
#' @param grid The number of points for which omega is calculated on each
#'   locus. Should be significantly lower than the locus length.
#' @param binary The path of the binary for OmegaPlus. If set to "automatic",
#'   coala will try to find a binary called "OmegaPlus" using the PATH
#'   environment variable.
#' @return A data frame listing of locus, genetic position and the
#'   calculated omega value.
#' @export
#' @template summary_statistics
#' @examples
#' \dontrun{
#' model <- coal_model(20, 1, 50000) +
#'   feat_recombination(50) +
#'   feat_mutation(1000) +
#'   feat_selection(strength_A = 1000, time = 0.03) +
#'   sumstat_omega()
#' stats <- simulate(model)
#' plot(stats$omega$omega, type = "l")}
sumstat_omega <- function(name = "omega", min_win = 100, max_win = 1000,
                          grid = 1000, binary = "automatic",
                          transformation = identity) {
  stat_omega_class$new(name, min_win, max_win, grid, binary, transformation)
}

has_omega <- function() {
  !is.null(search_executable("OmegaPlus", envir_var = "OMEGAPLUS"))
}
