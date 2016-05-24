#' @include model.R
#' @importFrom R6 R6Class
locus_class <- R6Class("locus", inherit = model_part,
  private = list(
    number = 0,
    length = NA
  ),
  public = list(
    initialize = function(locus_length, locus_number = 1) {
      assert_that(is.numeric(locus_length))
      assert_that(length(locus_length) == 1 || length(locus_length) == 5)
      assert_that(all(locus_length >= 0))
      private$length <- locus_length

      assert_that(is.numeric(locus_number))
      assert_that(locus_number > 0)
      assert_that(length(locus_number) == 1)
      private$number <- locus_number
    },
    get_length = function(trios = FALSE) {
      if (!trios) return(sum(private$length))

      if (length(private$length) == 1) {
        locus_length <- c(0, 0, private$length, 0, 0)
      } else if (length(private$length) == 5) {
        locus_length <- private$length
      } else stop("Failed to get locus length")

      names(locus_length) <- c("length_l", "length_il", "length_m",
                               "length_ir", "length_r")
      locus_length
    },
    get_number = function() private$number,
    print = function() {
      cat(self$get_number(), ifelse(self$get_number() == 1, "locus", "loci"),
          "of length", self$get_length(), "\n")
    }
  )
)


is.locus <- function(locus) any("locus" == class(locus))


#' Loci
#'
#' This functions adds one or more loci to a model. A locus is a continuous
#' stretch of DNA of a given length. All loci are simulated independently of each
#' other, and are genetically unlinked. A model can contain a
#' large number of different loci created with \code{locus_single}. This will,
#' however, slow down the simulation. For performance reasons, it is
#' better to add the same number of loci with averaged length using
#' \code{locus_averaged} if this simplification is justifiable. Both can also be
#' combined in a single model. In the results,
#' the summary statistics for the loci are returned in order in which they
#' are added to the model.
#'
#' @param length The length of the locus in base pairs.
#' @seealso For adding three loci which are linked to each other:
#'  \code{\link{locus_trio}}
#' @examples
#' # A model with one locus of length 1005 bp:
#' coal_model(10) + locus_single(1005)
#' # This is equivalent to:
#' coal_model(10, 1, 1005)
#'
#' # A model can contain multiple loci:
#' coal_model(5) + locus_single(100) + locus_single(200) + locus_single(300)
#' # Or more efficient with averaged length:
#' coal_model(5) + locus_averaged(3, 200)
#' # Or equivalently:
#' coal_model(5, 3, 200)
#'
#' # Single and averaged loci can also be combined arbitrarily:
#' coal_model(15) + locus_averaged(10, 150) + locus_single(250)
#' coal_model(15, 10, 150) + locus_single(250) + locus_averaged(10, 350)
#' @name locus
#' @aliases loci
NULL

#' @describeIn locus Adds a single locus.
#' @export
locus_single <- function(length) {
  locus_class$new(length, 1)
}


#' @describeIn locus Adds multiple loci with equal length.
#' @param number The number of loci to add.
#' @export
locus_averaged <- function(number, length) {
  locus_class$new(round(length), number)
}


#' Locus Trios
#'
#' This functions adds a group of three loci to the model that are genetically
#' linked to each other. They are still unlinked to all other loci or locus trios
#' in the model. Simulating linked loci that are far apart from each other can
#' be very slow. Please mind that mutation and recombination rates for locus
#' trios are rates per trio and not per locus, i.e. they account for mutations
#' that occur on the tree loci and the sequences in-between them together.
#'
#' @inheritParams locus
#' @param locus_length An integer vector of length 3, giving the length of each
#'   of the three loci (left, middle and right).
#' @param distance A vector of two, giving the distance between left and middle,
#'   and middle an right locus, in base pairs.
#' @export
#' @seealso For adding unlinked loci: \code{\link{locus}}
#' @examples
#' # A model with one locus trio
#' coal_model(25) +
#'   locus_trio(locus_length=c(1250, 1017, 980), distance=c(257, 814))
#'
#' # Ten identical locus trios:
#' coal_model(25) +
#'   locus_trio(locus_length=c(1250, 1017, 980), distance=c(257, 814), number = 10)
#'
#' # Two different ones:
#' coal_model(25) +
#'   locus_trio(locus_length=c(1000, 500, 900), distance=c(200, 400)) +
#'   locus_trio(locus_length=c(700, 500, 800), distance=c(350, 150))
locus_trio <- function(locus_length = c(left = 1000,
                                        middle = 1000,
                                        right = 1000),
                       distance = c(left_middle = 500,
                                    middle_right = 500),
                       number = 1) {

  stopifnot(length(locus_length) == 3)
  if (!is.null(names(locus_length))) {
    locus_length <- locus_length[c("left", "middle", "right")]
  }
  stopifnot(length(distance) == 2)
  if (!is.null(names(distance))) {
    distance <- distance[c("left_middle", "middle_right")]
  }

  locus_class$new(locus_length = c(locus_length, distance)[c(1, 4, 2, 5, 3)],
                  locus_number = number)
}


# Converts a position on the middle locus to the relative position
# on the simulated stretch
conv_middle_to_trio_pos <- function(pos, model, locus,
                                    relative_out = TRUE, relative_in = TRUE) {

  llm <- get_locus_length_matrix(model)
  group <- get_locus_group(model, locus)

  if (relative_in) pos <- pos * llm[group, 3]
  pos <- pos + llm[group, 1] + llm[group, 2]
  if (relative_out) pos <- pos / sum(llm[group, 1:5])

  names(pos) <- NULL
  pos
}
