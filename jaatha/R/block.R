block_class <- R6::R6Class("Block",
  private = list(border = NULL),
  public = list(
    initialize = function(border, cut) {
      assert_that(is.matrix(border))
      assert_that(ncol(border) == 2)
      assert_that(nrow(border) >= 1)
      assert_that(all(border[, 1] < border[, 2]))
      if (cut)  {
        border[border < 0] <- 0
        border[border > 1] <- 1
      }
      assert_that(all(border >= 0 & border <= 1))
      private$border <- border
    },
    get_border = function() private$border,
    print = function() print(private$border),
    print_border = function(jaatha) {
      lower <- denormalize(private$border[, 1], jaatha)
      upper <- denormalize(private$border[, 2], jaatha)
      paste0(round(lower, 3), "-", round(upper, 3), collapse = " x ")
    },
    includes = function(point) {
      all(private$border[, 1] - 1e-15 <= point &
            point <= private$border[, 2] + 1e-15)
    },
    get_middle = function() {
      m <- (private$border[, 2] - private$border[, 1]) / 2 +
        private$border[, 1]
      names(m) <- rownames(private$border)
      m
    },
    get_corners = function() {
      corners <- expand.grid(lapply(1:nrow(private$border), function(i) {
        private$border[i, , drop = FALSE] #nolint
      }), KEEP.OUT.ATTRS = FALSE) #nolint
      colnames(corners) <- rownames(private$border)
      as.matrix(corners)
    },
    sample_pars = function(number, add_corners = FALSE) {
      "Generates random parameter combinations inside the block's range"
      assert_that(is.count(number))
      assert_that(is_single_logical(add_corners))

      # Sample random simulation parameters
      par_number <- nrow(self$get_border())
      random_pars <- matrix(runif(par_number * number,
                                  min = self$get_border()[, 1],
                                  max = self$get_border()[, 2]),
                            number, par_number, byrow = TRUE)

      # Add corners if requested
      if (add_corners) random_pars <- rbind(random_pars, self$get_corners())
      assert_that(all(apply(random_pars, 1, self$includes)))
      random_pars
    },
    get_interior = function(x = 0.1) {
      "removes the outer `x` percent at each side of the block "
      border <- self$get_border()
      outer <- apply(border, 1, diff) * x
      border + cbind(outer, -outer)
    }
  )
)

create_block <- function(border, cut = FALSE) block_class$new(border, cut)
