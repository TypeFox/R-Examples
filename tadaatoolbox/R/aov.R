#' Tadaa, anova!
#'
#' @param formula Formula for model, passed to \code{aov}.
#' @param data Data for model.
#' @param show_effect_size If \code{TRUE} (default), effect size is appended as a column.
#' @param print Print method. Passed to \link[pixiedust]{sprinkle_print_method} as of now.
#' @return A \code{dust} object, depending on \code{print}.
#' @export
#' @family Tadaa-functions
#' @import stats
#' @examples
#' tadaa_aov(stunzahl ~ jahrgang, data = ngo)
tadaa_aov <- function(formula, data = NULL, show_effect_size = TRUE, print = "console"){

  model <- broom::tidy(aov(formula = formula, data = data))

  if (show_effect_size) {
    resid <- dplyr::last(model$sumsq)
    model$part.eta.sq <- model$sumsq / (resid + model$sumsq)
  }

  output <- pixiedust::dust(model)
  output <- pixiedust::sprinkle_colnames(output, statistic = "F")
  output <- pixiedust::sprinkle(output, col = 6, fn = quote(pixiedust::pvalString(value)))
  output <- pixiedust::sprinkle(output, col = 3:4, round = 2)
  output <- pixiedust::sprinkle(output, round = 3)

  if (!(print %in% c("console", "hmtl", "markdown"))) {
    stop("Print method must be 'console', 'html' or, 'markdown'")
  }

  return(pixiedust::sprinkle_print_method(output, print_method = print))
}
