#' Tadaa, t-Test!
#'
#' @param data A \code{data.frame}.
#' @param response The response variable (dependent).
#' @param group The group variable, usually a \code{factor}.
#' @param direction Test direction, like \code{alternative} in \link{t.test}.
#' @param paired If \code{TRUE}, a paired t.test is performed with approproate power calculation.
#' @param na.rm If \code{TRUE} (default), missing values are dropped.
#' @inheritParams tadaa_aov
#' @return A \code{data.frame}, optionally markdown'd
#' @import pixiedust
#' @import stats
#' @importFrom car leveneTest
#' @family Tadaa-functions
#' @export
#' @examples
#' df <- data.frame(x = runif(100), y = sample(c("A", "B"), 100, TRUE))
#' tadaa_t.test(df, x, y)
#'
#' df <- data.frame(x = runif(100), y = c(rep("A", 50), rep("B", 50)))
#' tadaa_t.test(df, x, y, paired = TRUE)
tadaa_t.test <- function(data, response, group, direction = "two.sided",
                         paired = FALSE, na.rm = TRUE, print = "console") {

  response <- deparse(substitute(response))
  group    <- deparse(substitute(group))

  # Check the type of the group
  if (is.factor(data[[group]])) {
    groups <- levels(data[[group]])
  } else {
    groups <- unique(data[[group]])
  }

  # Subset groups of response
  x <- data[data[[group]] == groups[1], ][[response]]
  y <- data[data[[group]] == groups[2], ][[response]]

  # Kick out NAs if specified
  if (na.rm) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }

  # Get n for each group
  n1   <- length(x)
  n2   <- length(y)

  # levene
  levene    <- broom::tidy(car::leveneTest(data[[response]], data[[group]], center = "mean"))
  var.equal <-  ifelse(levene$p.value[[1]] <= .1, FALSE, TRUE)

  # t.test
  test <- broom::tidy(t.test(x = x, y = y, direction = direction,
                             paired = paired, var.equal = var.equal))

  # Additions
  test$d     <- effect_size_t(data = data, response = response, group = group, na.rm = na.rm)
  if (paired) {
    test$power <- pwr::pwr.t.test(n = n1, d = test$d, alternative = direction, type = "paired")$power
  } else {
    test$power <- pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = test$d, alternative = direction)$power
  }

  output <- pixiedust::dust(test)
  output <- pixiedust::sprinkle_colnames(output,
                                         statistic = "t", p.value = "p", parameter = "df",
                                         conf.low = "conf_low", conf.high = "conf_high")

  if ("estimate" %in% output$body$col_name) {
    output <- pixiedust::sprinkle_colnames(output, estimate = "Differenz")
  }
  if ("estimate1" %in% output$body$col_name) {
    output <- pixiedust::sprinkle_colnames(output, estimate1 = groups[[1]], estimate2 = groups[[2]])
  }

  output <- pixiedust::sprinkle(output, cols = "p.value", fn = quote(pixiedust::pvalString(value)))
  output <- pixiedust::sprinkle(output, round = 3)

  if (!(print %in% c("console", "hmtl", "markdown"))) {
    stop("Print method must be 'console', 'html' or, 'markdown'")
  }

  return(pixiedust::sprinkle_print_method(output, print_method = print))
}
