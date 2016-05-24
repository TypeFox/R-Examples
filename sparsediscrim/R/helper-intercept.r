#' Removes the intercept term from a formula if it is included
#'
#' Often, we prefer not to have an intercept term in a model, but user-specified
#' formulas might have included the intercept term. In this case, we wish to
#' update the formula but without the intercept term. This is especially true in
#' numerous classification models, where errors and doom can occur if an
#' intercept is included in the model.
#' 
#' @param formula a model formula to remove its intercept term
#' @param data data frame
#' @return formula with no intercept term
#' @examples
#' iris_formula <- formula(Species ~ .)
#' sparsediscrim:::no_intercept(iris_formula, data = iris)
no_intercept <- function(formula, data) {
  # The 'terms' must be collected in case the dot (.) notation is used
  update(formula(terms(formula, data = data)), . ~ . - 1)
}
