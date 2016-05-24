#' Calculates the parameters
#' \code{parameters} calculates the parameters
#' @keywords internal
#' @export
parameters <- function(d, x, k, n, psyfunguesslapses, funname,
                    parini, pariniset, guess, lapses, optimization, groups) {

  d %>% do(one_parameters(., x, k, n, psyfunguesslapses, funname,
                  parini, pariniset, guess, lapses, optimization, groups))
}
