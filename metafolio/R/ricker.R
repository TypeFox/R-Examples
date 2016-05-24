#' A simple Ricker model
#'
#' @param spawners Spawner abundance
#' @param a Ricker productivity parameter. Recruits are e^a at the origin.
#' @param b Ricker density dependent parameter.
#' @return Returns the number of recruits.
#' @export
#' @examples
#' S <- seq(100, 1000, length.out = 100)
#' R <- ricker(S, a = 1.9, b = 900)
#' plot(S, R)

ricker <- function(spawners, a, b) {
  spawners * exp(a * (1 - spawners / b))
}
