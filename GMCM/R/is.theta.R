#' Check if parameters are valid
#'
#' Function to check whether the argument is coherent and in the correct
#' format.
#'
#' @param theta A list on the \code{theta}-form described in
#'   \code{\link{rtheta}}
#' @return \code{logical}. Returns \code{TRUE} if \code{theta} is coherent and
#'   in the correct format. Otherwise, the function returns \code{FALSE} with
#'   an accompanying warning message of the problem.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}
#' @examples
#' theta1 <- rtheta()  # Create a random correctly formatted theta
#' is.theta(theta1)
#'
#' theta2 <- rtheta(d = 3, m = 5)
#' theta2$m <- 6  # m is now incoherent with the number of components
#' is.theta(theta2)
#'
#' theta3 <- rtheta(d = 4, m = 2)
#' theta3$sigma$comp1[1, 2] <- 0  # Making the covariance matrix non-symmetric
#' is.theta(theta3)
#'
#' theta4 <- rtheta(d = 10, m = 10)
#' theta4$sigma$comp1[1, 1] <- 0  # Destroy positive semi-definiteness
#' is.theta(theta4)
#'
#' theta5 <- rtheta()
#' names(theta5) <- c("m", "d", "prop", "mu", "sigmas") # Incorrect names
#' is.theta(theta5)
#' @export
is.theta <- function(theta) {
  # Testing structure of theta
  if (!is.list(theta) | length(theta) != 5) {
    warning("theta is not a list of length 5")
    return(FALSE)
  }
  if (!is.list(theta[[4]])) {
    warning("theta[[4]] is not a list")
    return(FALSE)
  }
  if (!is.list(theta[[5]])) {
    warning("theta[[5]] is not a list")
    return(FALSE)
  }
  for (i in 1:2) {
    if (!is.numeric(theta[[i]]) | !length(theta[[i]]) == 1) {
      warning("theta[[",i,"]] is not a numeric vector of length 1")
      return(FALSE)
    }
  }
  # Testing mixture proportions
  if (length(theta[[3]]) != theta[[1]]) {
    warning("theta[[3]] is not a vector of length ", theta[[1]], " as defined",
            " by theta[[1]]")
    return(FALSE)
  }
  if (!all.equal(sum(theta[[3]]), 1)) {
    warning("The mixture proportions theta[[3]] does not sum to 1.")
    return(FALSE)
  }
  # Testing mean vectors
  if (!all(sapply(theta[[4]], length) ==  theta[[2]])) {
    warning("The length of the vectors in theta[[4]] does not equal ",
            theta[[2]], " as defined in theta[[2]].")
    return(FALSE)
  }
  # Testing covariance matrices
  if (length(theta[[5]]) != theta[[1]]) {
    warning("theta[[5]] is not a list of length ", theta[[1]], " as given by",
            " theta[[1]].")
    return(FALSE)
  }
  if (!all(c(sapply(theta[[5]], dim)) ==  theta[[2]])) {
    warning("The covariance matrices in theta[[5]] does not have dimensions ",
            theta[[2]], " times ", theta[[2]], " as given by theta[[2]].")
    return(FALSE)
  }
  if (!all(sapply(theta[[5]], isSymmetric))) {
    warning("Not all covariance matrices are symmetric.")
    return(FALSE)
  }
  is.PosDef <- function(x) {
    all(eigen(x)$values >= 0)
  }
  if (!all(sapply(theta[[5]], is.PosDef))) {
    warning("Not all covariance matrices are postive semi-definite.")
    return(FALSE)
  }
  if (!identical(names(theta), c("m", "d", "pie", "mu", "sigma"))) {
    warning('names(theta) does not equal c("m", "d", "pie", "mu", "sigma")')
    return(FALSE)
  }
  return(TRUE)
}
