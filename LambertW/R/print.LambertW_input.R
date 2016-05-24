#' @rdname LambertW_input_output-methods
#' @description
#' \code{print.LambertW_input} prints an overview of the input object.
#' @export

print.LambertW_input <- function(x, ...) {
  if (x$user.defined) {
    cat(" Note: This is a user-defined Lambert W x F distribution. \n")
    cat(" * * * * * * * * \n")
  }
  cat(" Input distribution: ")
  cat(x$distname)
  cat("\n with parameters: ")
  if (x$user.defined) {
    beta.names <- names(x$beta)
  } else {
    beta.names <- get_beta_names(x$distname)
  }
  cat(paste(paste0(beta.names, "=", round(x$beta, 3)), collapse = ", "))
  cat("\n")
} 