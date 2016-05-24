#' @rdname LambertW_input_output-methods
#' @description
#' \code{print.LambertW_output} prints an overview of the output object.
#' @export

print.LambertW_output <- function(x, ...) {
  
  if (is.null(x$theta$gamma)) {
    x$theta$gamma <- 0
  }
  if (is.null(x$theta$delta)) {
    x$theta$delta <- 0
  }
  if (is.null(x$theta$alpha)) {
    x$theta$alpha <- 1
  }
  cat(" Input distribution: ")
  distname.input <- gsub("Lambert W x ", "", x$distname)
  cat(distname.input)
  cat("\n")
  # cat('Lambert W type ('h' same tails; 'hh' different tails; 's' skewed):
  # ')
  cat(" Output distribution: ")
  pre.tex <- NULL
  if (x$type == "s") {
    pre.tex <- "skewed"
  } else if (x$type == "h") {
    pre.tex <- "heavy-tail (one parameter)"
  } else if (x$type == "hh") {
    pre.tex <- "heavy-tail (two parameters)"
  }
  cat(paste(pre.tex, x$distname.with.beta))
  cat("\n with (input) parameters: ")
  cat(paste(paste0(get_beta_names(distname.input), "=", round(x$theta$beta, 3)), collapse = ", "))
  cat("\n and transformation parameters: ")
  if (x$type == "s") {
    cat("gamma=", round(x$theta$gamma, 3), sep = "")
  } else if (x$type == "h") {
    cat("delta =", round(x$theta$delta, 3))
  } else if (x$type == "hh") {
    cat(paste0(c("delta_l", " delta_r"), "=", round(x$theta$delta, 3)), sep = " and")
  }
  cat("\n")
  if (x$theta$alpha != 1) {
    cat("alpha =", round(x$theta$alpha, 3), sep = "")
    cat("\n")
  }
} 