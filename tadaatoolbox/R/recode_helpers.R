#' Easily cluster numeric vectors in likert-like classes
#'
#' @param x Vector to be clustered.
#' @param classes Number of classes. Defaults to 3, can also be 5.
#' @param method How should the classes be calculated? Defaults to \code{quantiles}, can also be
#' \code{means} for mean and standard deviation.
#' @return An ordered \code{factor} with \code{classes} levels. And descriptive labels.
#' @importFrom car recode
#' @export
#' @examples
#' tadaa_likertize(x = runif(100, 0, 10), classes = 3, method = "quantiles")
#' tadaa_likertize(x = runif(100, 0, 10), classes = 3, method = "meansd")
tadaa_likertize <- function(x, classes = 3, method = "quantiles"){

  if (classes == 3) {
    if (method == "quantiles") {
      quantiles <- quantile(x, (1:2) / classes)
      recodes <- paste0("lo:", quantiles[[1]], " = 1; ",
                        quantiles[[1]], ":", quantiles[[2]], " = 2; ",
                        quantiles[[2]], ":hi = 3")
    } else if (method == "meansd") {
      recodes <- paste0("lo:", mean(x) - sd(x), " = 1; ",
                        mean(x) -  sd(x), ":", mean(x) + sd(x), " = 2; ",
                        mean(x) + sd(x), ":hi = 3")
    }
      xx <- car::recode(x, recodes = recodes)
      xx <- factor(xx, labels = c("niedrig", "mittel", "hoch"), ordered = TRUE)
      return(xx)
  } else if (classes == 5) {
    stop("not yet implemented")
  }
}

#' Convenience functions for interval recodes
#'
#' Get recode assigments for even intervals of discrete numeric values compatible
#' with \link[car]{recode}.
#'
#' @param from,to A \code{numeric value} for the beginning and the end of the interval.
#' @param width The width of the interval, e.g. 5 (default) for intervals 0-5.
#'
#' @return A \code{character} vector of recode assignments compatible with \link[car]{recode}
#' @export
#' @examples
#' \dontrun{
#' x       <- round(runif(100, 0, 100), 0)
#' recodes <- generate_recodes(0, 100, 10)
#'
#' library(car)
#' recode(x, recodes = recodes)
#'}
generate_recodes <- function(from, to, width = 5){
  paste(sapply(seq(from, to - width, 2 * width), function(x) {
    first  <- paste0(x, ":", x + width, "='", x, "-", x + width, "'")
    second <- paste0(x + width + 1, ":", x + (2 * width), "='", x + width + 1, "-", x + (2 * width), "'")
    paste(first, second, sep = "; ", collapse = "; ")
  }), collapse = "; ")
}

#' Convenience functions for interval recodes
#'
#' Get interval labels for even intervals of discrete numeric values compatible
#' with \link{cut}.
#'
#' @param from,to A \code{numeric value} for the beginning and the end of the interval.
#' @param width The width of the interval, e.g. 5 (default) for intervals 0-5.
#'
#' @return A \code{character} vector of interval labels compatible with \link{cut}
#' @export
#' @examples
#' \dontrun{
#' x       <- round(runif(100, 0, 100), 0)
#' labels  <- interval_labels(0, 100, 10)
#'
#' cut(x, breaks = seq(0, 100, 10), labels = labels)
#'}
interval_labels <- function(from, to, width = 5){
  labs <- lapply(seq(from, to - width, 2 * width), function(x) {
    c(paste0(x, "-", x + width), paste0(x + width + 1, "-",  x + (2 * width)))
  })
  return(unlist(labs))
}
