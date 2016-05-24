#' Generate a matrix of similar characters
#'
#' This function prints a matrix of characters which are very similar to each
#' other.
#' @param x a character vector of length 2 (usually two similar characters)
#' @param n the total number of characters in the matrix
#' @param nrow the number of rows
#' @return a character matrix on the screen
#' @author Yihui Xie <\url{http://yihui.name}>
#' @export
#' @examples
#' char_gen()
#'
#' char_gen(c('O', 'Q'))
char_gen = function(x = c("V", "W"), n = 300, nrow = 10) {
  cat(apply(matrix(sample(c(x[1], rep(x[2], n - 1))), nrow),
            1, paste, collapse = ""), sep = "\n")
}
