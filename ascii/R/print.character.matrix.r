##' div
##'
##' @keywords internal
##' @param x x
##' @param n n
div <- function(x, n = 2) {
  xx <- floor(x/n)
  c(xx, x - xx)
}

##' rep.char
##'
##' @keywords internal
##' @param x x
##' @param times times
rep.char <- function(x, times = 1) {
  paste(rep(x, times), collapse = "")
}

##' align.table
##'
##' @keywords internal
##' @param x x
##' @param align align
##' @param space space
align.table <- function(x, align, space = 1) {
  align <- expand(align, nrow(x), ncol(x), drop = FALSE)
  nchar <- apply(x, 2, nchar)
  max.nchar <- expand(apply(expand(nchar, nrow(x), ncol(x), drop = FALSE), 2, max) + space*2, nrow(x), ncol(x), drop = FALSE)
  diff.nchar <- max.nchar - nchar

  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (align[i, j] == "l")
        x[i, j] <- paste(x[i, j], paste(rep(" ", diff.nchar[i, j]), collapse = ""), collapse = "", sep = "")
      if (align[i, j] == "r")
        x[i, j] <- paste(paste(rep(" ", diff.nchar[i, j]), collapse = ""), x[i, j], collapse = "", sep = "")
      if (align[i, j] == "c") {
        sides <- div(diff.nchar[i, j])
        x[i, j] <- paste(paste(rep(" ", sides[1]), collapse = ""), x[i, j], paste(rep(" ", sides[2]), collapse = ""), collapse = "", sep = "")
      }
    }
  }
  x
}

##' print.character.matrix
##'
##' @param x x
##' @param vsep vsep
##' @param before_vsep before_vsep
##' @param after_vsep after_vsep
##' @param hsep hsep
##' @param csep csep
##' @param before_cell_content before_cell_content
##' @param after_cell_content after_cell_content
##' @param line_separator line_separator
##' @param line_separator_pos line_separator_pos
##' @param justify justify
##' @param right_alignment right_alignment
##' @param print print
##' @keywords internal
print.character.matrix <- function(x, vsep = "|", before_vsep = "", after_vsep = "", hsep = "-", csep = "+", before_cell_content = " ", after_cell_content = " ", line_separator = TRUE, line_separator_pos = NULL, justify = "l", space = 0, right_alignment = FALSE, print = TRUE) {

  # after et before cell_content
  x <- paste.matrix(before_cell_content, x, after_cell_content, sep = "", transpose.vector = TRUE)
  x <- align.table(x, justify, space)

  # dim
  nrowx <- nrow(x)
  ncolx <- ncol(x)

  # vseps
  if (is.null(dim(vsep))) {
    vseps <- expand(vsep, nrowx, ncolx+1, drop = FALSE)
  } else {
    vseps <- expand(vsep, nrowx, ncolx+1, what = "", drop = FALSE)
  }
  before_vseps <- expand(before_vsep, nrow(vseps), ncol(vseps), what = "")
  after_vseps <- expand(after_vsep, nrow(vseps), ncol(vseps), what = "")
  final_vseps <- paste.matrix(before_vseps, vseps, after_vseps, sep = "")

  # nchar
  ncharvseps <- nchar(vseps)
  ncharx <- nchar(x)
  ncharbefore <- nchar(before_vseps)
  ncharafter <- nchar(after_vseps)
  nchartot <- ncharafter[, 1:(ncol(ncharafter)-1)] + ncharx + ncharbefore[, -1]
  nchartotmax <- apply(nchartot, 2, max)

  # rows
  lines <- interleave.matrix(final_vseps, x, byrow = FALSE)
  lines <- matrix(lines[!is.na(lines)], nrow = nrowx)
  lines <- paste.matrix(lines, sep = "", collapse = "")

  # row separators
  row_lines <- NULL
  if (line_separator) {
    hseps <- expand(hsep, nrowx+1, ncolx)
    hseps <- t(apply(hseps, 1, function(x) Vectorize(rep.char, c("x", "times"))(x, nchartotmax)))
    if (ncolx == 1)
      hseps <- t(hseps)
    cseps <- expand(csep, nrowx+1, ncolx+1)
    row_lines <- paste.matrix(interleave.matrix(cseps, hseps, byrow = FALSE), collapse = "", byrow = FALSE)

    if (!is.null(line_separator_pos)) {
      line_separator_pos <- sort(unique(ifelse(line_separator_pos >= 0, line_separator_pos + 1, nrowx + 1 + line_separator_pos + 1)))
      line_separator_pos <- line_separator_pos[line_separator_pos <= nrowx + 1]
      row_lines[-line_separator_pos] <- NA
    }
  }

  # all together
  results <- interleave(row_lines, lines)
  results <- results[!is.na(results)]

  # right alignment of the table
  if (right_alignment) {
    ncharres <- nchar(results)
    nchardiff <- sapply(mapply(rep, times = max(ncharres) - ncharres, x = " "), paste, collapse = "")
    results <- paste(nchardiff, results, sep = "")
  }

  # cat
  if (print)
    cat(results, sep = "\n")
  invisible(results)
}
