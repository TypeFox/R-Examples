booktabs <- function(x, ...) {
  if(!require(xtable, quietly = TRUE))
      stop("xtable not installed.")
  UseMethod("booktabs")
}

booktabs.default <- function(x, align = c("@{}l", rep("r", ncol(x)-1),
                                          "r@{}"), ...) {
  if(!require(xtable))
      stop("xtable not installed.")
  x <- xtable(x, align = align, ...)
  class(x) <- c("booktabs", class(x))
  x
}

booktabs.xtable <- function(x, ...) {
  if(!require(xtable))
      stop("xtable not installed.")
  a <- attr(x, "align")
  if(substring(a[1], 1, 3) != "@{}") {
    a[1] <- paste0("@{}", a[1])
    a[length(a)] <- paste0(a[length(a)], "@{}")
    attr(x, "align") <- a
  }
  class(x) <- c("booktabs", class(x))
  x
}

print.booktabs <- function(x, ...) {
  cat("\\begin{center}\n")
  xtable::print.xtable(x,
                       floating = FALSE,
                       hline.after = NULL,
                       booktabs = TRUE,
                       add.to.row = list(pos = list(-1, 0, nrow(x)),
                                         command = c('\\toprule ',
                                                     '\\midrule ',
                                                     '\\bottomrule '))) -> res
  cat("\\end{center}\n")
  invisible(res)
}
