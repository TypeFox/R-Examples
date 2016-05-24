print.bounds <- function(x, digits = max(2, getOption("digits") - 4), ...) { 
  cat("\nDeterministic bounds:\n\n")
  for(i in 1:length(x$bounds)){
    cat(names(x$bounds)[i], "\n")
    print.default(format(x$bounds[[i]], digits = digits),
                print.gap = 2, quote = FALSE, ...)
    cat("\n")
  }
  invisible(x)
}
