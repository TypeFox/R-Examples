print.out_all <- function(x, ...) {

  n.out <- 0
  for (ilist in 1:length(x)) {
    if (substr(names(x[ilist]),1,4) == "out_") {  # then print component

      n.out <- n.out + 1
      if (nchar(x[ilist]) > 0) {
        if (n.out == 1) {
          if (is.null(options()$knitr.in.progress)) cat("\n")
        }
        else 
          cat("\n\n")

        if (substr(names(x[ilist]),5,9) != "title")
          for (i in 1:length(x[[ilist]])) cat(x[[ilist]][i], "\n")
        else
          cat(x[[ilist]]) 
      }

    }
  }

  cat("\n")

}
