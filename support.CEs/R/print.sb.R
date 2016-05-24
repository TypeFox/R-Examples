print.sb <- function(x, ...) 
{
  if (x$sb$input$operation == "construct") {
    nalt <- x$design.information$nalternatives

    cat("\nChoice sets:\n")
    for (i in 1:nalt) {
      cat("alternative", i, "in each choice set\n")
      print(x$alternatives[[i]], ...)
      cat("\n")
    }

    cat("Candidate design:\n")
    dimnames(x$candidate) <- list(c(1:nrow(x$candidate)),
                                  c(1:ncol(x$candidate)))
    print(x$candidate, ...)
    cat("\n")

    cat("Design information:\n")
    cat("number of blocks =", x$design.information$nblocks, "\n")
    cat("number of questions per block =", x$design.information$nquestions, "\n")
    cat("number of alternatives per choice set =", x$design.information$nalternatives, "\n")
    cat("number of attributes per alternative =", x$design.information$nattributes, "\n")

    cat("\nMessage from the website 'Discrete Choice Experiments':\n")
    cat(x$sb$messages, "\n\n")

    cat(x$sb$version, "\n\n")
  } else {
    cat("\nResults of checking choice sets\n\n")

    cat(x$sb$messages, "\n\n")

    cat(x$sb$version, "\n\n")
  }

   invisible(x)
}
