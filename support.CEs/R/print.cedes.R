print.cedes <- function(x, ...)
{
# Name: print.cedes
# Title: print() for S3 class "cedes"
# Arguments:
#  x            An object of S3 class "cedes"
#  ...          Arguments passed to format()


# set variable

  nalt <- x$design.information$nalternatives # nunber of alternatives per choice set

# display choice sets

  cat("\nChoice sets:\n")
  for(i in 1:nalt) {  # ecah alternative
    cat("alternative", i, "in each choice set\n")
    print(x$alternatives[[i]], ...)
    cat("\n")
  }

# display chandidate design

  cat("Candidate design:\n")
  print(x$candidate, ...)
  cat("\n")

# display basic infomation regarding choice sets

  cat("Design information:\n")
  cat("number of blocks =", x$design.information$nblocks, "\n")
  cat("number of questions per block =", x$design.information$nquestions, "\n")
  cat("number of alternatives per choice set =", x$design.information$nalternatives, "\n")
  cat("number of attributes per alternative =", x$design.information$nattributes, "\n\n")

  invisible(x)
}

