##' Print labels
##'
##' A function to print standardised formatted labels without having
##' messy codes in the functions.
##'
##' @param label The label to be printed
##' @param span Whether the dash should span the whole width of the
##' screen(80 characters)
##' @param width The width of the screen.
##' @export
##' @return The formatted print

## TODO (Michael): Need to wrap the label
printLab = function(label, span = FALSE, width = getOption("width")){
  nc = nchar(label)
  sides = (width - nc)/2 - 3
  if(span){
    pre = paste(c("\n\n", rep("-", width), "\n"), collapse = "")
    post = paste(c("\n", rep("-", width), "\n\n"), collapse = "")
  } else {
    pre = paste(c("\n\n", rep(" ", sides), rep("-", nc + 6),
                   rep(" ", sides), "\n"), collapse = "")
    post = paste(c("\n", rep(" ", sides), rep("-", nc + 6),
                    rep(" ", sides), "\n\n"), collapse = "")
  }
  sandwich = paste(c(rep(" ", sides), "** ", label, " **",
                      rep(" ", sides)), collapse = "")
  cat(paste(pre, sandwich, post, sep = ""))
}
