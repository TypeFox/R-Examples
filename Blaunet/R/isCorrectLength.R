isCorrectLength <-
function(arg, argLength = 1) { #spits out an error if argument is of wrong length. returns true if there are no errors (NULL or correct length)
  if (!is.null(arg)){
    if (length(arg) != argLength) {sprintf('Option %s takes argument of at most length %s', arg, argLength)}
    else {TRUE}
  }
  else {TRUE}
}
