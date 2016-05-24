# This code was contributed by Hadley Wickham

end_iteration <- function() stop('StopIteration', call.=FALSE)

iteration_has_ended <- function(e) {
  identical(conditionMessage(e), 'StopIteration')
}

new_iterator <- function(nextElem, ...) {
  structure(list(nextElem=nextElem, ...), class=c('abstractiter', 'iter'))
}

is.iterator <- function(x) inherits(x, 'iter')
