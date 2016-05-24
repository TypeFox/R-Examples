# This files provides general urilires for the package


#' @import stringr
NULL


# collapse
collapse <- function(x,collapse='') paste(x,collapse = collapse)

# qw 
qw <- function (...) as.character(match.call())[-1]




#   `%|%` =
#   function(left, right){
#     subsright = substitute(right)
#     lazyright = lazy(right)
#     if(".." %in% all.vars(subsright))
#       lazy_eval(lazyright, list(.. = left))
#     else {
#       if(is.call(subsright)) {
#         lsr = as.list(subsright)
#         do.call(
#           as.character(
#             lsr[1]),
#           c(substitute(left),
#             lsr[-1]),
#           envir = lazyright$env)}
#       else{
#         if(is.function(right))
#           right(left)
#         else
#           stop("Don't know how to pipe THAT!")}}}
