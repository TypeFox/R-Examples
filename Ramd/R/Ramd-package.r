#' Most R projects are a collection of loosely organized scripts that include
#' each other through haphazard `source` calls. However, this often pollutes
#' the global namespace and makes it easy for the developer to be stranded in
#' a limbo in which they are unsure of the present state of the program.
#' The result is unnecessary time wasted in restarting and re-executing
#' past work.
#' 
#' Ramd aims to solve this problem by explicitly declaring dependencies.
#' See the \code{define} function.
#'
#' @name Ramd
#' @docType package
NULL
