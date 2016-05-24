# Copyright Seth Wenchel 2015

#' @param  .data A tbl or something that can be coerced into one
#' @param ... conditions that will be passed to dplyr::filter
#' @param false_fun A function or functional that will be applied to the data that doesn't pass the supplied filters (the scion)
#' @param false_name optional, the name of the object to which the scion will be assigned.
#' @param false_env optional, the environment into which the scion will be assigned. If specified, false_name must also be specified.
#'        If unspecified (default), scions will be placed into the internal package environment.
#' @return A tbl whose rows have passed the stated conditions
#' @details \code{.data} will be split into two chunks based on the conditions. The scion will be passed through \code{false_fun} and then either placed on
#' the package's internal stack or assigned as specified by \code{false_name} and \code{false_env}.
#' @examples
#' library(dplyr)
#' aframe <- data.frame(zed = runif(100))
#' set_to_zero <- . %>% mutate(zed = 0)
#' aframe %>% scion(zed >0.5, false_fun=set_to_zero) %>% mutate(zed=1) %>% graft
#'
#' @export
#' @title scion
#' @author Seth Wenchel
#' @importFrom dplyr filter
#' @importFrom dplyr setdiff
#' @importFrom magrittr %>%
scion <- function(.data, ..., false_fun, false_name, false_env){
  .data %>% dplyr::filter(...) -> tru
  falls <- dplyr::setdiff(.data, tru)
  if(!missing(false_fun))
    falls <- false_fun(falls)
  if(!missing(false_name)){
    if(!missing(false_env))
      assign(false_name, falls, envir = false_env)
    else
      assign(false_name, falls, envir = .pkgenv)
  }
  else
    if(!missing(false_env))
      stop("false_env specified but not false_name.")
  .push(falls)
  return(tru)
}


#' @param .data A tbl or something that can be coerced into one
#' @param data2 A tbl or something that can be coerced into one
#' @param combine_fun optional, A function that will combine two tbls such as full_join or bind_rows
#' @description Graft one dataset onto another
#' @details Graft requires two data objects.  The first must be provided by the user. The second can either be passed
#' in or automatically pulled off of the package's internal stack of scions. These will be combined accoring to the following rules in order:
#' \itemize{
#'   \item If either dataset has zero rows, the other dataset will be returned.
#'   \item If combine_fun is specifed, \code{combine_fun(.data, data2)} will be called
#'   \item If all column names match, a row bind will occur
#'   \item If at least some column names match, a full join will occur
#'   \item If both have the same number of rows a column bind will be performed
#'  }
#' @return A single tbl object
#' @export
#' @author Seth Wenchel
#' @title Graft
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr intersect
#' @importFrom dplyr full_join
graft <- function(.data, combine_fun, data2){
  if(missing(data2))
    data2 <- .pop()

  if(length(.data)==0)
    return(data2)
  if(length(data2)==0)
    return(.data)
  if(!missing(combine_fun))
    return(combine_fun(.data,data2))
  if(names(.data)==names(data2))
    return(dplyr::bind_rows(.data, data2))
  if(length(dplyr::intersect(names(.data),names(data2)))>0)
    return(dplyr::full_join(.data, data2))
  if(nrow(.data)==nrow(data2))
    return(dplyr::bind_cols(.data, data2))


  # same # columns & names <=> rbind
  # same # rows & no names match <=> cbind
  # at least one column names matches <=> full_join
}

#' Remove all objects from the stack by deleting them from memory.
#' @export
clear_stack <- function(){
  vars <-  ls(envir = .pkgenv)
  rm(list = vars[which(substr(vars, 1,10)=="stack_obj_")], envir = .pkgenv)
  assign("objects_on_stack", 0,.pkgenv)
}



