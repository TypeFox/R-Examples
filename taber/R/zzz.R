# Copyright Seth Wenchel 2015

.pkgenv <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname){
  assign("objects_on_stack", 0,.pkgenv)
}



.push <- function(.data){
  count <- get("objects_on_stack",envir = .pkgenv)+1
  assign(paste0("stack_obj_",count),.data, envir = .pkgenv)
  assign("objects_on_stack", (count),.pkgenv)
}

.pop <- function(){
  count <- get("objects_on_stack",envir = .pkgenv)
  if(count <1) stop("No more scions on stack. Perhaps you meant to specify the \"data2\" argument to graft()?")
  .data <- get(paste0("stack_obj_",count), envir = .pkgenv)
  assign("objects_on_stack", (count-1),.pkgenv)
  return(.data)
}


#' @title See what's on the stack
#' @description This is primarily to help with debugging.
#' @param x optional string. If supplied it should match the name of an object in the package enviroment.
#' The value of the corresponding variable will be returned. If missing, a list of all objects in the package enviroment.
#' @note Note that \code{\link{graft}} does not delete objects from the environment. See \code{\link{clear_stack}} for this behavior.
#' @export
stack_view<- function(x){
  if(missing(x)) ls(envir = .pkgenv)
  else get(x, envir =.pkgenv)
}
